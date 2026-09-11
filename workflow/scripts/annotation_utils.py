"""Shared helpers for deciding which cell-type annotations are meaningful.

The names below are SpaceBlocks' canonical ``adata.obs`` columns. In particular,
``external_annotation.column`` may name any column in the source metadata TSV;
``annotate_cells.py`` copies that source column to ``cell_type_external`` before
these helpers are called.
"""

from __future__ import annotations


ACTIVE_ANNOTATION_COLUMNS_KEY = "active_annotation_columns"
CANONICAL_ANNOTATION_COLUMNS = (
    "cell_type_tsv",
    "cell_type_ingest",
    "cell_type_external",
    "cell_type_refined",
)
_PLACEHOLDER_LABELS = {
    "", "<na>", "n/a", "na", "nan", "none", "null", "unannotated",
}


def non_placeholder_annotation_mask(labels):
    """Return a boolean mask selecting real labels from a pandas Series."""
    normalized = labels.astype("string").str.strip().str.casefold()
    return labels.notna() & ~normalized.isin(_PLACEHOLDER_LABELS)


def normalize_annotation_labels(labels, placeholder="Unannotated"):
    """Return stripped, object-backed labels with placeholders normalized.

    Pandas' nullable string dtype is useful while normalizing because it preserves
    missing values through string operations.  Do not expose that dtype to callers,
    though: nullable-string categories use a newer h5ad encoding that older AnnData
    readers cannot consume.
    """
    normalized = labels.astype("string").str.strip()
    normalized = normalized.where(
        non_placeholder_annotation_mask(labels), placeholder
    )
    return normalized.astype(object)


def meaningful_annotation_labels(adata, column):
    """Return distinct, non-placeholder labels from an ``obs`` column."""
    if not column or column not in adata.obs.columns:
        return []
    labels = adata.obs[column]
    return (
        labels[non_placeholder_annotation_mask(labels)]
        .astype(str)
        .str.strip()
        .unique()
        .tolist()
    )


def has_meaningful_annotation(adata, column, min_labels=1):
    """Whether ``column`` contains at least ``min_labels`` real labels."""
    return len(meaningful_annotation_labels(adata, column)) >= min_labels


def record_active_annotation_columns(adata, candidates=CANONICAL_ANNOTATION_COLUMNS):
    """Store and return canonical annotation columns carrying at least one real label.

    Placeholder labels are ignored for this availability decision only. Thus, a
    column containing both ``Unannotated`` and real cell types remains active and
    retains its unannotated cells in downstream plots.
    """
    active = [column for column in candidates if has_meaningful_annotation(adata, column)]
    adata.uns[ACTIVE_ANNOTATION_COLUMNS_KEY] = active
    return active


def active_annotation_columns(adata, candidates=CANONICAL_ANNOTATION_COLUMNS):
    """Return usable annotation columns, with a fallback for older h5ad files.

    New outputs record the active columns explicitly. Objects written before that
    metadata existed are inspected semantically, so an all-``Unannotated`` column
    never becomes a report panel merely because it is present in ``obs``.
    """
    candidates = list(candidates)
    if ACTIVE_ANNOTATION_COLUMNS_KEY not in adata.uns:
        return [
            column
            for column in candidates
            if has_meaningful_annotation(adata, column)
        ]

    recorded = adata.uns[ACTIVE_ANNOTATION_COLUMNS_KEY]
    if recorded is None:
        recorded = []
    elif isinstance(recorded, str):
        recorded = [recorded]
    else:
        recorded = [str(column) for column in recorded]
    return [
        column
        for column in candidates
        if column in recorded and has_meaningful_annotation(adata, column)
    ]
