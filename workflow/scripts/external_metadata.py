"""Read cell metadata without changing the user-provided TSV."""
from __future__ import annotations

from collections.abc import Iterable
from os import PathLike

import pandas as pd


_MISSING_LABELS = {"", "na", "n/a", "nan", "none", "null", "unannotated"}


def optional_input_path(value) -> str:
    """Return a path from an optional Snakemake input value."""
    if value is None:
        return ""
    if isinstance(value, (str, PathLike)):
        return str(value)
    return str(value[0]) if len(value) else ""


def read_cell_metadata(
    path: str,
    *,
    id_column: str = "cell_id",
    required_columns: Iterable[str] = (),
) -> pd.DataFrame:
    """Load a TSV and use ``cell_id`` as its index when that column is present.

    All input columns are retained. Empty and duplicate cell identifiers are rejected
    because either case would make downstream joins ambiguous.
    """
    metadata = pd.read_csv(path, sep="\t", index_col=0, comment="#", dtype=str)
    if metadata.empty:
        raise ValueError(f"Cell metadata is empty: {path}")

    missing_columns = [column for column in required_columns if column not in metadata]
    if missing_columns:
        raise ValueError(
            f"Required column(s) {missing_columns} not found in {path}; "
            f"available columns: {list(metadata.columns)}"
        )

    if id_column in metadata.columns:
        cell_ids = metadata[id_column].astype("string").str.strip()
    else:
        cell_ids = pd.Series(metadata.index, index=metadata.index, dtype="string").str.strip()

    missing_ids = cell_ids.isna() | cell_ids.eq("")
    if bool(missing_ids.any()):
        raise ValueError(
            f"Cell metadata {path} contains {int(missing_ids.sum())} empty "
            f"'{id_column}' value(s)."
        )

    duplicate_ids = cell_ids.duplicated(keep=False)
    if bool(duplicate_ids.any()):
        examples = cell_ids[duplicate_ids].drop_duplicates().head(5).tolist()
        raise ValueError(
            f"Cell metadata {path} contains {int(duplicate_ids.sum())} rows with "
            f"duplicate '{id_column}' values (examples: {examples})."
        )

    metadata.index = pd.Index(cell_ids.astype(str), name=id_column)
    return metadata


def normalized_labels(metadata: pd.DataFrame, label_column: str) -> pd.Series:
    """Return trimmed labels, representing blank/unannotated values as missing."""
    labels = metadata[label_column].astype("string").str.strip()
    missing = labels.isna() | labels.str.casefold().isin(_MISSING_LABELS)
    return labels.mask(missing)


def annotated_cell_ids(metadata: pd.DataFrame, label_column: str) -> set[str]:
    """Return IDs whose external annotation is non-empty and not ``Unannotated``."""
    return set(normalized_labels(metadata, label_column).dropna().index.astype(str))
