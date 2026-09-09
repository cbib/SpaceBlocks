################################################################################
## pseudobulk_de.R – explicit DESeq2 Wald contrasts and optional omnibus LRT
################################################################################

# Redirect logs before parsing parameters so early failures are captured.
log_out <- tryCatch(snakemake@log[["out"]], error = function(e) "")
log_err <- tryCatch(snakemake@log[["err"]], error = function(e) "")
if (!is.null(log_out) && nzchar(log_out)) {
  dir.create(dirname(log_out), recursive = TRUE, showWarnings = FALSE)
  log_out_con <- file(log_out, open = "wt")
  sink(log_out_con, split = TRUE)
}
if (!is.null(log_err) && nzchar(log_err)) {
  dir.create(dirname(log_err), recursive = TRUE, showWarnings = FALSE)
  log_err_con <- file(log_err, open = "wt")
  sink(log_err_con, type = "message")
}

message("=== pseudobulk_de.R starting ===")
suppressPackageStartupMessages({
  library(DESeq2)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
})

agg_dir <- snakemake@input[["agg_dir"]]
results_dir <- snakemake@output[["results_dir"]]
analysis_name <- as.character(snakemake@params[["analysis_name"]])
min_replicates <- as.integer(snakemake@params[["min_replicates"]])
de_n_genes <- as.integer(snakemake@params[["de_n_genes"]])
padj_thr <- as.numeric(snakemake@params[["padj_threshold"]])
lfc_thr <- as.numeric(snakemake@params[["lfc_threshold"]])
res_dpi <- as.integer(snakemake@params[["dpi"]])
group_by_columns <- as.character(snakemake@params[["group_by_columns"]])
paired_by <- as.character(snakemake@params[["paired_by"]])
covariates <- as.character(snakemake@params[["covariates"]])
contrast_names <- as.character(snakemake@params[["contrast_names"]])
contrast_numerators <- as.character(snakemake@params[["contrast_numerators"]])
contrast_denominators <- as.character(snakemake@params[["contrast_denominators"]])
lrt_enabled <- isTRUE(as.logical(snakemake@params[["lrt_enabled"]]))
paired_by <- paired_by[nzchar(paired_by)]
covariates <- covariates[nzchar(covariates)]
condition_level_order <- tryCatch(
  as.character(snakemake@params[["condition_level_order"]]),
  error = function(e) character(0)
)
condition_level_order <- unique(condition_level_order[nzchar(condition_level_order)])
px_scale <- res_dpi / 150

if (!(length(contrast_names) == length(contrast_numerators) &&
      length(contrast_names) == length(contrast_denominators))) {
  stop("Contrast parameter lists have inconsistent lengths.")
}

# Reconstruct configured palettes from parallel lists.
group_colors <- tryCatch({
  nms <- as.character(snakemake@params[["group_color_names"]])
  vals <- as.character(snakemake@params[["group_color_values"]])
  if (length(nms) > 0 && length(nms) == length(vals)) setNames(vals, nms)
  else character(0)
}, error = function(e) character(0))

extra_annot_columns <- tryCatch(
  as.character(snakemake@params[["extra_annot_columns"]]),
  error = function(e) character(0)
)
extra_annot_columns <- extra_annot_columns[nzchar(extra_annot_columns)]
extra_palettes <- tryCatch({
  cols <- as.character(snakemake@params[["extra_anno_col_names"]])
  vals <- as.character(snakemake@params[["extra_anno_values"]])
  hex <- as.character(snakemake@params[["extra_anno_colors"]])
  palettes <- list()
  if (length(cols) > 0 && length(cols) == length(vals) && length(vals) == length(hex)) {
    for (column in unique(cols)) {
      idx <- which(cols == column)
      palettes[[column]] <- setNames(hex[idx], vals[idx])
    }
  }
  palettes
}, error = function(e) list())

message("  analysis: ", analysis_name)
message("  group_by: [", paste(group_by_columns, collapse = ", "), "]")
message("  paired_by: ", ifelse(length(paired_by), paired_by, "<none>"))
message("  covariates: [", paste(covariates, collapse = ", "), "]")
message("  contrasts: [", paste(contrast_names, collapse = ", "), "]")
message("  LRT: ", lrt_enabled)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


safe_name <- function(value) {
  safe <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(value))
  safe <- gsub("^_+|_+$", "", safe)
  ifelse(nzchar(safe), safe, "contrast")
}


heatmap_condition_factor <- function(values) {
  observed <- unique(as.character(values))
  preferred <- condition_level_order[condition_level_order %in% observed]
  factor(as.character(values), levels = unique(c(preferred, observed)))
}


formula_for <- function(condition_col, paired = character(0), covars = character(0),
                        include_condition = TRUE) {
  terms <- c(covars, paired)
  if (include_condition) terms <- c(terms, condition_col)
  terms <- unique(terms[nzchar(terms)])
  if (length(terms) == 0) return(~1)
  as.formula(paste("~", paste(sprintf("`%s`", terms), collapse = " + ")))
}


assert_full_rank <- function(meta, design_formula, label) {
  matrix <- model.matrix(design_formula, data = meta)
  rank <- qr(matrix)$rank
  if (rank < ncol(matrix)) {
    stop(
      label, " model matrix is not full rank (rank ", rank, "/", ncol(matrix),
      "). Check for confounded covariates or pairing variables."
    )
  }
  invisible(matrix)
}


build_top_annotation <- function(meta, condition_col) {
  annotation_args <- list(Group = as.character(meta[[condition_col]]))
  annotation_colors <- list()
  present_groups <- unique(as.character(meta[[condition_col]]))
  available <- group_colors[names(group_colors) %in% present_groups]
  if (length(available) > 0) annotation_colors[["Group"]] <- available

  for (column in extra_annot_columns) {
    if (column %in% colnames(meta)) {
      values <- as.character(meta[[column]])
      annotation_args[[column]] <- values
      levels <- unique(values)
      palette <- setNames(rep("#cccccc", length(levels)), levels)
      configured <- extra_palettes[[column]]
      if (!is.null(configured)) {
        common <- intersect(names(configured), levels)
        palette[common] <- configured[common]
      }
      annotation_colors[[column]] <- palette
    }
  }
  if (length(annotation_colors) > 0) annotation_args[["col"]] <- annotation_colors
  annotation_args[["show_legend"]] <- TRUE
  do.call(HeatmapAnnotation, annotation_args)
}


draw_complex_heatmaps <- function(mat_scaled, meta, condition_col, contrast_name, out_dir) {
  finite_rows <- apply(mat_scaled, 1, function(values) all(is.finite(values)))
  mat_scaled <- mat_scaled[finite_rows, , drop = FALSE]
  if (nrow(mat_scaled) < 2) return(invisible(NULL))

  condition <- heatmap_condition_factor(meta[[condition_col]])
  order_idx <- order(condition)
  mat_ordered <- mat_scaled[, order_idx, drop = FALSE]
  meta_ordered <- meta[order_idx, , drop = FALSE]
  split_factor <- droplevels(condition[order_idx])
  annotation <- build_top_annotation(meta_ordered, condition_col)
  color_function <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

  n_row <- nrow(mat_ordered)
  n_col <- ncol(mat_ordered)
  body_mm <- max(n_row * 5, 40)
  body_height <- grid::unit(body_mm, "mm")
  device_height <- round((body_mm / 25.4 + 2.6) * res_dpi)

  tryCatch({
    heatmap <- Heatmap(
      mat_ordered, name = "Z-score", col = color_function,
      top_annotation = annotation, height = body_height,
      cluster_columns = TRUE, cluster_rows = TRUE,
      show_row_names = TRUE, show_column_names = TRUE,
      row_names_gp = grid::gpar(fontsize = 7), column_names_gp = grid::gpar(fontsize = 7),
      column_title = paste0("Top DE genes — ", contrast_name)
    )
    png(
      file.path(out_dir, paste0("heatmap_unsplit_", safe_name(contrast_name), ".png")),
      width = round(max(800, n_col * 50) * px_scale),
      height = device_height, res = res_dpi
    )
    draw(heatmap, merge_legend = TRUE)
    dev.off()
  }, error = function(e) message("      Unsplit heatmap failed: ", conditionMessage(e)))

  tryCatch({
    heatmap <- Heatmap(
      mat_ordered, name = "Z-score", col = color_function,
      top_annotation = annotation, column_split = split_factor,
      height = body_height, cluster_columns = FALSE,
      cluster_column_slices = FALSE, cluster_rows = TRUE,
      show_row_names = TRUE, show_column_names = TRUE,
      row_names_gp = grid::gpar(fontsize = 7), column_names_gp = grid::gpar(fontsize = 7),
      column_title = paste0("Top DE genes — ", contrast_name, " (split)")
    )
    png(
      file.path(out_dir, paste0("heatmap_split_", safe_name(contrast_name), ".png")),
      width = round(max(900, n_col * 55) * px_scale),
      height = device_height, res = res_dpi
    )
    draw(heatmap, merge_legend = TRUE)
    dev.off()
  }, error = function(e) message("      Split heatmap failed: ", conditionMessage(e)))
}


write_volcano <- function(results_df, contrast_name, prefix, out_dir) {
  volcano <- results_df %>%
    filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
    mutate(
      sig_group = case_when(
        padj < padj_thr & log2FoldChange > lfc_thr ~ "Up",
        padj < padj_thr & log2FoldChange < -lfc_thr ~ "Down",
        TRUE ~ "NS"
      ),
      neg_log10p = -log10(pmax(padj, 1e-300))
    )
  labels <- bind_rows(
    volcano %>% filter(sig_group == "Up") %>% arrange(desc(log2FoldChange)) %>% head(10),
    volcano %>% filter(sig_group == "Down") %>% arrange(log2FoldChange) %>% head(10)
  )
  plot <- ggplot(volcano, aes(x = log2FoldChange, y = neg_log10p, color = sig_group)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_color_manual(values = c("NS" = "grey70", "Up" = "#e41a1c", "Down" = "#377eb8"), name = "") +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", color = "grey40", linewidth = 0.3) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", color = "grey40", linewidth = 0.3) +
    ggrepel::geom_text_repel(
      data = labels, aes(label = gene), size = 2.5, max.overlaps = 20,
      color = "black", segment.size = 0.2, segment.color = "grey50"
    ) +
    labs(title = contrast_name, subtitle = paste0(prefix, " — Wald test"),
         x = "log2 Fold Change", y = "-log10(padj)") +
    theme_bw()
  ggsave(
    file.path(out_dir, paste0("volcano_", safe_name(contrast_name), ".png")),
    plot = plot, width = 9, height = 7, dpi = res_dpi
  )
}


write_de_heatmap <- function(results_df, vsd, meta, condition_col, contrast_name, out_dir) {
  significant <- results_df %>% filter(padj < padj_thr)
  top_up <- significant %>%
    filter(log2FoldChange > 0) %>% arrange(desc(log2FoldChange)) %>%
    head(de_n_genes) %>% pull(gene)
  top_down <- significant %>%
    filter(log2FoldChange < 0) %>% arrange(log2FoldChange) %>%
    head(de_n_genes) %>% pull(gene)
  genes <- unique(c(top_up, top_down))
  genes <- genes[genes %in% rownames(assay(vsd))]
  if (length(genes) >= 2) {
    scaled <- t(scale(t(assay(vsd)[genes, , drop = FALSE])))
    draw_complex_heatmaps(scaled, meta, condition_col, contrast_name, out_dir)
  }
}


fit_wald <- function(counts, meta, condition_col, paired = character(0)) {
  design_formula <- formula_for(condition_col, paired, covariates, TRUE)
  assert_full_rank(meta, design_formula, "Wald")
  counts_t <- t(counts)
  counts_t <- counts_t[, rownames(meta), drop = FALSE]
  dds <- DESeqDataSetFromMatrix(countData = counts_t, colData = meta, design = design_formula)
  dds <- DESeq(dds, test = "Wald")
  vsd <- tryCatch(
    vst(dds, blind = FALSE),
    error = function(e) {
      message("    vst failed, using normTransform: ", conditionMessage(e))
      normTransform(dds)
    }
  )
  list(dds = dds, vsd = vsd)
}


run_lrt <- function(counts, meta, condition_col, sample_col, prefix, out_base) {
  if (!lrt_enabled) return(invisible(NULL))
  message("    Explicit omnibus LRT …")

  unit_col <- if (length(paired_by)) paired_by else sample_col
  replicate_counts <- meta %>%
    mutate(.unit = as.character(.data[[unit_col]])) %>%
    group_by(.data[[condition_col]]) %>%
    summarise(n = n_distinct(.unit), .groups = "drop")
  valid <- as.character(replicate_counts[[condition_col]][replicate_counts$n >= min_replicates])
  lrt_meta <- meta[as.character(meta[[condition_col]]) %in% valid, , drop = FALSE]
  lrt_counts <- counts[rownames(lrt_meta), , drop = FALSE]
  lrt_meta[[condition_col]] <- droplevels(factor(lrt_meta[[condition_col]]))

  if (length(paired_by)) {
    informative <- lrt_meta %>%
      mutate(.row = rownames(lrt_meta)) %>%
      group_by(.data[[paired_by]]) %>%
      filter(n_distinct(.data[[condition_col]]) >= 2) %>%
      pull(.row)
    lrt_meta <- lrt_meta[informative, , drop = FALSE]
    lrt_counts <- lrt_counts[informative, , drop = FALSE]
    lrt_meta[[paired_by]] <- droplevels(factor(lrt_meta[[paired_by]]))
  }

  lrt_meta[[condition_col]] <- droplevels(factor(lrt_meta[[condition_col]]))
  if (nlevels(lrt_meta[[condition_col]]) < 2) {
    message("    LRT skipped: fewer than two non-excluded levels have sufficient replicates.")
    lrt_dir <- file.path(out_base, "LRT")
    dir.create(lrt_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines("Insufficient replicated levels for LRT.", file.path(lrt_dir, "SKIPPED_insufficient.txt"))
    return(invisible(NULL))
  }

  full_formula <- formula_for(condition_col, paired_by, covariates, TRUE)
  reduced_formula <- formula_for(condition_col, paired_by, covariates, FALSE)
  assert_full_rank(lrt_meta, full_formula, "LRT full")
  assert_full_rank(lrt_meta, reduced_formula, "LRT reduced")
  counts_t <- t(lrt_counts)
  counts_t <- counts_t[, rownames(lrt_meta), drop = FALSE]
  dds <- DESeqDataSetFromMatrix(countData = counts_t, colData = lrt_meta, design = full_formula)
  dds <- DESeq(dds, test = "LRT", reduced = reduced_formula)
  result <- results(dds)
  # An LRT p-value tests the group term as a whole. DESeq2 also returns one
  # coefficient's fold change, which is not an omnibus effect size and is easy to
  # misinterpret when there are more than two levels, so omit it here.
  result_df <- as.data.frame(result) %>%
    rownames_to_column("gene") %>%
    select(gene, baseMean, stat, pvalue, padj) %>%
    arrange(padj)
  lrt_dir <- file.path(out_base, "LRT")
  dir.create(lrt_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(result_df, file.path(lrt_dir, "LRT_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  significant <- rownames(result)[which(result$padj < padj_thr)]
  message("    LRT: ", length(significant), " significant genes")
  if (length(significant) < 10) return(invisible(NULL))

  tryCatch({
    transformed <- tryCatch(vst(dds, blind = FALSE), error = function(e) normTransform(dds))
    matrix <- assay(transformed)[significant, , drop = FALSE]
    clusters <- DEGreport::degPatterns(
      matrix, metadata = lrt_meta, time = condition_col, plot = FALSE
    )
    if (!is.null(clusters$df)) {
      normalized <- clusters$normalized
      old_levels <- sort(unique(normalized$cluster))
      remap <- setNames(seq_along(old_levels), as.character(old_levels))
      normalized$cluster <- as.integer(remap[as.character(normalized$cluster)])
      cluster_df <- clusters$df
      cluster_df$cluster <- as.integer(remap[as.character(cluster_df$cluster)])
      cluster_df <- cluster_df %>% arrange(cluster, genes)

      plot <- DEGreport::degPlotCluster(
        normalized, time = condition_col, color = condition_col, points = TRUE
      ) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
      if (length(group_colors) > 0) {
        plot <- plot + scale_color_manual(values = group_colors)
      }
      ggsave(file.path(lrt_dir, "DEGpatterns_groups.png"), plot = plot,
             width = 14, height = 10, dpi = max(res_dpi, 500L))

      grouped_genes <- split(cluster_df$genes, paste0("group_", cluster_df$cluster))
      max_length <- max(lengths(grouped_genes))
      gene_table <- as.data.frame(lapply(grouped_genes, function(genes) {
        c(genes, rep("", max_length - length(genes)))
      }))
      write.table(gene_table, file.path(lrt_dir, "gene_groups.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }, error = function(e) message("    DEGpatterns failed: ", conditionMessage(e)))
}


run_de_for_prefix <- function(prefix, condition_col, sample_col, out_base) {
  counts_file <- file.path(agg_dir, "matrices", paste0("counts_", prefix, ".tsv"))
  metadata_file <- file.path(agg_dir, "metadata", paste0("metadata_", prefix, ".tsv"))
  if (!file.exists(counts_file) || !file.exists(metadata_file)) return(NULL)
  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

  counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
  meta <- read.delim(metadata_file, row.names = 1, check.names = FALSE)
  counts <- round(as.matrix(counts))
  storage.mode(counts) <- "integer"
  common <- intersect(rownames(counts), rownames(meta))
  counts <- counts[common, , drop = FALSE]
  meta <- meta[common, , drop = FALSE]

  required <- unique(c(condition_col, sample_col, paired_by, covariates))
  missing <- setdiff(required, colnames(meta))
  if (length(missing) > 0) stop("Missing pseudobulk metadata columns: ", paste(missing, collapse = ", "))
  complete <- complete.cases(meta[, required, drop = FALSE])
  counts <- counts[complete, , drop = FALSE]
  meta <- meta[complete, , drop = FALSE]
  meta[[condition_col]] <- factor(meta[[condition_col]])
  for (column in unique(c(paired_by, covariates))) meta[[column]] <- factor(meta[[column]])

  unit_col <- if (length(paired_by)) paired_by else sample_col
  unit_groups <- meta[, c(unit_col, condition_col), drop = FALSE]
  duplicate_units <- duplicated(unit_groups) | duplicated(unit_groups, fromLast = TRUE)
  if (any(duplicate_units)) {
    examples <- unique(apply(unit_groups[duplicate_units, , drop = FALSE], 1, paste, collapse = " / "))
    stop(
      "More than one pseudobulk observation exists for an independent unit/group: ",
      paste(head(examples, 5), collapse = ", "),
      ". Choose an aggregation aligned with group_by; repeated observations are not independent."
    )
  }
  if (!length(paired_by)) {
    groups_per_sample <- tapply(
      as.character(meta[[condition_col]]), as.character(meta[[sample_col]]),
      function(values) length(unique(values))
    )
    if (any(groups_per_sample > 1)) {
      stop(
        "At least one sample contributes multiple comparison levels. Configure paired_by: ",
        sample_col, " for a paired design."
      )
    }
  }

  counts <- counts[, colSums(counts) > 0, drop = FALSE]
  if (nrow(counts) == 0 || ncol(counts) == 0) {
    writeLines("No usable counts.", file.path(out_base, "SKIPPED_no_counts.txt"))
    return(NULL)
  }

  message("  '", prefix, "': ", nrow(counts), " observations, ", ncol(counts), " genes")
  tryCatch(
    run_lrt(counts, meta, condition_col, sample_col, prefix, out_base),
    error = function(e) {
      message("    LRT failed: ", conditionMessage(e))
      lrt_dir <- file.path(out_base, "LRT")
      dir.create(lrt_dir, recursive = TRUE, showWarnings = FALSE)
      writeLines(conditionMessage(e), file.path(lrt_dir, "ERROR.txt"))
    }
  )

  pairwise_dir <- file.path(out_base, "pairwise")
  dir.create(pairwise_dir, recursive = TRUE, showWarnings = FALSE)
  summary_rows <- list()

  # Unpaired analyses share a single fit across all sufficiently replicated levels.
  shared_fit <- NULL
  shared_meta <- NULL
  if (!length(paired_by) && length(contrast_names) > 0) {
    replicate_counts <- meta %>%
      mutate(.sample = as.character(.data[[sample_col]])) %>%
      group_by(.data[[condition_col]]) %>%
      summarise(n = n_distinct(.sample), .groups = "drop")
    valid <- as.character(replicate_counts[[condition_col]][replicate_counts$n >= min_replicates])
    shared_meta <- meta[as.character(meta[[condition_col]]) %in% valid, , drop = FALSE]
    shared_meta[[condition_col]] <- droplevels(factor(shared_meta[[condition_col]]))
    if (nlevels(shared_meta[[condition_col]]) >= 2) {
      shared_counts <- counts[rownames(shared_meta), , drop = FALSE]
      shared_fit <- fit_wald(shared_counts, shared_meta, condition_col)
    }
  }

  for (index in seq_along(contrast_names)) {
    contrast_name <- contrast_names[index]
    numerator <- contrast_numerators[index]
    denominator <- contrast_denominators[index]
    message("    Contrast: ", contrast_name, " (", numerator, " vs ", denominator, ")")

    tryCatch({
      if (length(paired_by)) {
        pair_meta <- meta[as.character(meta[[condition_col]]) %in% c(numerator, denominator), , drop = FALSE]
        pair_meta[[condition_col]] <- droplevels(factor(pair_meta[[condition_col]]))
        complete_units <- pair_meta %>%
          mutate(.row = rownames(pair_meta)) %>%
          group_by(.data[[paired_by]]) %>%
          filter(all(c(numerator, denominator) %in% as.character(.data[[condition_col]]))) %>%
          pull(.row)
        pair_meta <- pair_meta[complete_units, , drop = FALSE]
        pair_meta[[condition_col]] <- droplevels(factor(pair_meta[[condition_col]]))
        pair_meta[[paired_by]] <- droplevels(factor(pair_meta[[paired_by]]))
        n_replicates <- n_distinct(pair_meta[[paired_by]])
        if (n_replicates < min_replicates) {
          stop("only ", n_replicates, " complete pairs; need ", min_replicates)
        }
        pair_counts <- counts[rownames(pair_meta), , drop = FALSE]
        fit <- fit_wald(pair_counts, pair_meta, condition_col, paired_by)
        fit_meta <- pair_meta
      } else {
        if (is.null(shared_fit) ||
            !all(c(numerator, denominator) %in% levels(shared_meta[[condition_col]]))) {
          stop("one or both levels have fewer than ", min_replicates, " independent samples")
        }
        fit <- shared_fit
        fit_meta <- shared_meta
      }

      result <- results(
        fit$dds,
        contrast = c(condition_col, numerator, denominator),
        alpha = padj_thr
      )
      result_df <- as.data.frame(result) %>%
        rownames_to_column("gene") %>% arrange(desc(log2FoldChange))
      output_name <- safe_name(contrast_name)
      write.table(result_df, file.path(pairwise_dir, paste0(output_name, ".tsv")),
                  sep = "\t", quote = FALSE, row.names = FALSE)

      n_sig <- sum(result_df$padj < padj_thr, na.rm = TRUE)
      n_up <- sum(result_df$padj < padj_thr & result_df$log2FoldChange > 0, na.rm = TRUE)
      n_down <- sum(result_df$padj < padj_thr & result_df$log2FoldChange < 0, na.rm = TRUE)
      summary_rows[[contrast_name]] <- tibble(
        contrast = contrast_name, numerator = numerator, denominator = denominator,
        n_tested = nrow(result_df), n_significant = n_sig, n_up = n_up, n_down = n_down,
        error = NA_character_
      )
      write_volcano(result_df, contrast_name, prefix, pairwise_dir)
      write_de_heatmap(result_df, fit$vsd, fit_meta, condition_col, contrast_name, pairwise_dir)
    }, error = function(e) {
      message("      Contrast skipped: ", conditionMessage(e))
      summary_rows[[contrast_name]] <<- tibble(
        contrast = contrast_name, numerator = numerator, denominator = denominator,
        n_tested = 0L, n_significant = 0L, n_up = 0L, n_down = 0L,
        error = conditionMessage(e)
      )
    })
  }

  if (length(summary_rows) > 0) {
    bind_rows(summary_rows) %>%
      write.table(file.path(out_base, "pairwise_summary.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
  }
  invisible(summary_rows)
}


manifest_file <- file.path(agg_dir, "manifest.tsv")
if (!file.exists(manifest_file)) {
  writeLines("No pseudobulk matrices found.", file.path(results_dir, "SKIPPED_no_manifest.txt"))
} else {
  manifest <- read.delim(manifest_file)
  message("Manifest: ", nrow(manifest), " subgroup(s)")
  for (index in seq_len(nrow(manifest))) {
    row <- manifest[index, ]
    prefix <- row$prefix
    out_base <- if (prefix == "pooled") results_dir else file.path(results_dir, prefix)
    tryCatch(
      run_de_for_prefix(prefix, row$condition_col, row$sample_col, out_base),
      error = function(e) {
        message("  FAILED for '", prefix, "': ", conditionMessage(e))
        dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
        writeLines(conditionMessage(e), file.path(out_base, "ERROR.txt"))
      }
    )
  }
}

message("=== Pseudobulk DE complete ===")
if (exists("log_err_con")) { sink(type = "message"); close(log_err_con) }
if (exists("log_out_con")) { sink(); close(log_out_con) }
