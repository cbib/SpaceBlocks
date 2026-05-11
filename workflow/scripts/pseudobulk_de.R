################################################################################
## pseudobulk_de.R – DESeq2 Wald + LRT, EnhancedVolcano, ComplexHeatmap,
##                   DEGpatterns
################################################################################

# ── Redirect logs FIRST — before anything else ──────────────────────────────
# This ensures any crash during param parsing is captured in the log files.
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
message("Loading libraries …")

suppressPackageStartupMessages({
  library(DESeq2)
  library(EnhancedVolcano)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(tidyverse)
  library(DEGreport)
})

message("Libraries loaded OK")

# ── Snakemake params (with diagnostics) ──────────────────────────────────────
message("Parsing snakemake params …")

agg_dir        <- snakemake@input[["agg_dir"]]
results_dir    <- snakemake@output[["results_dir"]]
message("  agg_dir:     ", agg_dir)
message("  results_dir: ", results_dir)

min_replicates <- as.integer(snakemake@params[["min_replicates"]])
de_n_genes     <- as.integer(snakemake@params[["de_n_genes"]])
padj_thr       <- as.numeric(snakemake@params[["padj_threshold"]])
lfc_thr        <- as.numeric(snakemake@params[["lfc_threshold"]])
message("  min_reps: ", min_replicates, ", de_n_genes: ", de_n_genes,
        ", padj: ", padj_thr, ", lfc: ", lfc_thr)

# Region levels
message("  Parsing region_levels …")
region_levels <- tryCatch({
  rl <- snakemake@params[["region_levels"]]
  if (is.null(rl)) character(0) else as.character(rl)
}, error = function(e) {
  message("  WARNING: Failed to parse region_levels: ", conditionMessage(e))
  character(0)
})
message("  region_levels: [", paste(region_levels, collapse = ", "), "]")

# Region colors — reconstruct named vector from parallel lists
message("  Parsing region_colors …")
region_colors <- tryCatch({
  rc_names  <- snakemake@params[["region_color_names"]]
  rc_values <- snakemake@params[["region_color_values"]]
  message("    rc_names type: ", class(rc_names), " length: ", length(rc_names))
  message("    rc_values type: ", class(rc_values), " length: ", length(rc_values))
  if (length(rc_names) > 0 && length(rc_names) == length(rc_values)) {
    setNames(as.character(rc_values), as.character(rc_names))
  } else {
    character(0)
  }
}, error = function(e) {
  message("  WARNING: Failed to parse region_colors: ", conditionMessage(e))
  character(0)
})
message("  region_colors: ", paste(names(region_colors), region_colors, sep = "=", collapse = ", "))

message("Params parsed OK")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


# ── ComplexHeatmap helper ────────────────────────────────────────────────────
draw_complex_heatmaps <- function(mat_scaled, meta, condition_col, top_genes,
                                   contrast_name, out_dir, region_lvls, region_cols) {
  # Reorder columns by region level
  if (length(region_lvls) > 0) {
    col_order_factor <- factor(meta[[condition_col]], levels = region_lvls, ordered = TRUE)
  } else {
    col_order_factor <- factor(meta[[condition_col]])
  }
  ord <- order(col_order_factor)
  mat_ordered <- mat_scaled[, ord, drop = FALSE]
  meta_ordered <- meta[ord, , drop = FALSE]
  split_factor <- factor(meta_ordered[[condition_col]],
                         levels = if (length(region_lvls) > 0) region_lvls else unique(meta_ordered[[condition_col]]),
                         ordered = TRUE)

  # Color mapping for annotation
  avail_cols <- region_cols[names(region_cols) %in% levels(split_factor)]
  if (length(avail_cols) > 0) {
    col_anno <- HeatmapAnnotation(
      Region = meta_ordered[[condition_col]],
      col = list(Region = avail_cols),
      show_legend = TRUE
    )
  } else {
    col_anno <- HeatmapAnnotation(
      Region = meta_ordered[[condition_col]],
      show_legend = TRUE
    )
  }

  col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

  # 1. Unsplit heatmap
  tryCatch({
    ht_unsplit <- Heatmap(
      mat_ordered,
      name = "Z-score",
      col = col_fun,
      top_annotation = col_anno,
      cluster_columns = FALSE,
      cluster_rows = TRUE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_gp = gpar(fontsize = 7),
      column_names_gp = gpar(fontsize = 7),
      column_title = paste0("Top DE genes — ", contrast_name),
      heatmap_legend_param = list(title = "Z-score")
    )
    png(file.path(out_dir, paste0("heatmap_unsplit_", contrast_name, ".png")),
        width = max(800, ncol(mat_ordered) * 50),
        height = max(400, nrow(mat_ordered) * 20), res = 150)
    draw(ht_unsplit, merge_legend = TRUE)
    dev.off()
  }, error = function(e) message("      Unsplit heatmap failed: ", conditionMessage(e)))

  # 2. Split heatmap by region level
  tryCatch({
    ht_split <- Heatmap(
      mat_ordered,
      name = "Z-score",
      col = col_fun,
      top_annotation = col_anno,
      column_split = split_factor,
      cluster_columns = FALSE,
      cluster_column_slices = FALSE,
      cluster_rows = TRUE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_gp = gpar(fontsize = 7),
      column_names_gp = gpar(fontsize = 7),
      column_title = paste0("Top DE genes — ", contrast_name, " (split)"),
      heatmap_legend_param = list(title = "Z-score")
    )
    png(file.path(out_dir, paste0("heatmap_split_", contrast_name, ".png")),
        width = max(900, ncol(mat_ordered) * 55),
        height = max(400, nrow(mat_ordered) * 20), res = 150)
    draw(ht_split, merge_legend = TRUE)
    dev.off()
  }, error = function(e) message("      Split heatmap failed: ", conditionMessage(e)))
}


# ── Helper: run DE on one pseudobulk subgroup ────────────────────────────────
run_de_for_prefix <- function(prefix, condition_col, sample_col, out_base) {

  counts_file   <- file.path(agg_dir, "matrices", paste0("counts_", prefix, ".tsv"))
  metadata_file <- file.path(agg_dir, "metadata", paste0("metadata_", prefix, ".tsv"))

  if (!file.exists(counts_file) || !file.exists(metadata_file)) {
    message("  Skipping '", prefix, "': files not found.")
    return(NULL)
  }

  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

  counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
  meta   <- read.delim(metadata_file, row.names = 1, check.names = FALSE)

  counts <- round(counts)
  counts <- as.matrix(counts)
  storage.mode(counts) <- "integer"

  common <- intersect(rownames(counts), rownames(meta))
  if (length(common) == 0) {
    message("  '", prefix, "': no common samples. Skipping.")
    return(NULL)
  }
  counts <- counts[common, , drop = FALSE]
  meta   <- meta[common, , drop = FALSE]

  # Filter to configured region levels if provided
  if (length(region_levels) > 0) {
    keep <- meta[[condition_col]] %in% region_levels
    if (sum(!keep) > 0) {
      dropped <- unique(meta[[condition_col]][!keep])
      message("  Filtering out region levels not in config: ", paste(dropped, collapse = ", "))
    }
    counts <- counts[keep, , drop = FALSE]
    meta   <- meta[keep, , drop = FALSE]
    meta[[condition_col]] <- factor(meta[[condition_col]], levels = region_levels)
  } else {
    meta[[condition_col]] <- factor(meta[[condition_col]])
  }

  meta[[condition_col]] <- droplevels(meta[[condition_col]])
  conditions <- levels(meta[[condition_col]])

  cond_counts <- table(meta[[condition_col]])
  valid_conds <- names(cond_counts[cond_counts >= min_replicates])
  if (length(valid_conds) < 2) {
    message("  '", prefix, "': < 2 conditions with >= ", min_replicates, " replicates. Skipping.")
    writeLines("Insufficient replicates.", file.path(out_base, "SKIPPED_insufficient.txt"))
    return(NULL)
  }
  keep_samples <- meta[[condition_col]] %in% valid_conds
  counts <- counts[keep_samples, , drop = FALSE]
  meta   <- meta[keep_samples, , drop = FALSE]
  meta[[condition_col]] <- droplevels(meta[[condition_col]])
  conditions <- levels(meta[[condition_col]])

  gene_sums <- colSums(counts)
  counts <- counts[, gene_sums > 0, drop = FALSE]

  message("  '", prefix, "': ", nrow(counts), " samples, ", ncol(counts), " genes, ",
          length(conditions), " conditions: ", paste(conditions, collapse = " -> "))

  counts_t <- t(counts)

  if (!identical(colnames(counts_t), rownames(meta))) {
    meta <- meta[colnames(counts_t), , drop = FALSE]
  }

  design_formula <- as.formula(paste("~", condition_col))
  dds <- DESeqDataSetFromMatrix(
    countData = counts_t, colData = meta, design = design_formula
  )

  # ── LRT ────────────────────────────────────────────────────────────────
  message("    LRT …")
  dds_lrt <- NULL
  tryCatch({
    dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)
    res_lrt <- results(dds_lrt)
    res_lrt_df <- as.data.frame(res_lrt) %>%
      rownames_to_column("gene") %>%
      arrange(desc(log2FoldChange))
    write.table(res_lrt_df, file.path(out_base, "LRT_results.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    n_sig <- sum(res_lrt_df$padj < padj_thr, na.rm = TRUE)
    message("    LRT: ", n_sig, " significant genes (padj < ", padj_thr, ")")
  }, error = function(e) {
    message("    LRT failed: ", conditionMessage(e))
  })

  # ── Pairwise Wald ──────────────────────────────────────────────────────
  message("    Wald pairwise tests …")
  dds_wald <- DESeq(dds, test = "Wald")

  vsd <- tryCatch(
    vst(dds_wald, blind = FALSE),
    error = function(e) {
      message("    vst failed, using normTransform: ", conditionMessage(e))
      normTransform(dds_wald)
    }
  )

  pairwise_dir <- file.path(out_base, "pairwise")
  dir.create(pairwise_dir, recursive = TRUE, showWarnings = FALSE)
  summary_rows <- list()
  pairs <- combn(conditions, 2, simplify = FALSE)

  for (pair in pairs) {
    cond_a <- pair[1]
    cond_b <- pair[2]
    contrast_name <- paste0(cond_a, "_vs_", cond_b)
    message("      Contrast: ", contrast_name)

    tryCatch({
      res <- results(dds_wald, contrast = c(condition_col, cond_a, cond_b))
      res_df <- as.data.frame(res) %>%
        rownames_to_column("gene") %>%
        arrange(desc(log2FoldChange))

      write.table(res_df, file.path(pairwise_dir, paste0(contrast_name, ".tsv")),
                  sep = "\t", quote = FALSE, row.names = FALSE)

      n_sig  <- sum(res_df$padj < padj_thr, na.rm = TRUE)
      n_up   <- sum(res_df$padj < padj_thr & res_df$log2FoldChange > 0, na.rm = TRUE)
      n_down <- sum(res_df$padj < padj_thr & res_df$log2FoldChange < 0, na.rm = TRUE)
      message("      ", n_sig, " sig (up:", n_up, ", down:", n_down, ")")

      summary_rows[[contrast_name]] <- tibble(
        contrast = contrast_name, n_tested = nrow(res_df),
        n_significant = n_sig, n_up = n_up, n_down = n_down
      )

      # EnhancedVolcano
      tryCatch({
        p <- EnhancedVolcano(
          res_df, lab = res_df$gene,
          x = "log2FoldChange", y = "padj",
          pCutoff = padj_thr, FCcutoff = lfc_thr,
          title = contrast_name,
          subtitle = paste0(prefix, " — Wald test"),
          legendPosition = "right"
        )
        ggsave(file.path(pairwise_dir, paste0("volcano_", contrast_name, ".png")),
               plot = p, width = 10, height = 8, dpi = 300)
      }, error = function(e) {
        message("      Volcano failed: ", conditionMessage(e))
      })

      # ComplexHeatmap
      tryCatch({
        sig_genes <- res_df %>% filter(padj < padj_thr)
        top_up   <- sig_genes %>% filter(log2FoldChange > 0) %>% head(de_n_genes) %>% pull(gene)
        top_down <- sig_genes %>% filter(log2FoldChange < 0) %>% head(de_n_genes) %>% pull(gene)
        top_genes <- c(top_up, top_down)
        top_genes <- top_genes[top_genes %in% rownames(assay(vsd))]

        if (length(top_genes) >= 2) {
          mat <- assay(vsd)[top_genes, , drop = FALSE]
          mat_scaled <- t(scale(t(mat)))

          draw_complex_heatmaps(
            mat_scaled, meta, condition_col, top_genes,
            contrast_name, pairwise_dir, conditions, region_colors
          )
        }
      }, error = function(e) {
        message("      Heatmap failed: ", conditionMessage(e))
      })

    }, error = function(e) {
      message("      Wald test failed: ", conditionMessage(e))
      summary_rows[[contrast_name]] <- tibble(
        contrast = contrast_name, n_tested = 0,
        n_significant = 0, n_up = 0, n_down = 0, error = conditionMessage(e)
      )
    })
  }

  if (length(summary_rows) > 0) {
    bind_rows(summary_rows) %>%
      write.table(file.path(out_base, "pairwise_summary.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # ── DEGpatterns ────────────────────────────────────────────────────────
  if (!is.null(dds_lrt)) {
    tryCatch({
      res_lrt_all <- results(dds_lrt)
      sig_genes <- rownames(res_lrt_all)[which(res_lrt_all$padj < padj_thr)]

      if (length(sig_genes) >= 10) {
        message("    DEGpatterns on ", length(sig_genes), " LRT-significant genes …")
        sig_mat <- assay(vsd)[sig_genes, , drop = FALSE]

        deg_dir <- file.path(out_base, "DEGpatterns")
        dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)

        clusters <- degPatterns(
          sig_mat, metadata = meta, time = condition_col, plot = FALSE
        )

        if (!is.null(clusters$df)) {
          cluster_df <- clusters$df %>% arrange(cluster, genes)
          write.table(cluster_df, file.path(deg_dir, "gene_clusters.tsv"),
                      sep = "\t", quote = FALSE, row.names = FALSE)

          # Base plot from degPlotCluster
          p <- degPlotCluster(
            clusters$normalized, time = condition_col,
            color = condition_col, points = TRUE
          )

          # Apply custom colors if available
          if (length(region_colors) > 0) {
            p <- p +
              ggplot2::aes(col = .data[[condition_col]]) +
              ggplot2::scale_color_manual(values = region_colors)
          }

          # Add loess trend line and rotate x-axis labels
          p <- p +
            ggplot2::geom_smooth(
              mapping = ggplot2::aes(
                x = .data[[condition_col]],
                y = value,
                group = 1
              ),
              method = "loess",
              color = "black",
              se = FALSE,
              linewidth = 1.2
            ) +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
            )

          ggsave(file.path(deg_dir, "DEGpatterns_clusters.png"),
                 plot = p, width = 14, height = 10, dpi = 500)
          message("    DEGpatterns: ", length(unique(cluster_df$cluster)), " clusters")
        }
      } else {
        message("    Too few LRT-significant genes for DEGpatterns (", length(sig_genes), ")")
      }
    }, error = function(e) {
      message("    DEGpatterns failed: ", conditionMessage(e))
    })
  }

  return(summary_rows)
}


# ── Main ─────────────────────────────────────────────────────────────────────
manifest_file <- file.path(agg_dir, "manifest.tsv")

if (!file.exists(manifest_file)) {
  message("No manifest.tsv found.")
  writeLines("No pseudobulk matrices found.", file.path(results_dir, "SKIPPED_no_manifest.txt"))
  if (exists("log_err_con")) { sink(type = "message"); close(log_err_con) }
  if (exists("log_out_con")) { sink(); close(log_out_con) }
  quit(save = "no", status = 0)
}

manifest <- read.delim(manifest_file)
message("Manifest: ", nrow(manifest), " subgroups")

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, ]
  prefix        <- row$prefix
  condition_col <- row$condition_col
  sample_col    <- row$sample_col

  message("\n--- Processing: ", prefix, " ---")

  out_base <- if (prefix == "pooled") results_dir else file.path(results_dir, prefix)
  tryCatch(
    run_de_for_prefix(prefix, condition_col, sample_col, out_base),
    error = function(e) {
      message("  FAILED for '", prefix, "': ", conditionMessage(e))
      dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
      writeLines(conditionMessage(e), file.path(out_base, "ERROR.txt"))
    }
  )
}

message("\n=== Pseudobulk DE complete ===")

if (exists("log_err_con")) { sink(type = "message"); close(log_err_con) }
if (exists("log_out_con")) { sink(); close(log_out_con) }
