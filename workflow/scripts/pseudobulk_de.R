################################################################################
## pseudobulk_de.R – DESeq2 Wald + LRT, ComplexHeatmap, DEGpatterns,
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
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(tidyverse)
  library(ggrepel)
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
# Plot resolution (config param); defaults to 300 if the rule doesn't pass it.
.dpi_param     <- snakemake@params[["dpi"]]
res_dpi        <- if (is.null(.dpi_param)) 300L else as.integer(.dpi_param)
px_scale       <- res_dpi / 150                  # keep physical size, scale pixels
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

# ── Design annotation columns + palettes (parallel lists → named vectors) ─────
extra_annot_columns <- tryCatch(as.character(snakemake@params[["extra_annot_columns"]]),
                           error = function(e) character(0))
extra_annot_columns <- extra_annot_columns[nzchar(extra_annot_columns)]
extra_palettes <- tryCatch({
  dc <- as.character(snakemake@params[["extra_anno_col_names"]])
  dv <- as.character(snakemake@params[["extra_anno_values"]])
  dh <- as.character(snakemake@params[["extra_anno_colors"]])
  pals <- list()
  if (length(dc) > 0 && length(dc) == length(dv) && length(dv) == length(dh)) {
    for (u in unique(dc)) {
      idx <- which(dc == u)
      pals[[u]] <- setNames(dh[idx], dv[idx])
    }
  }
  pals
}, error = function(e) {
  message("  WARNING: Failed to parse design palettes: ", conditionMessage(e))
  list()
})
message("  extra_annot_columns: [", paste(extra_annot_columns, collapse = ", "), "]")

# Build a ComplexHeatmap top annotation: Region (existing behaviour) plus one
# track per design column present in `meta_df`, coloured from extra_palettes with
# a grey (#cccccc) fallback for values absent from the config palette. Design
# columns are annotation-only — they are NOT added to the DESeq2 model.
build_top_annotation <- function(meta_df, condition_col, region_cols) {
  anno_args <- list(Region = meta_df[[condition_col]])
  anno_cols <- list()
  avail <- region_cols[names(region_cols) %in% unique(as.character(meta_df[[condition_col]]))]
  if (length(avail) > 0) anno_cols[["Region"]] <- avail
  for (dc in extra_annot_columns) {
    if (dc %in% colnames(meta_df)) {
      vals <- as.character(meta_df[[dc]])
      anno_args[[dc]] <- vals
      lv   <- unique(vals)
      full <- setNames(rep("#cccccc", length(lv)), lv)
      pal  <- extra_palettes[[dc]]
      if (!is.null(pal)) {
        common <- intersect(names(pal), lv)
        if (length(common) > 0) full[common] <- pal[common]
      }
      anno_cols[[dc]] <- full
    }
  }
  if (length(anno_cols) > 0) anno_args[["col"]] <- anno_cols
  anno_args[["show_legend"]] <- TRUE
  do.call(HeatmapAnnotation, anno_args)
}

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

  # Drop genes that became all-NaN after scaling (zero variance) — they
  # collapse the layout and break row clustering.
  finite_rows <- apply(mat_ordered, 1, function(z) all(is.finite(z)))
  if (any(!finite_rows)) {
    message("      Dropping ", sum(!finite_rows),
            " zero-variance gene(s) from heatmap")
    mat_ordered <- mat_ordered[finite_rows, , drop = FALSE]
  }
  if (nrow(mat_ordered) < 2) {
    message("      <2 finite genes left — skipping heatmap")
    return(invisible(NULL))
  }

  # Pin the heatmap BODY so each gene row has a fixed physical height.
  n_row <- nrow(mat_ordered); n_col <- ncol(mat_ordered)
  row_mm   <- 5
  body_mm  <- max(n_row * row_mm, 40)
  body_h   <- unit(body_mm, "mm")
  dev_h_px <- round((body_mm / 25.4 + 2.6) * res_dpi)  # +2.6in for title/anno/legend

  # Region + design-column annotation tracks (grey fallback for missing values)
  col_anno <- build_top_annotation(meta_ordered, condition_col, region_cols)

  col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

  # 1. Unsplit heatmap
  tryCatch({
    ht_unsplit <- Heatmap(
      mat_ordered,
      name = "Z-score",
      col = col_fun,
      top_annotation = col_anno,
      height = body_h,
      cluster_columns = TRUE,
      cluster_rows = TRUE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_gp = gpar(fontsize = 7),
      column_names_gp = gpar(fontsize = 7),
      column_title = paste0("Top DE genes — ", contrast_name),
      heatmap_legend_param = list(title = "Z-score")
    )
    png(file.path(out_dir, paste0("heatmap_unsplit_", contrast_name, ".png")),
        width = round(max(800, n_col * 50) * px_scale),
        height = dev_h_px, res = res_dpi)
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
      height = body_h,
      cluster_columns = FALSE,
      cluster_column_slices = TRUE,
      cluster_rows = TRUE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_gp = gpar(fontsize = 7),
      column_names_gp = gpar(fontsize = 7),
      column_title = paste0("Top DE genes — ", contrast_name, " (split)"),
      heatmap_legend_param = list(title = "Z-score")
    )
    png(file.path(out_dir, paste0("heatmap_split_", contrast_name, ".png")),
        width = round(max(900, n_col * 55) * px_scale),
        height = dev_h_px, res = res_dpi)
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
    # LRT results saved inside DEGpatterns dir (created later)
    lrt_save_dir <- file.path(out_base, "DEGpatterns")
    dir.create(lrt_save_dir, recursive = TRUE, showWarnings = FALSE)
    write.table(res_lrt_df, file.path(lrt_save_dir, "LRT_results.tsv"),
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

      # Custom volcano plot (grey=NS, red=up, blue=down, top genes labelled)
      tryCatch({
        volc_df <- res_df %>%
          filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
          mutate(
            sig_group = case_when(
              padj >= padj_thr ~ "NS",
              log2FoldChange > lfc_thr ~ "Up",
              log2FoldChange < -lfc_thr ~ "Down",
              TRUE ~ "NS"
            ),
            neg_log10p = -log10(pmax(padj, 1e-300))
          )

        # Top genes to label
        top_up <- volc_df %>% filter(sig_group == "Up") %>%
          arrange(desc(log2FoldChange)) %>% head(10)
        top_down <- volc_df %>% filter(sig_group == "Down") %>%
          arrange(log2FoldChange) %>% head(10)
        label_genes <- bind_rows(top_up, top_down)

        p <- ggplot(volc_df, aes(x = log2FoldChange, y = neg_log10p, color = sig_group)) +
          geom_point(size = 0.8, alpha = 0.6) +
          scale_color_manual(values = c("NS" = "grey70", "Up" = "#e41a1c", "Down" = "#377eb8"),
                             name = "") +
          geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", color = "grey40", linewidth = 0.3) +
          geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", color = "grey40", linewidth = 0.3) +
          ggrepel::geom_text_repel(
            data = label_genes, aes(label = gene),
            size = 2.5, max.overlaps = 20, color = "black",
            segment.size = 0.2, segment.color = "grey50"
          ) +
          labs(title = contrast_name, subtitle = paste0(prefix, " — Wald test"),
               x = "log2 Fold Change", y = "-log10(padj)") +
          theme_bw() +
          theme(legend.position = "right")

        ggsave(file.path(pairwise_dir, paste0("volcano_", contrast_name, ".png")),
               plot = p, width = 9, height = 7, dpi = res_dpi)
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
          # ── Remap cluster ids to a contiguous 1..N (degPatterns skips ids)
          norm <- clusters$normalized
          old_levels <- sort(unique(norm$cluster))
          remap <- setNames(seq_along(old_levels), as.character(old_levels))
          # keep cluster as INTEGER — degPlotCluster joins on it internally and a
          # factor here triggers "Can't join <factor> with <integer>".
          norm$cluster <- as.integer(remap[as.character(norm$cluster)])

          cluster_df <- clusters$df
          cluster_df$cluster <- as.integer(remap[as.character(cluster_df$cluster)])
          cluster_df <- cluster_df %>% arrange(cluster, genes)

          # Base plot from degPlotCluster (now facetted by 1..N)
          p <- degPlotCluster(
            norm, time = condition_col, color = condition_col, points = TRUE
          )

          # Apply custom colors if available
          if (length(region_colors) > 0) {
            p <- p +
              ggplot2::aes(col = .data[[condition_col]]) +
              ggplot2::scale_color_manual(values = region_colors)
          }

          # Loess trend + rotated x labels
          p <- p +
            ggplot2::geom_smooth(
              mapping = ggplot2::aes(x = .data[[condition_col]], y = value, group = 1),
              method = "loess", color = "black", se = FALSE, linewidth = 1.2
            ) +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
            )

          ggsave(file.path(deg_dir, "DEGpatterns_groups.png"),   # renamed
                 plot = p, width = 14, height = 10, dpi = max(res_dpi, 500L))
          message("    DEGpatterns: ", length(old_levels), " groups (renumbered 1..",
                  length(old_levels), ")")

          # Gene groups TSV — use the remapped, contiguous group numbers
          tryCatch({
            grouped_genes <- split(cluster_df$genes,
                                   paste0("group_", cluster_df$cluster))
            max_len <- max(sapply(grouped_genes, length))
            gene_df <- as.data.frame(do.call(cbind, lapply(grouped_genes, function(x) {
              c(x, rep("", max_len - length(x)))
            })))
            colnames(gene_df) <- gsub("[^A-Za-z0-9_]", "_", colnames(gene_df))
            colnames(gene_df) <- gsub("_+", "_", colnames(gene_df))
            colnames(gene_df) <- gsub("_$", "", colnames(gene_df))
            write.table(gene_df, file.path(deg_dir, "gene_groups.tsv"),
                        sep = "\t", quote = FALSE, row.names = FALSE)
            message("    Gene groups saved: ", ncol(gene_df), " groups")
          }, error = function(e) {
            message("    Gene groups export failed: ", conditionMessage(e))
          })

          # ── Per-group heatmaps: the genes of every DEGpatterns group ─────
          # Reuses the DE heatmap helper (unsplit + region-split, z-scored vsd),
          # written alongside the group plots in DEGpatterns/.
          tryCatch({
            for (g in sort(unique(cluster_df$cluster))) {
              genes_g <- cluster_df$genes[cluster_df$cluster == g]
              genes_g <- genes_g[genes_g %in% rownames(assay(vsd))]
              if (length(genes_g) >= 2) {
                mat_g_scaled <- t(scale(t(assay(vsd)[genes_g, , drop = FALSE])))
                draw_complex_heatmaps(
                  mat_g_scaled, meta, condition_col, genes_g,
                  paste0("group_", g), deg_dir, conditions, region_colors
                )
              }
            }
            message("    DEGpatterns per-group heatmaps written")
          }, error = function(e) {
            message("    DEGpatterns per-group heatmaps failed: ", conditionMessage(e))
          })
        }
      } else {
        message("    Too few LRT-significant genes for DEGpatterns (", length(sig_genes), ")")
      }
    }, error = function(e) {
      message("    DEGpatterns failed: ", conditionMessage(e))
    })
  }

  # ── One-vs-rest comparisons + unique marker identification ──────────
  message("    One-vs-rest comparisons …")
  ovr_dir <- file.path(out_base, "one_vs_rest")
  dir.create(ovr_dir, recursive = TRUE, showWarnings = FALSE)

  ovr_results <- list()  # store results for cross-comparison

  for (focal in conditions) {
    focal_safe <- gsub("[^A-Za-z0-9_]", "_", focal)
    message("      ", focal, " vs rest")

    tryCatch({
      # Create binary factor: focal vs rest
      meta_ovr <- meta
      meta_ovr[["ovr_group"]] <- ifelse(meta_ovr[[condition_col]] == focal, focal, "rest")
      meta_ovr[["ovr_group"]] <- factor(meta_ovr[["ovr_group"]], levels = c("rest", focal))

      dds_ovr <- DESeqDataSetFromMatrix(
        countData = counts_t, colData = meta_ovr, design = ~ ovr_group
      )
      dds_ovr <- DESeq(dds_ovr, test = "Wald")
      res_ovr <- results(dds_ovr, contrast = c("ovr_group", focal, "rest"))

      res_ovr_df <- as.data.frame(res_ovr) %>%
        rownames_to_column("gene") %>%
        arrange(desc(log2FoldChange))

      write.table(res_ovr_df, file.path(ovr_dir, paste0(focal_safe, "_vs_rest.tsv")),
                  sep = "\t", quote = FALSE, row.names = FALSE)

      n_sig <- sum(res_ovr_df$padj < padj_thr, na.rm = TRUE)
      n_up  <- sum(res_ovr_df$padj < padj_thr & res_ovr_df$log2FoldChange > 0, na.rm = TRUE)
      message("      ", n_sig, " sig (", n_up, " up)")

      # Store for cross-comparison
      ovr_results[[focal]] <- res_ovr_df

      # Volcano plot (same style as pairwise)
      tryCatch({
        volc_df <- res_ovr_df %>%
          filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
          mutate(
            sig_group = case_when(
              padj >= padj_thr ~ "NS",
              log2FoldChange > lfc_thr ~ "Up",
              log2FoldChange < -lfc_thr ~ "Down",
              TRUE ~ "NS"
            ),
            neg_log10p = -log10(pmax(padj, 1e-300))
          )
        top_up <- volc_df %>% filter(sig_group == "Up") %>% arrange(desc(log2FoldChange)) %>% head(10)
        top_down <- volc_df %>% filter(sig_group == "Down") %>% arrange(log2FoldChange) %>% head(10)
        label_genes <- bind_rows(top_up, top_down)

        p <- ggplot(volc_df, aes(x = log2FoldChange, y = neg_log10p, color = sig_group)) +
          geom_point(size = 0.8, alpha = 0.6) +
          scale_color_manual(values = c("NS" = "grey70", "Up" = "#e41a1c", "Down" = "#377eb8"), name = "") +
          geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", color = "grey40", linewidth = 0.3) +
          geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", color = "grey40", linewidth = 0.3) +
          ggrepel::geom_text_repel(data = label_genes, aes(label = gene), size = 2.5,
                                    max.overlaps = 20, color = "black",
                                    segment.size = 0.2, segment.color = "grey50") +
          labs(title = paste0(focal, " vs rest"), x = "log2 Fold Change", y = "-log10(padj)") +
          theme_bw() + theme(legend.position = "right")

        ggsave(file.path(ovr_dir, paste0("volcano_", focal_safe, "_vs_rest.png")),
               plot = p, width = 9, height = 7, dpi = res_dpi)
      }, error = function(e) {
        message("      OVR volcano failed: ", conditionMessage(e))
      })

      # Top DE gene heatmap
      tryCatch({
        sig_genes_ovr <- res_ovr_df %>% filter(padj < padj_thr)
        top_up_genes   <- sig_genes_ovr %>% filter(log2FoldChange > 0) %>% head(de_n_genes) %>% pull(gene)
        top_down_genes <- sig_genes_ovr %>% filter(log2FoldChange < 0) %>% head(de_n_genes) %>% pull(gene)
        top_genes_ovr <- c(top_up_genes, top_down_genes)
        top_genes_ovr <- top_genes_ovr[top_genes_ovr %in% rownames(assay(vsd))]

        if (length(top_genes_ovr) >= 2) {
          draw_complex_heatmaps(
            t(scale(t(assay(vsd)[top_genes_ovr, , drop = FALSE]))),
            meta, condition_col, top_genes_ovr,
            paste0(focal_safe, "_vs_rest"), ovr_dir, conditions, region_colors
          )
        }
      }, error = function(e) {
        message("      OVR heatmap failed: ", conditionMessage(e))
      })

    }, error = function(e) {
      message("      OVR test failed for '", focal, "': ", conditionMessage(e))
    })
  }

  # ── Unique marker identification ───────────────────────────────────
  # A gene is a unique marker for level X if it is:
  # - significantly UP in X vs rest (padj < thr, log2FC > lfc_thr)
  # - NOT significantly UP in any other level vs rest
  if (length(ovr_results) >= 2) {
    message("    Identifying unique markers …")
    tryCatch({
      # For each level, get the set of significantly upregulated genes
      up_gene_sets <- lapply(ovr_results, function(df) {
        df %>% filter(padj < padj_thr & log2FoldChange > lfc_thr) %>% pull(gene)
      })

      unique_markers <- list()
      for (lvl in names(up_gene_sets)) {
        others <- setdiff(names(up_gene_sets), lvl)
        other_genes <- unique(unlist(up_gene_sets[others]))
        unique_to_lvl <- setdiff(up_gene_sets[[lvl]], other_genes)

        if (length(unique_to_lvl) > 0) {
          # Get the DE stats for these genes
          lvl_df <- ovr_results[[lvl]] %>%
            filter(gene %in% unique_to_lvl) %>%
            arrange(desc(log2FoldChange))
          unique_markers[[lvl]] <- lvl_df
          message("      ", lvl, ": ", nrow(lvl_df), " unique markers")
        } else {
          message("      ", lvl, ": 0 unique markers")
        }
      }

      # Save as combined TSV
      if (length(unique_markers) > 0) {
        all_markers <- bind_rows(unique_markers, .id = "level") %>%
          arrange(level, desc(log2FoldChange))
        write.table(all_markers, file.path(ovr_dir, "unique_markers.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)

        # Also save as column-format (one column per level, padded)
        max_len <- max(sapply(unique_markers, nrow))
        marker_cols <- lapply(unique_markers, function(df) {
          c(df$gene, rep("", max_len - nrow(df)))
        })
        marker_df <- as.data.frame(marker_cols)
        colnames(marker_df) <- gsub("[^A-Za-z0-9_]", "_", colnames(marker_df))
        write.table(marker_df, file.path(ovr_dir, "unique_markers_by_level.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
      }

      # Summary barplot: number of unique markers per level
      if (length(unique_markers) > 0) {
        counts_bar <- tibble(
          level = factor(names(unique_markers), levels = conditions),
          n_markers = sapply(unique_markers, nrow)
        )
        p <- ggplot(counts_bar, aes(x = level, y = n_markers, fill = level)) +
          geom_col(show.legend = FALSE) +
          geom_text(aes(label = n_markers), vjust = -0.3, size = 3) +
          labs(title = "Unique upregulated markers per level",
               x = "", y = "Number of unique markers") +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        ggsave(file.path(ovr_dir, "unique_markers_barplot.png"),
               plot = p, width = max(6, length(unique_markers) * 1.2), height = 5, dpi = res_dpi)
      }

      # Heatmap: ALL unique marker genes (rows), samples (cols), row-split by
      # the level each gene is unique to. May be tall — that is intended.
      if (length(unique_markers) > 0) {
        tryCatch({
          # gene -> level map, preserving level order in `conditions`
          gene_level <- unlist(lapply(names(unique_markers), function(lvl)
            setNames(rep(lvl, nrow(unique_markers[[lvl]])),
                     unique_markers[[lvl]]$gene)))
          genes_all <- names(gene_level)
          genes_all <- genes_all[genes_all %in% rownames(assay(vsd))]

          if (length(genes_all) >= 2) {
            mat <- assay(vsd)[genes_all, , drop = FALSE]
            mat_scaled <- t(scale(t(mat)))

            # order columns by region level
            col_ord <- order(factor(meta[[condition_col]],
                                    levels = if (length(conditions) > 0) conditions
                                             else unique(meta[[condition_col]])))
            mat_scaled <- mat_scaled[, col_ord, drop = FALSE]
            meta_ord   <- meta[col_ord, , drop = FALSE]

            # drop zero-variance genes (all-NaN rows after scaling)
            finite_rows <- apply(mat_scaled, 1, function(z) all(is.finite(z)))
            mat_scaled  <- mat_scaled[finite_rows, , drop = FALSE]

            if (nrow(mat_scaled) >= 2) {
              row_split <- factor(gene_level[rownames(mat_scaled)],
                                  levels = if (length(conditions) > 0) conditions
                                           else unique(gene_level))

              col_anno <- build_top_annotation(meta_ord, condition_col, region_colors)

              col_fun  <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
              n_row    <- nrow(mat_scaled); n_col <- ncol(mat_scaled)
              res_png  <- res_dpi
              extra_in <- 3.0                              # title + anno + legend
              max_px   <- 30000                            # png/cairo device ceiling
              row_mm   <- 3                                # desired mm per gene row
              body_mm  <- max(n_row * row_mm, 40)
              dev_h_px <- round((body_mm / 25.4 + extra_in) * res_png)
              if (dev_h_px > max_px) {                      # too tall for the device:
                dev_h_px <- max_px                          # cap and shrink rows to fit
                body_mm  <- (max_px / res_png - extra_in) * 25.4
              }
              row_fs   <- if (n_row > 400) 3 else 5         # smaller font when dense

              ht_uniq <- Heatmap(
                mat_scaled, name = "Z-score", col = col_fun,
                top_annotation = col_anno,
                height = unit(body_mm, "mm"),
                row_split = row_split,
                cluster_row_slices = FALSE,
                cluster_columns = FALSE,
                cluster_rows = TRUE,
                show_row_names = (n_row <= 1500),           # labels unreadable beyond this
                show_column_names = TRUE,
                row_names_gp = gpar(fontsize = row_fs),
                column_names_gp = gpar(fontsize = 7),
                row_title_gp = gpar(fontsize = 8),
                column_title = "Unique markers per level (all genes)",
                heatmap_legend_param = list(title = "Z-score")
              )
              png(file.path(ovr_dir, "unique_markers_heatmap.png"),
                  width = round(max(900, n_col * 55) * px_scale),
                  height = dev_h_px, res = res_png)
              draw(ht_uniq, merge_legend = TRUE)
              dev.off()
              message("    Unique markers heatmap: ", n_row, " genes")
            }
          }
        }, error = function(e) {
          message("    Unique markers heatmap failed: ", conditionMessage(e))
        })
      }

    }, error = function(e) {
      message("    Unique marker identification failed: ", conditionMessage(e))
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
