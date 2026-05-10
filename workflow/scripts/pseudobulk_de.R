################################################################################
## pseudobulk_de.R – DESeq2 Wald + LRT, EnhancedVolcano, DEGpatterns
##
## Reads pseudobulk count matrices from aggregated/{matrices,metadata}/.
## All result tables ordered by decreasing abs(log2FoldChange).
################################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(EnhancedVolcano)
  library(pheatmap)
  library(RColorBrewer)
  library(tidyverse)
  library(DEGreport)
})

# ── Snakemake params ─────────────────────────────────────────────────────────
agg_dir        <- snakemake@input[["agg_dir"]]
results_dir    <- snakemake@output[["results_dir"]]
min_replicates <- as.integer(snakemake@params[["min_replicates"]])
de_n_genes     <- as.integer(snakemake@params[["de_n_genes"]])
padj_thr       <- as.numeric(snakemake@params[["padj_threshold"]])
lfc_thr        <- as.numeric(snakemake@params[["lfc_threshold"]])

# ── Redirect logs safely ─────────────────────────────────────────────────────
log_out <- snakemake@log[["out"]]
log_err <- snakemake@log[["err"]]

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

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

message("=== Pseudobulk DE (R) ===")
message("  agg_dir:     ", agg_dir)
message("  results_dir: ", results_dir)
message("  min_reps:    ", min_replicates)
message("  de_n_genes:  ", de_n_genes)
message("  padj_thr:    ", padj_thr)
message("  lfc_thr:     ", lfc_thr)


# ── Helper: run DE on one pseudobulk subgroup ────────────────────────────────
run_de_for_prefix <- function(prefix, condition_col, sample_col, out_base) {

  # Files are now in matrices/ and metadata/ subdirectories
  counts_file   <- file.path(agg_dir, "matrices", paste0("counts_", prefix, ".tsv"))
  metadata_file <- file.path(agg_dir, "metadata", paste0("metadata_", prefix, ".tsv"))

  if (!file.exists(counts_file) || !file.exists(metadata_file)) {
    message("  Skipping '", prefix, "': files not found.")
    message("    Looked for: ", counts_file)
    message("    Looked for: ", metadata_file)
    return(NULL)
  }

  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

  # Read data
  counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
  meta   <- read.delim(metadata_file, row.names = 1, check.names = FALSE)

  # Ensure integer counts
  counts <- round(counts)
  counts <- as.matrix(counts)
  storage.mode(counts) <- "integer"

  # Align
  common <- intersect(rownames(counts), rownames(meta))
  if (length(common) == 0) {
    message("  '", prefix, "': no common sample names between counts and metadata. Skipping.")
    return(NULL)
  }
  counts <- counts[common, , drop = FALSE]
  meta   <- meta[common, , drop = FALSE]

  # Ensure condition is a factor
  meta[[condition_col]] <- factor(meta[[condition_col]])
  conditions <- levels(meta[[condition_col]])

  # Check replicates
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

  # Remove zero-count genes
  gene_sums <- colSums(counts)
  counts <- counts[, gene_sums > 0, drop = FALSE]

  message("  '", prefix, "': ", nrow(counts), " samples, ", ncol(counts), " genes, ",
          length(conditions), " conditions: ", paste(conditions, collapse = ", "))

  # Transpose: DESeq2 expects genes as rows
  counts_t <- t(counts)

  # Verify alignment
  if (!identical(colnames(counts_t), rownames(meta))) {
    message("  WARNING: sample name mismatch after transpose. Attempting reorder.")
    meta <- meta[colnames(counts_t), , drop = FALSE]
  }

  # ── DESeq2 dataset ──────────────────────────────────────────────────────
  design_formula <- as.formula(paste("~", condition_col))
  dds <- DESeqDataSetFromMatrix(
    countData = counts_t,
    colData   = meta,
    design    = design_formula
  )

  # ── 1. Holistic LRT ────────────────────────────────────────────────────
  message("    LRT …")
  dds_lrt <- NULL
  tryCatch({
    dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)
    res_lrt <- results(dds_lrt)
    res_lrt_df <- as.data.frame(res_lrt) %>%
      rownames_to_column("gene") %>%
      arrange(desc(abs(log2FoldChange)))
    write.table(res_lrt_df, file.path(out_base, "LRT_results.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    n_sig <- sum(res_lrt_df$padj < padj_thr, na.rm = TRUE)
    message("    LRT: ", n_sig, " significant genes (padj < ", padj_thr, ")")
  }, error = function(e) {
    message("    LRT failed: ", conditionMessage(e))
  })

  # ── 2. Pairwise Wald tests ─────────────────────────────────────────────
  message("    Wald pairwise tests …")
  dds_wald <- DESeq(dds, test = "Wald")

  # Normalised counts for heatmaps
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
        arrange(desc(abs(log2FoldChange)))

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

      # ── EnhancedVolcano ──────────────────────────────────────────────
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

      # ── Top-N heatmap ────────────────────────────────────────────────
      tryCatch({
        sig_genes <- res_df %>% filter(padj < padj_thr)
        top_up   <- sig_genes %>% filter(log2FoldChange > 0) %>% head(de_n_genes) %>% pull(gene)
        top_down <- sig_genes %>% filter(log2FoldChange < 0) %>% head(de_n_genes) %>% pull(gene)
        top_genes <- c(top_up, top_down)
        top_genes <- top_genes[top_genes %in% rownames(assay(vsd))]

        if (length(top_genes) >= 2) {
          mat <- assay(vsd)[top_genes, , drop = FALSE]
          mat_scaled <- t(scale(t(mat)))

          annotation_col <- data.frame(
            condition = meta[[condition_col]],
            row.names = colnames(mat)
          )

          png(file.path(pairwise_dir, paste0("heatmap_top", de_n_genes, "_", contrast_name, ".png")),
              width = max(800, ncol(mat) * 60), height = max(400, length(top_genes) * 25), res = 150)
          pheatmap(
            mat_scaled,
            annotation_col = annotation_col,
            cluster_rows = TRUE, cluster_cols = TRUE,
            show_rownames = TRUE, show_colnames = TRUE,
            main = paste0("Top ", de_n_genes, " DE genes — ", contrast_name),
            fontsize_row = 8, fontsize_col = 8
          )
          dev.off()
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

  # Save pairwise summary
  if (length(summary_rows) > 0) {
    bind_rows(summary_rows) %>%
      write.table(file.path(out_base, "pairwise_summary.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # ── 3. DEGpatterns on LRT significant genes ─────────────────────────────
  if (!is.null(dds_lrt)) {
    tryCatch({
      res_lrt_all <- results(dds_lrt)
      sig_genes <- rownames(res_lrt_all)[which(res_lrt_all$padj < padj_thr)]

      if (length(sig_genes) >= 10) {
        message("    DEGpatterns on ", length(sig_genes), " LRT-significant genes …")
        rlog_mat <- assay(vsd)
        sig_mat  <- rlog_mat[sig_genes, , drop = FALSE]

        deg_dir <- file.path(out_base, "DEGpatterns")
        dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)

        clusters <- degPatterns(
          sig_mat,
          metadata = meta,
          time = condition_col,
          plot = FALSE
        )

        if (!is.null(clusters$df)) {
          cluster_df <- clusters$df %>% arrange(cluster, genes)
          write.table(cluster_df, file.path(deg_dir, "gene_clusters.tsv"),
                      sep = "\t", quote = FALSE, row.names = FALSE)

          p <- degPlotCluster(
            clusters$normalized,
            time = condition_col,
            color = condition_col,
            points = TRUE
          )
          ggsave(file.path(deg_dir, "DEGpatterns_clusters.png"),
                 plot = p, width = 12, height = 8, dpi = 500)
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


# ── Main: read manifest and process each subgroup ────────────────────────────
manifest_file <- file.path(agg_dir, "manifest.tsv")

if (!file.exists(manifest_file)) {
  message("No manifest.tsv found. Nothing to process.")
  writeLines("No pseudobulk matrices found.", file.path(results_dir, "SKIPPED_no_manifest.txt"))
  # Close sinks before exit
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

  out_base <- if (prefix == "pooled") {
    results_dir
  } else {
    file.path(results_dir, prefix)
  }

  run_de_for_prefix(prefix, condition_col, sample_col, out_base)
}

message("\n=== Pseudobulk DE complete ===")

# Close sinks properly (reverse order)
if (exists("log_err_con")) { sink(type = "message"); close(log_err_con) }
if (exists("log_out_con")) { sink(); close(log_out_con) }
