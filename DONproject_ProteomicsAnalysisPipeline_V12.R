# =====================================================================
# DON project proteomics -- V11 FINAL pipeline
# ---------------------------------------------------------------------
# Primary analysis:
#   * limma differential expression, Don vs DCTRL
#   * ALL 24 runs used; technical replicates handled by
#     duplicateCorrelation (repeated measures within bio-sample)
#   * Batch (Mouse Harvest Order, Day1 vs Day2, from the DefinedBatch
#     workbook) enters the design as a covariate:  ~ batch + group
#   * Proteins kept if quantified in >=50% of bio-samples within at
#     least one group (>=2 of 4)
#   * FDR tier      BH adj.p < 0.05
#
# Visualization ONLY (never used for inference):
#   * ComBat batch-adjusted matrix -> PCA before/after, at the 24-run
#     level and the 8-bio-sample level. ComBat output is NEVER passed
#     to lmFit; the DE test always runs on the raw imputed intensities
#     with batch in the design.
# =====================================================================

# ---- User-editable paths -------------------------------------------
INPUT_TSV   <- "/mnt/user-uploads/matrix.Don.tsv"   # DIA-NN pg_matrix export
BATCH_XLSX  <- "/mnt/user-uploads/DonProject_Proteomics_DefinedBatch.xlsx"
OUT_PDF     <- "/workspace/don/DonProject_QC_Report_V11_Final.pdf"
OUT_CSV     <- "/mnt/results/DonProject_Significant_Proteins_V11_Final.csv"

library(limma)
library(sva)        # ComBat for the visualization matrix ONLY
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggrepel)
library(readxl)
library(grid)
library(gridExtra)

set.seed(123)

LOGFC_CUTOFF <- 0.58

add_text_page <- function(title, lines) {
  grid::grid.newpage()
  grid::grid.text(
    label = paste(c(title, "", lines), collapse = "\n"),
    x = 0.05, y = 0.95, just = c("left", "top"),
    gp = grid::gpar(fontsize = 12, fontfamily = "mono")
  )
}

# ===============================
# 1. Load DIA-NN protein table
# ===============================
expr_raw <- read.delim(INPUT_TSV, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

expr <- expr_raw %>%
  select(
    Protein.Group, Protein.Names, Genes,
    N.Sequences, N.Proteotypic.Sequences,
    starts_with("L:\\clinprot\\Astral_RAW_DATA\\Fatemeh\\Project_1\\")
  ) %>%
  mutate(
    Genes = ifelse(is.na(Genes) | Genes == "", Protein.Group, Genes)
  ) %>%
  mutate(across(starts_with("L:\\"), ~ ifelse(. == 0, NA, .)))

expr_clean <- expr %>%
  group_by(Protein.Group) %>%
  summarise(across(where(is.numeric), median, na.rm = TRUE), .groups = "drop") %>%
  distinct(Protein.Group, .keep_all = TRUE)

expr_mat_raw <- expr_clean %>%
  column_to_rownames("Protein.Group") %>%
  select(starts_with("L:\\clinprot\\")) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

cat("Protein groups loaded:", nrow(expr_mat_raw), "\n")

# ===============================
# 2. Metadata + batch
# ===============================
meta_run <- data.frame(sample_name = colnames(expr_mat_raw)) %>%
  mutate(
    raw_file = sub(".*[\\\\/]", "", sample_name),   # portable basename
    group = case_when(
      grepl("DCTRL", raw_file, ignore.case = TRUE) ~ "DCTRL",
      grepl("Don",  raw_file, ignore.case = TRUE)  ~ "Don",
      TRUE ~ "Unknown"),
    bio_id = sub("_\\d\\.raw$", "", raw_file) %>%
      gsub("20260130_Fatemeh_Kidney_", "", .) %>%
      make.names()
  )

batch_raw <- read_excel(BATCH_XLSX, sheet = 1, col_names = FALSE)
batch_tbl <- batch_raw %>%
  select(sample = `...1`, harvest_order = `...3`) %>%
  mutate(sample = str_trim(as.character(sample)),
         harvest_order = suppressWarnings(as.integer(harvest_order))) %>%
  filter(!is.na(harvest_order)) %>%
  distinct(sample, harvest_order)

meta_run$batch <- batch_tbl$harvest_order[match(meta_run$bio_id, batch_tbl$sample)]
meta_run$batch <- factor(meta_run$batch, levels = c(1, 2), labels = c("Day1", "Day2"))
stopifnot(!any(is.na(meta_run$batch)))

cat("\n=== Batch (Harvest Order) x group, run level ===\n")
print(table(meta_run$batch, meta_run$group))

# ===============================
# 3. log2 + Shapiro-gated median-centering
# ===============================
expr_log <- log2(expr_mat_raw + 1)

shapiro_p <- sapply(colnames(expr_log), function(s) {
  x <- expr_log[, s]; x <- x[!is.na(x) & is.finite(x)]
  if (length(x) > 5000) x <- sample(x, 5000)
  shapiro.test(x)$p.value
})
needs_normalization <- !all(shapiro_p > 0.05)
cat("\nShapiro-gated normalization needed:", needs_normalization, "\n")
if (needs_normalization) {
  medians <- apply(expr_log, 2, median, na.rm = TRUE)
  expr_norm <- sweep(expr_log, 2, medians, "-")
} else {
  expr_norm <- expr_log
}

# ===============================
# 4. Validity filter (bio-sample level)
# ===============================
expr_long <- expr_norm %>%
  as.data.frame() %>%
  rownames_to_column("protein") %>%
  pivot_longer(-protein, names_to = "sample_name", values_to = "intensity") %>%
  left_join(meta_run %>% select(sample_name, bio_id, group), by = "sample_name")

expr_collapsed <- expr_long %>%
  group_by(protein, bio_id, group) %>%
  summarise(intensity = median(intensity, na.rm = TRUE), .groups = "drop")

mat8 <- expr_collapsed %>%
  select(protein, bio_id, intensity) %>%
  pivot_wider(names_from = bio_id, values_from = intensity) %>%
  column_to_rownames("protein") %>%
  as.matrix()

meta8 <- expr_collapsed %>% distinct(bio_id, group) %>%
  arrange(match(bio_id, colnames(mat8))) %>% as.data.frame()
mat8 <- mat8[, meta8$bio_id]
meta8$batch <- batch_tbl$harvest_order[match(meta8$bio_id, batch_tbl$sample)]
meta8$batch <- factor(meta8$batch, levels = c(1, 2), labels = c("Day1", "Day2"))

valid_mat <- !is.na(mat8)
grp_of_col <- setNames(meta8$group, meta8$bio_id)
n_valid_D <- rowSums(valid_mat[, grp_of_col[colnames(mat8)] == "DCTRL"])
n_valid_T <- rowSums(valid_mat[, grp_of_col[colnames(mat8)] == "Don"])

keep <- (n_valid_D >= 2) | (n_valid_T >= 2)   # >=50% in >=1 group
cat("Proteins passing validity filter (>=2 of 4 bio-samples in >=1 group):",
    sum(keep), "of", length(keep), "\n")


# ===============================
# 5. Differential expression
#    (24 runs, duplicateCorrelation, ~ batch + group)
# ===============================
design <- model.matrix(~ batch + group, data = meta_run)
stopifnot(qr(design)$rank == ncol(design))

corfit <- duplicateCorrelation(imp24, design, block = meta_run$bio_id)
cat("\nConsensus within-bio-sample correlation (tech reps):",
    round(corfit$consensus.correlation, 3), "\n")

fit  <- lmFit(imp24, design, block = meta_run$bio_id,
              correlation = corfit$consensus.correlation)
fit2 <- eBayes(fit, trend = FALSE, robust = TRUE)
tt   <- topTable(fit2, coef = "groupDon", number = Inf) %>%
  rownames_to_column("protein_group")

results_final <- tt %>%
  left_join(
    expr_raw %>% janitor::clean_names() %>%
      distinct(protein_group, .keep_all = TRUE) %>%
      select(protein_group, genes, first_protein_description),
    by = "protein_group"
  ) %>%
  mutate(
    gene_symbol = coalesce(sub(";.*", "", genes), protein_group),
    significant_FDR      = adj.P.Val < 0.05 & abs(logFC) > LOGFC_CUTOFF,
    significant_nominal  = P.Value   < 0.05 & abs(logFC) > LOGFC_CUTOFF
  )

n_tested <- nrow(results_final)
n_fdr    <- sum(results_final$significant_FDR)
n_nom    <- sum(results_final$significant_nominal)
n_up     <- sum(results_final$significant_FDR & results_final$logFC > 0)
n_down   <- sum(results_final$significant_FDR & results_final$logFC < 0)

cat("\nTotal proteins tested:", n_tested, "\n")
cat("FDR-tier DEPs (adj.p < 0.05 & |log2FC| >", LOGFC_CUTOFF, "):", n_fdr,
    " (up:", n_up, " down:", n_down, ")\n")
cat("Nominal-tier DEPs (p < 0.05 & |log2FC| >", LOGFC_CUTOFF, "):", n_nom,
    " -- uncorrected, interpretation only\n")
cat("Model: limma ~ batch + group on 24 runs (duplicateCorrelation);\n",
    "ComBat used for the PCA visualization ONLY. 4 bio-replicates per\n",
    "group, each a pool of 2 animals -- state this wherever n is quoted.\n")

write.csv(results_final, OUT_CSV, row.names = FALSE)

# ===============================
# 6. Visualization: ComBat PCA (BEFORE/AFTER) -- visualization ONLY
# ===============================
# ComBat removes the batch (day) component while preserving group.
# The result feeds PCA plots ONLY.
combat24 <- ComBat(dat = imp24, batch = meta_run$batch,
                   mod = model.matrix(~ group, data = meta_run))
combat8  <- ComBat(dat = imp8,  batch = meta8$batch,
                   mod = model.matrix(~ group, data = meta8))

pca_panel <- function(mat, meta, title) {
  pca <- prcomp(t(mat), scale. = FALSE)
  ve  <- pca$sdev^2 / sum(pca$sdev^2)
  df  <- data.frame(pca$x[, 1:2], group = meta$group, batch = meta$batch,
                    sample = sub(".*Kidney_", "", colnames(mat)))
  ggplot(df, aes(PC1, PC2, color = group, shape = batch)) +
    geom_point(size = 4, alpha = 0.75) +
    ggrepel::geom_text_repel(aes(label = sample), size = 3,
                              box.padding = 0.5, point.padding = 0.4,
                              max.overlaps = Inf, segment.color = "grey40",
                              min.segment.length = 0) +
    xlab(paste0("PC1 (", round(ve[1]*100, 1), "%)")) +
    ylab(paste0("PC2 (", round(ve[2]*100, 1), "%)")) +
    ggtitle(title) +
    theme_minimal()
}

pdf(OUT_PDF, width = 10, height = 7)

print(pca_panel(imp24,    meta_run, "PCA run level - BEFORE ComBat (raw imputed)"))
print(pca_panel(combat24, meta_run, "PCA run level - AFTER ComBat (visualization only)"))
print(pca_panel(imp8,     meta8,    "PCA bio-sample level - BEFORE ComBat (raw imputed)"))
print(pca_panel(combat8,  meta8,    "PCA bio-sample level - AFTER ComBat (visualization only)"))

print(
  ggplot(data.frame(p = tt$P.Value), aes(x = p)) +
    geom_histogram(bins = 40, fill = "grey30", color = "white") +
    ggtitle("p-value histogram - limma ~ batch + group (24 runs)") +
    xlab("raw p-value") + theme_minimal()
)

completeness <- tibble(protein_group = rownames(mat8),
                       n_valid_DCTRL = n_valid_D[keep],
                       n_valid_Don   = n_valid_T[keep]) %>%
  count(n_valid_DCTRL, n_valid_Don)
print(
  ggplot(completeness, aes(n_valid_DCTRL, n_valid_Don, fill = n)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c() +
    xlab("# valid DCTRL bio-samples (of 4)") +
    ylab("# valid Don bio-samples (of 4)") +
    ggtitle("Protein quantification completeness") +
    theme_minimal()
)

add_text_page(
  "V11 final pipeline -- summary",
  c(
    "Input: DIA-NN pg_matrix (24 runs = 8 bio-samples x 3 tech reps)",
    "Batch: Mouse Harvest Order (Day1/Day2), covariate in the design",
    "Model: limma ~ batch + group, 24 runs, duplicateCorrelation",
    "ComBat: visualization (PCA) ONLY -- never used for inference",
    "",
    paste("Proteins tested:", n_tested),
    paste("FDR-tier DEPs (adj.p<0.05 & |log2FC|>0.58):", n_fdr,
          paste0("(", n_up, " up / ", n_down, " down)")),
    paste("Nominal-tier DEPs (p<0.05, uncorrected):", n_nom),
    "",
    paste("Consensus tech-rep correlation:",
          round(corfit$consensus.correlation, 3))
  )
)

dev.off()
cat("\nQC PDF saved to:", OUT_PDF, "\n")
cat("DE CSV saved to:", OUT_CSV, "\n")
cat("Done.\n")
