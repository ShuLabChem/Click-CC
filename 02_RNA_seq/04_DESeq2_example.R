library(DESeq2)
library(dplyr)

# Example drug and input files
args <- commandArgs(trailingOnly = TRUE)
drug_name <- args[[1]]

count_file <- paste0(
  "/cluster2/home/xiaorui/project/Click_EpiChem/",
  "Bulk_RNA_20260513/1.pipeline/RNA_count_and_length.txt"
)
annotation_file <- paste0(
  "/cluster2/home/xiaorui/reference/hg19/",
  "protein_coding_genes.txt"
)
output_dir <- file.path(
  "/cluster2/home/xiaorui/project/Click_EpiChem/",
  "RNA_seq_pipeline_from_existing/example_results",
  drug_name
)

# Sample columns used for this drug example. Numeric sample names in the input
control_samples <- c("X71", "X72", "X73")
drug_samples <- c("X81", "X82", "X83")
biotin_samples <- c("X11", "X12", "X13")
n3_samples <- c("X21", "X22", "X23")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read the featureCounts matrix and retain protein-coding genes.
counts <- read.table(count_file, sep = "\t", header = TRUE)

protein_coding <- read.table(
  annotation_file,
  sep = "\t",
  header = FALSE
)
colnames(protein_coding) <- c("Geneid", "Gene")

counts <- merge(counts, protein_coding, by = "Geneid") %>%
  filter(!is.na(Gene)) %>%
  select(-Length) %>%
  group_by(Gene) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>%
  as.data.frame()

rownames(counts) <- counts$Gene
counts <- counts %>% select(-Gene)

selected_samples <- c(
  control_samples,
  drug_samples,
  biotin_samples,
  n3_samples
)
missing_samples <- setdiff(selected_samples, colnames(counts))

if (length(missing_samples) > 0) {
  stop(
    "The following samples are absent from the count matrix: ",
    paste(missing_samples, collapse = ", ")
  )
}

run_comparison <- function(
    treatment_samples,
    treatment_name,
    output_name
) {
  comparison_samples <- c(control_samples, treatment_samples)
  comparison_groups <- c(
    rep("DMSO", length(control_samples)),
    rep(treatment_name, length(treatment_samples))
  )

  col_data <- data.frame(
    row.names = comparison_samples,
    condition = factor(
      comparison_groups,
      levels = c("DMSO", treatment_name)
    )
  )

  dds <- DESeqDataSetFromMatrix(
    countData = counts[, comparison_samples, drop = FALSE],
    colData = col_data,
    design = ~condition
  )

  dds <- dds[rowSums(DESeq2::counts(dds)) > 10, ]
  dds <- DESeq(dds)

  result <- results(
    dds,
    contrast = c("condition", treatment_name, "DMSO")
  )
  result_df <- as.data.frame(result) %>%
    arrange(padj, pvalue)

  write.csv(
    result_df,
    file.path(output_dir, output_name),
    quote = FALSE
  )
}

# Drug, biotin-labelled drug, and N3-labelled drug versus DMSO.
run_comparison(
  treatment_samples = drug_samples,
  treatment_name = drug_name,
  output_name = "DESeq2_results_DMSO_vs_drug.csv"
)

run_comparison(
  treatment_samples = biotin_samples,
  treatment_name = paste0(drug_name, "-biotin"),
  output_name = "DESeq2_results_DMSO_vs_biotin.csv"
)

run_comparison(
  treatment_samples = n3_samples,
  treatment_name = paste0(drug_name, "-N3"),
  output_name = "DESeq2_results_DMSO_vs_N3.csv"
)
