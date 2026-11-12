library(dplyr)
library(ggplot2)

# Use the same drug name and output directory as 04_DESeq2_example.R.
args <- commandArgs(trailingOnly = TRUE)
drug_name <- args[[1]]

result_dir <- file.path(
  "/cluster2/home/xiaorui/project/Click_EpiChem/",
  "RNA_seq_pipeline_from_existing/example_results",
  drug_name
)

drug_file <- file.path(
  result_dir,
  "DESeq2_results_DMSO_vs_drug.csv"
)
biotin_file <- file.path(
  result_dir,
  "DESeq2_results_DMSO_vs_biotin.csv"
)
n3_file <- file.path(
  result_dir,
  "DESeq2_results_DMSO_vs_N3.csv"
)

plot_lfc_correlation <- function(
    probe_file,
    probe_name,
    output_prefix
) {
  probe <- read.csv(probe_file) %>%
    select(X, log2FoldChange) %>%
    rename(probe_lfc = log2FoldChange)

  parent <- read.csv(drug_file) %>%
    select(X, log2FoldChange) %>%
    rename(drug_lfc = log2FoldChange)

  plot_data <- inner_join(probe, parent, by = "X") %>%
    filter(!is.na(probe_lfc), !is.na(drug_lfc))

  pearson_r <- cor(
    plot_data$probe_lfc,
    plot_data$drug_lfc,
    method = "pearson"
  )

  plot_data <- plot_data %>%
    mutate(
      color = case_when(
        probe_lfc < -1 | drug_lfc < -1 ~ "#4DBBD5FF",
        probe_lfc > 1 | drug_lfc > 1 ~ "#E64B35FF",
        TRUE ~ "black"
      )
    )

  correlation_plot <- ggplot(
    plot_data,
    aes(probe_lfc, drug_lfc, fill = color)
  ) +
    geom_rect(
      aes(xmin = -Inf, xmax = -1, ymin = -Inf, ymax = -1),
      fill = "#DDDDEE",
      alpha = 0.1
    ) +
    geom_rect(
      aes(xmin = 1, xmax = Inf, ymin = 1, ymax = Inf),
      fill = "#F6E3D7",
      alpha = 0.1
    ) +
    geom_point(size = 0.8, shape = 21, color = "black", stroke = 0.01) +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    geom_hline(
      yintercept = c(-1, 1),
      color = "#00A087FF",
      linetype = "dashed"
    ) +
    geom_vline(
      xintercept = c(-1, 1),
      color = "#00A087FF",
      linetype = "dashed"
    ) +
    scale_fill_manual(values = c("#4DBBD5FF", "#E64B35FF", "black")) +
    scale_x_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 2)) +
    scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 2)) +
    labs(
      x = paste0("LFC_", probe_name, "_vs_DMSO"),
      y = paste0("LFC_", drug_name, "_vs_DMSO"),
      subtitle = paste0("Pearson Cor = ", round(pearson_r, 3))
    ) +
    guides(fill = "none") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 20, hjust = 1),
      plot.subtitle = element_text(size = 10, vjust = 1)
    )

  write.csv(
    plot_data,
    file.path(result_dir, paste0(output_prefix, "_plot_data.csv")),
    row.names = FALSE
  )

  ggsave(
    file.path(result_dir, paste0(output_prefix, ".pdf")),
    plot = correlation_plot,
    width = 4,
    height = 4
  )
}

plot_lfc_correlation(
  probe_file = biotin_file,
  probe_name = paste0(drug_name, "-biotin"),
  output_prefix = paste0(drug_name, "_biotin_vs_", drug_name)
)

plot_lfc_correlation(
  probe_file = n3_file,
  probe_name = paste0(drug_name, "-N3"),
  output_prefix = paste0(drug_name, "_N3_vs_", drug_name)
)
