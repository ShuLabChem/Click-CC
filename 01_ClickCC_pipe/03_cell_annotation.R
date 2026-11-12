#!/usr/bin/env Rscript
# Existing paired ATAC/drug analysis with sample names supplied as metadata.
# Usage: Rscript 03_cell_annotation.R INPUT_DIR CELLS.tsv ANNOTATION.rds OUTPUT_DIR [DATASET_ID]
args <- commandArgs(trailingOnly = TRUE)

library(Seurat)
library(Signac)
library(rtracklayer)
library(ggplot2)
library(tidyverse)
library(scDblFinder)
library(patchwork)
library(qs)

options(future.globals.maxSize = 10 * 1024^3)
set.seed(1234)
data_dir <- normalizePath(args[[1]], mustWork = TRUE)
cell_metadata <- read.delim(args[[2]], stringsAsFactors = FALSE, check.names = FALSE)
annotation_rds <- normalizePath(args[[3]], mustWork = TRUE)
output_dir <- args[[4]]
dataset_id <- if (length(args) >= 5) args[[5]] else "dataset"
timepoint_levels <- c("Day1", "Day12")
required <- c("barcode", "sample_id", "timepoint")
if (!all(required %in% names(cell_metadata))) stop("Cell metadata requires: barcode, sample_id, timepoint")
if (anyDuplicated(cell_metadata$barcode)) stop("Cell barcodes must be unique")
if (anyNA(cell_metadata[, required]) || any(!nzchar(as.matrix(cell_metadata[, required])))) stop("Cell metadata contains missing values")
if (any(!cell_metadata$timepoint %in% timepoint_levels)) stop("Timepoints must match the retained Day1/Day12 QC settings")
rownames(cell_metadata) <- cell_metadata$barcode
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)
dir.create("Rds", showWarnings = FALSE)

# ATAC: peak matrices, QC and doublet removal
group_mapping <- setNames(list(data_dir), dataset_id)

peaks_list <- list()

for (sample_id in names(group_mapping)) {
  path_prefix <- group_mapping[[sample_id]]
  search_pattern <- file.path(path_prefix, "bulk", "04_macs2", "*.ATAC_peaks.narrowPeak")
  peak_path <- Sys.glob(search_pattern)

  if(length(peak_path) > 0){
    message(paste("Reading peaks for:", sample_id))
    peaks_list[[sample_id]] <- import(peak_path)
  } else {
    warning(paste("Peak file not found for:", sample_id, "-", peak_path))
  }
}

message("Reducing peaks...")
peak_union <- Signac::reduce(do.call(c, unlist(peaks_list, use.names = FALSE)))
peakwidths <- width(peak_union)
chromosomes <- seqnames(peak_union)
standard_chroms <- c(paste0("chr", 1:22), "chrX", "chrY")
peak_union <- peak_union[peakwidths < 10000 & peakwidths > 20 & chromosomes %in% standard_chroms]
message(paste("Final union peak set size:", length(peak_union)))
annotations <- readRDS(annotation_rds)

seu_list <- list()
for (sample_id in names(group_mapping)) {
  path_prefix <- group_mapping[[sample_id]]

  message(paste("Processing sample:", sample_id))
  frag_path <- file.path(path_prefix, "sinto", "06_fragTab", "ATAC_merge.sort.tsv.gz")

  if (!file.exists(frag_path)) {
    stop(paste("Fragment file not found:", frag_path))
  }
  frags <- CreateFragmentObject(path = frag_path)
  counts <- FeatureMatrix(fragments = frags, features = peak_union)
  assay_obj <- CreateChromatinAssay(counts = counts, fragments = frags)
  seu_obj <- CreateSeuratObject(assay_obj, assay = "ATAC")
  seu_obj$orig.ident <- sample_id
  seu_obj$dataset <- sample_id
  Annotation(seu_obj[["ATAC"]]) <- annotations
  seu_obj <- NucleosomeSignal(seu_obj)
  seu_obj <- TSSEnrichment(seu_obj, fast = FALSE)
  seu_list[[sample_id]] <- seu_obj
}

seu_obj <- seu_list[[dataset_id]]
cells_use <- intersect(Cells(seu_obj), cell_metadata$barcode)
if (!length(cells_use)) stop("No input barcodes match the cell metadata")
seu_obj <- subset(seu_obj, cells = cells_use)
seu_obj$timepoint <- factor(cell_metadata[Cells(seu_obj), "timepoint"], levels = timepoint_levels)
seu_obj$sample_id <- cell_metadata[Cells(seu_obj), "sample_id"]
table(seu_obj$timepoint)
seu_list[[dataset_id]] <- seu_obj
qsave(seu_list,"Rds/seu_list.qs")

dir.create("Rds")
dir.create("QC_plots/ATAC",recursive = T)

for (sample_id in names(seu_list)) {
  seu_obj <- seu_list[[sample_id]]
  p <- VlnPlot(object = seu_obj, features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"), pt.size = 0, ncol = 4)
  ggsave(filename = paste0("QC_plots/ATAC/", sample_id, "_QC_vlnplot.pdf"),plot = p,width = 8,height = 5)
}

seu_obj <- subset(seu_obj,nCount_ATAC < 5000 & nFeature_ATAC < 10000)

p <- VlnPlot(
  object = seu_obj,
  features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"),
  group.by = "timepoint",
  pt.size = 0,
  ncol = 4
)
ggsave("QC_plots/ATAC/QC_vlnplot_by_TreatTime.pdf", plot = p, width = 12, height = 5)

qc_cutoffs_tp <- list(
  Day1  = list(min_count = 100, max_count = 2000, min_feature = 50, max_feature = 2000, min_tss = 2, max_tss = 15, max_nucleosome = 1),
  Day12 = list(min_count = 100, max_count = 2000, min_feature = 50, max_feature = 2000, min_tss = 2, max_tss = 10, max_nucleosome = 1.5)
)

md <- seu_obj@meta.data

cells_keep <- c()

for (tp in names(qc_cutoffs_tp)) {
  th <- qc_cutoffs_tp[[tp]]

  keep_tp <- rownames(md)[
    md$timepoint == tp &
    md$nCount_ATAC > th$min_count &
    md$nCount_ATAC < th$max_count &
    md$nFeature_ATAC > th$min_feature &
    md$nFeature_ATAC < th$max_feature &
    md$TSS.enrichment > th$min_tss &
    md$TSS.enrichment < th$max_tss &
    md$nucleosome_signal < th$max_nucleosome
  ]

  message(tp, ": ", sum(md$timepoint == tp), " -> ", length(keep_tp), " cells after QC")
  cells_keep <- c(cells_keep, keep_tp)
}

seu_filt <- subset(seu_obj, cells = unique(cells_keep))
filtered_list <- setNames(list(seu_filt), dataset_id)

p <- VlnPlot(
  object = seu_filt,
  features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"),
  group.by = "timepoint",
  pt.size = 0,
  ncol = 4
)
ggsave("QC_plots/ATAC/QC_vlnplot_after_filter_by_TreatTime.pdf", plot = p, width = 12, height = 5)

qsave(filtered_list, "Rds/ATAC_filtered_list.qs")

singlet_list <- list()

for (sample_id in names(filtered_list)) {
  seu_obj <- filtered_list[[sample_id]]
  DefaultAssay(seu_obj) <- "ATAC"

  sce <- as.SingleCellExperiment(seu_obj)

  sce <- scDblFinder(sce)

  seu_obj$scDblFinder.score <- colData(sce)$scDblFinder.score
  seu_obj$scDblFinder.class <- colData(sce)$scDblFinder.class

  seu_singlet <- subset(seu_obj, subset = scDblFinder.class == "singlet")
  singlet_list[[sample_id]] <- seu_singlet

  message(sample_id, ": ",
          ncol(seu_obj), " -> ", ncol(seu_singlet), " cells after doublet removal")
}

if (length(singlet_list) > 1) {
  ATAC_final <- merge(
    x = singlet_list[[1]],
    y = singlet_list[-1],
    add.cell.ids = names(singlet_list)
  )
} else {
  ATAC_final <- singlet_list[[1]]
}

print(ATAC_final)
dir.create("Rds")
qsave(ATAC_final,"Rds/ATAC_obj0.qs")
ATAC <- ATAC_final

# drug: peak matrices, QC and doublet removal
group_mapping <- setNames(list(data_dir), dataset_id)

peaks_list <- list()

for (sample_id in names(group_mapping)) {
  path_prefix <- group_mapping[[sample_id]]
  search_pattern <- file.path(path_prefix, "bulk", "04_macs2", "*.drug_peaks.narrowPeak")
  peak_path <- Sys.glob(search_pattern)

  if(length(peak_path) > 0){
    message(paste("Reading peaks for:", sample_id))
    peaks_list[[sample_id]] <- import(peak_path)
  } else {
    warning(paste("Peak file not found for:", sample_id, "-", peak_path))
  }
}

message("Reducing peaks...")
peak_union <- Signac::reduce(do.call(c, unlist(peaks_list, use.names = FALSE)))
peakwidths <- width(peak_union)
chromosomes <- seqnames(peak_union)
standard_chroms <- c(paste0("chr", 1:22), "chrX", "chrY")
peak_union <- peak_union[peakwidths < 10000 & peakwidths > 20 & chromosomes %in% standard_chroms]
message(paste("Final union peak set size:", length(peak_union)))

seu_list <- list()
for (sample_id in names(group_mapping)) {
  path_prefix <- group_mapping[[sample_id]]

  message(paste("Processing sample:", sample_id))
  frag_path <- file.path(path_prefix, "sinto", "06_fragTab", "drug_merge.sort.tsv.gz")
  if (!file.exists(frag_path)) {
    stop(paste("Fragment file not found:", frag_path))
  }
  frags <- CreateFragmentObject(path = frag_path)
  counts <- FeatureMatrix(fragments = frags, features = peak_union)
  assay_obj <- CreateChromatinAssay(counts = counts, fragments = frags)
  seu_obj <- CreateSeuratObject(assay_obj, assay = "drug")
  Annotation(seu_obj[["drug"]]) <- annotations
  seu_obj <- NucleosomeSignal(seu_obj)
  seu_obj <- TSSEnrichment(seu_obj, fast = FALSE)
  seu_obj$orig.ident <- sample_id
  seu_obj$dataset <- sample_id
  seu_list[[sample_id]] <- seu_obj
}

dir.create("QC_plots/drug",recursive = T)
for (sample_id in names(seu_list)) {
  seu_obj <- seu_list[[sample_id]]
  p <- VlnPlot(object = seu_obj, features = c("nCount_drug", "nFeature_drug", "TSS.enrichment", "nucleosome_signal"), pt.size = 0, ncol = 4)
  ggsave(filename = paste0("QC_plots/drug/", sample_id, "_QC_vlnplot.pdf"),plot = p,width = 8,height = 5)
}

seu_obj <- seu_list[[dataset_id]]
cells_use <- intersect(Cells(seu_obj), cell_metadata$barcode)
if (!length(cells_use)) stop("No input barcodes match the cell metadata")
seu_obj <- subset(seu_obj, cells = cells_use)
seu_obj$timepoint <- factor(cell_metadata[Cells(seu_obj), "timepoint"], levels = timepoint_levels)
seu_obj$sample_id <- cell_metadata[Cells(seu_obj), "sample_id"]
table(seu_obj$timepoint)
seu_list[[dataset_id]] <- seu_obj
qsave(seu_list,"Rds/seu_list_drug.qs")

dir.create("QC_plots/drug",recursive = T)

for (sample_id in names(seu_list)) {
  seu_obj <- seu_list[[sample_id]]
  p <- VlnPlot(object = seu_obj, features = c("nCount_drug", "nFeature_drug", "TSS.enrichment", "nucleosome_signal"), pt.size = 0, ncol = 4)
  ggsave(filename = paste0("QC_plots/drug/", sample_id, "_QC_vlnplot.pdf"),plot = p,width = 8,height = 5)
}

seu_obj <- subset(seu_obj,nCount_drug < 1000 & nFeature_drug < 1000)

p <- VlnPlot(
  object = seu_obj,
  features = c("nCount_drug", "nFeature_drug", "TSS.enrichment", "nucleosome_signal"),
  group.by = "timepoint",
  pt.size = 0,
  ncol = 4
)
ggsave("QC_plots/drug/QC_vlnplot_by_TreatTime.pdf", plot = p, width = 12, height = 5)

qc_cutoffs_tp <- list(
  Day1  = list(min_count = 10, max_count = 250, min_feature = 10, max_feature = 250, min_tss = 1, max_tss = 10, max_nucleosome = 1),
  Day12 = list(min_count = 10, max_count = 250, min_feature = 10, max_feature = 250, min_tss = 1, max_tss = 10, max_nucleosome = 1.5)
)

md <- seu_obj@meta.data

cells_keep <- c()

for (tp in names(qc_cutoffs_tp)) {
  th <- qc_cutoffs_tp[[tp]]

  keep_tp <- rownames(md)[
    md$timepoint == tp &
    md$nCount_drug > th$min_count &
    md$nCount_drug < th$max_count &
    md$nFeature_drug > th$min_feature &
    md$nFeature_drug < th$max_feature &
    md$TSS.enrichment > th$min_tss &
    md$TSS.enrichment < th$max_tss &
    md$nucleosome_signal < th$max_nucleosome
  ]

  message(tp, ": ", sum(md$timepoint == tp), " -> ", length(keep_tp), " cells after QC")
  cells_keep <- c(cells_keep, keep_tp)
}

seu_filt <- subset(seu_obj, cells = unique(cells_keep))
filtered_list <- setNames(list(seu_filt), dataset_id)

p <- VlnPlot(
  object = seu_filt,
  features = c("nCount_drug", "nFeature_drug", "TSS.enrichment", "nucleosome_signal"),
  group.by = "timepoint",
  pt.size = 0,
  ncol = 4
)
ggsave("QC_plots/drug/QC_vlnplot_after_filter_by_TreatTime.pdf", plot = p, width = 12, height = 5)

qsave(filtered_list, "Rds/drug_filtered_list.qs")

singlet_list <- list()

for (sample_id in names(filtered_list)) {
  seu_obj <- filtered_list[[sample_id]]
  DefaultAssay(seu_obj) <- "drug"

  sce <- as.SingleCellExperiment(seu_obj)

  sce <- scDblFinder(sce)

  seu_obj$scDblFinder.score <- colData(sce)$scDblFinder.score
  seu_obj$scDblFinder.class <- colData(sce)$scDblFinder.class

  seu_singlet <- subset(seu_obj, subset = scDblFinder.class == "singlet")
  singlet_list[[sample_id]] <- seu_singlet

  message(sample_id, ": ",
          ncol(seu_obj), " -> ", ncol(seu_singlet), " cells after doublet removal")
}

if (length(singlet_list) > 1) {
  drug_final <- merge(
    x = singlet_list[[1]],
    y = singlet_list[-1],
    add.cell.ids = names(singlet_list)
  )
} else {
  drug_final <- singlet_list[[1]]
}

print(drug_final)
qsave(drug_final,"Rds/drug_obj0.qs")
drug <- drug_final
# Merge paired assays and add gene activity
omics2_cells <- intersect(Cells(drug),Cells(ATAC))
length(omics2_cells)
library(data.table)
library(Seurat)

drug <- subset(drug,cells = omics2_cells)
ATAC <- subset(ATAC,cells = omics2_cells)
obj <- ATAC
obj[["drug"]] <-  CreateChromatinAssay(counts = GetAssayData(drug, slot = "counts"),fragments = Fragments(drug))
obj$timepoint <- factor(obj$timepoint, levels = timepoint_levels)

annotations <- readRDS(annotation_rds)
Annotation(obj[["ATAC"]]) <- annotations
Annotation(obj[["drug"]]) <- annotations
qsave(obj,"Rds/obj_merge_drug_ATAC.qs")

DefaultAssay(obj) <- "ATAC"
ga <- GeneActivity(obj)
obj[["ATAC_genes"]] <- CreateAssayObject(counts = ga)
obj <- NormalizeData(obj,assay = "ATAC_genes",normalization.method = "LogNormalize",scale.factor = median(obj@meta.data[["nCount_ATAC_genes"]]))

DefaultAssay(obj) <- "drug"
ga <- GeneActivity(obj)
obj[["drug_genes"]] <- CreateAssayObject(counts = ga)
obj <- NormalizeData(obj,assay = "drug_genes",normalization.method = "LogNormalize",scale.factor = median(obj@meta.data[["nCount_drug_genes"]]))

qsave(obj,"Rds/obj1.qs")

# ATAC dimensional reduction and integration

DefaultAssay(obj) <- "ATAC"
obj <- FindTopFeatures(obj, min.cutoff = 5)
obj <- RunTFIDF(obj)
obj <- RunSVD(obj,reduction.name = "lsi.atac", reduction.key = "LSIATAC_")
ggsave("ATAC_depthCor.png",DepthCor(obj,reduction = "lsi.atac"))

obj <- RunUMAP(object = obj, dims = 2:30, reduction = "lsi.atac", reduction.name = "umap.atac")
obj <- FindNeighbors(object = obj, reduction = 'lsi.atac', dims = 2:30)
obj <- FindClusters(object = obj, verbose = FALSE, algorithm = 3)
p <- DimPlot(obj,group.by = "timepoint") + DimPlot(obj)
ggsave("ATAC_umap.png",p,height = 7,width = 14)

obj.list <- SplitObject(obj, split.by = "sample_id")

for (i in 1:length(obj.list)) {
  DefaultAssay(obj.list[[i]]) <- "ATAC"
  obj.list[[i]] <- FindTopFeatures(obj.list[[i]], min.cutoff = 5)
  obj.list[[i]] <- RunTFIDF(obj.list[[i]])
  obj.list[[i]] <- RunSVD(obj.list[[i]], reduction.name = "lsi", reduction.key = "LSI_")
}

integration.anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  anchor.features = rownames(obj),
  reduction = "rlsi",
  assay = rep("ATAC",length(obj.list)),
  dims = 2:30
)

integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = obj[["lsi.atac"]],
  new.reduction.name = "integrated_lsi.atac",
  dims.to.integrate = 1:30
)

integrated <- RunUMAP(integrated, reduction = "integrated_lsi.atac", dims = 2:30)
qsave(integrated,"Rds/integrated.ATAC.qs")
p2 <- DimPlot(integrated, group.by = "timepoint") + DimPlot(integrated)
ggsave("ATAC_inte_umap.png",p2,width = 14,height =5)
obj[["integrated_lsi.atac"]] <- integrated[["integrated_lsi.atac"]]
qsave(obj,"Rds/obj.ATAC.qs")

# drug dimensional reduction and integration
DefaultAssay(obj) <- "drug"
obj <- FindTopFeatures(obj, min.cutoff = 5)
obj <- RunTFIDF(obj)
obj <- RunSVD(obj,reduction.name = "lsi.drug", reduction.key = "LSIDRUG_")
ggsave("drug_depthCor.png",DepthCor(obj,reduction = "lsi.drug"))

obj <- RunUMAP(object = obj, dims = 2:30, reduction = "lsi.drug", reduction.name = "umap.drug")
obj <- FindNeighbors(object = obj, reduction = 'lsi.drug', dims = 2:30)
obj <- FindClusters(object = obj, verbose = FALSE, algorithm = 3)
p <- DimPlot(obj,group.by = "timepoint") + DimPlot(obj)
ggsave("drug_umap.png",p,height = 7,width = 14)

qsave(obj,"Rds/tmp.qs")
obj.list <- SplitObject(obj, split.by = "sample_id")

for (i in 1:length(obj.list)) {
  DefaultAssay(obj.list[[i]]) <- "drug"
  obj.list[[i]] <- FindTopFeatures(obj.list[[i]], min.cutoff = 5)
  obj.list[[i]] <- RunTFIDF(obj.list[[i]])
  obj.list[[i]] <- RunSVD(obj.list[[i]], reduction.name = "lsi", reduction.key = "LSI_")
}

integration.anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  anchor.features = rownames(obj),
  reduction = "rlsi",
  assay = rep("drug",length(obj.list)),
  dims = 2:30
)

integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = obj[["lsi.drug"]],
  new.reduction.name = "integrated_lsi.drug",
  dims.to.integrate = 1:30,
  k.weight = 35
)
integrated <- RunUMAP(integrated, reduction = "integrated_lsi.drug", dims = 2:30)
qsave(integrated,"Rds/integrated.drug.qs")
p2 <- DimPlot(integrated, group.by = "orig.ident") + DimPlot(integrated, group.by = "timepoint") + DimPlot(integrated)
ggsave("drug_inte_umap.png",p2,height = 7,width = 14)
obj[["integrated_lsi.drug"]] <- integrated[["integrated_lsi.drug"]]
qsave(obj,"Rds/obj.integrated.drug.qs")
p2 <- DimPlot(obj, group.by = "orig.ident",reduction = "umap.drug") + DimPlot(obj, group.by = "timepoint",reduction = "umap.drug") + DimPlot(obj,reduction = "umap.drug")
ggsave("drug_umap.png",p2,height = 7,width = 14)

# Graph clustering
obj <- RunUMAP(obj, reduction = "integrated_lsi.atac", dims = c(2:30), reduction.name = "umap.inte.atac")
obj <- RunUMAP(obj, reduction = "integrated_lsi.drug", dims = c(2:30), reduction.name = "umap.inte.drug")
DefaultAssay(obj)<- "ATAC"
lsi_list <- names(obj@reductions)[grepl("lsi", names(obj@reductions)) & !grepl("wnn", names(obj@reductions))]
lsi_list <- lsi_list[lsi_list != "lsi"]
for (i in lsi_list){
  prefix  <- gsub("\\.", "_", i)
  nn_name <- paste0(prefix, "_nn")
  snn_name <- paste0(prefix, "_snn")
  obj <- FindNeighbors(obj, reduction = i, dims = 2:30, graph.name = c(nn_name, snn_name), verbose = FALSE)
  for (res in c(0.1,0.2,0.3,0.4,0.5,0.7)){
    obj <- FindClusters(obj, resolution = res, verbose = FALSE,graph.name = snn_name)
    obj[[paste0("cluster_", i ,"_res.",res)]] <- obj$seurat_clusters
  }
}
qsave(obj,"Rds/obj.integrated.total0.qs")

# Cluster markers for manual cell-type annotation.
annotation_cluster <- "cluster_integrated_lsi.atac_res.0.1"
Idents(obj) <- annotation_cluster
p <- DimPlot(obj, reduction = "umap.inte.atac", group.by = annotation_cluster)
ggsave("clusters.pdf", p, width = 7, height = 5)

library(COSG)
marker_cosg <- cosg(
  obj,
  groups='all',
  assay='ATAC_genes',
  slot='data',
  mu=100,
  n_genes_user=100)
symbols_list <-  as.list(as.data.frame(apply(marker_cosg$names,2,head,50)))

write.csv(marker_cosg$names, "cluster_marker_genes.csv", row.names = FALSE)
qsave(obj, "Rds/obj_clustered.qs")
