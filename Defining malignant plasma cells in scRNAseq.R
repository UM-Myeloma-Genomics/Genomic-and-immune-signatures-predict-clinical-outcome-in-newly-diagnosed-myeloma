# Define malignant plasma cell population in scRNAseq
# David Coffey
# February 22, 2023

library(Seurat)
library(tidyverse)
library(infercnv)
library(data.table)

# STEP 1: Run Seurat to QC, integrate, and annotate scRNAseq datasets

# Read in GEO gene matrix files
paths = list.files("GEO_Gene_Matrix_Directory/", full.names = TRUE, all.files = FALSE, recursive = FALSE, include.dirs = FALSE)
names = gsub(basename(paths), pattern = ".txt", replacement = "")

# Convert gene matrix files into a list of Seurat objects
sce = list()
i = 1
for(i in 1:length(names)){
  mtx = read.table(paths[i], row.names = 1)
  mtx = CreateSeuratObject(counts = mtx, project = names[i])
  sce = c(sce, list(mtx))
}

# Calculate nCount_RNA and nFeature_RNA
sce = lapply(X = sce, FUN = function(x) {
  calcn = as.data.frame(x = Seurat:::CalcN(object = x))
  colnames(x = calcn) = paste(colnames(x = calcn), "RNA", sep = '_')
  x = AddMetaData(object = x, metadata = calcn)
})

# Calculate percent mitochondrial genes if the sce contains genes matching the regular expression "^MT-"
sce = lapply(X = sce, FUN = function(x) {
  x = PercentageFeatureSet(object = x, pattern = '^MT-', col.name = "Percent_mt", assay = "RNA")
})

# Remove cells with < 500 UMI, < 250 genes, or > 20% mitochondrial genes
sce.filtered = lapply(X = sce, FUN = function(x) {
  cells.use = x[["nCount_RNA", drop = TRUE]] > 500 | x[["nFeature_RNA", drop = TRUE]] > 250 | x[["Percent_mt", drop = TRUE]] < 20
  x = x[, cells.use]
})

# Normalize
sce.normalized = lapply(X = sce.filtered, FUN = function(x) {
  x = SCTransform(x, verbose = FALSE, conserve.memory = TRUE)
})

# Feature selection
features = SelectIntegrationFeatures(object.list = sce.normalized)

# Run PCA
sce.normalized = lapply(X = sce.normalized, FUN = function(x) {
  x = RunPCA(x, features = features, verbose = FALSE)
})

# Integrate data sets selecting a representative reference sample from each of the 4 sample types: 
# 10 GSM4891350, RPMI-8226 cell line
# 16 GSM4891356, Baseline
# 38 GSM4891378, Cycle 4
# 152 GSM4891492, Cycle 10

anchors.integrate = FindIntegrationAnchors(object.list = sce.normalized, reference = c(10, 16, 38, 152), reduction = "rpca", dims = 1:50)
integrated = IntegrateData(anchorset = anchors.integrate, dims = 1:50)

# Find anchors between reference and query using a precomputed supervised PCA (spca) transformation (Azimuth reference datasets available here: https://azimuth.hubmapconsortium.org/references/)
reference = LoadH5Seurat("pbmc_multimodal.h5seurat")
integrated.normalized = SCTransform(integrated, verbose = FALSE, conserve.memory = TRUE)

anchors.mapping = FindTransferAnchors(
  reference = reference,
  query = integrated.normalized,
  recompute.residuals = FALSE,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

# MapQuery is a wrapper around three functions: TransferData, IntegrateEmbeddings, and ProjectUMAP. 
# TransferData is used to transfer cell type labels and impute the ADT values. 
# IntegrateEmbeddings and ProjectUMAP are used to project the query data onto the UMAP structure of the reference.
integrated.mapped = MapQuery(
  anchorset = anchors.mapping,
  query = integrated.normalized,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

# Add metadata
meta.data = integrated.mapped@meta.data
meta.data$barcode = rownames(meta.data)
meta.data$Sample.Name = gsub(meta.data$orig.ident, pattern = "_.*", replacement = "")
meta.data$azimuth = meta.data$predicted.celltype.l2

sample.metadata = read.csv("SraRunTable.csv")
sample.metadata = unique(sample.metadata[,c("Replicate_ID", "Sample.Name", "source_name", "Time_point", "disease_state", "Organ", "selection_marker", "Cohort", "Cell_Line", "Genotype")])

meta.data.new = right_join(sample.metadata, meta.data)
rownames(meta.data.new) = meta.data.new$barcode

integrated.mapped@meta.data = meta.data.new
save(integrated.mapped, file = "integrated.mapped.Rda")

# STEP 2: Run InferCNV to estimate CNV

# Load integrated and mapped seurat object
load("integrated.mapped.Rda")

# Subset seurat object into B cell and non B cell populations
integrated.mapped = integrated.mapped[,!(integrated.mapped@meta.data$labels %in% c("Pro-B_cell_CD34+","Pre-B_cell_CD34-"))]
integrated.mapped@meta.data$CellGroups = ifelse(integrated.mapped@meta.data$labels == "B_cell", "B cell", "Non B cells")

# Run InferCNV on each sample
samples = unique(integrated.mapped@meta.data$Sample.Name)
for(i in 1:length(samples)){
  integrated.mapped.sample = integrated.mapped[,integrated.mapped@meta.data$Sample.Name == samples[i]]
  
  # Import gene order file (data file containing the positions of each gene along each chromosome in the genome)
  gene.order = data.frame(fread("InferCNV gene order file by gene symbol.txt"))
  
  if(sum(integrated.mapped.sample@meta.data$CellGroups == "B cell") > 2){
    
    # Create destination directory
    system(paste("mkdir InferCNV/", samples[i], sep = ""))
    
    # Make invercnv object
    infercnv = CreateInfercnvObject(raw_counts_matrix=integrated.mapped.sample@assays$RNA@counts,
                                    annotations_file=data.frame(row.names = integrated.mapped.sample@meta.data$barcode, integrated.mapped.sample@meta.data$CellGroups),
                                    gene_order_file=gene.order,
                                    ref_group_names="Non B cells")
    
    # Run infercnv
    infercnv = infercnv::run(infercnv, 
                             cutoff = 0.1, 
                             out_dir = paste("InferCNV/", samples[i], sep = ""), 
                             cluster_by_groups = FALSE, 
                             denoise = TRUE, 
                             HMM = TRUE, 
                             analysis_mode = "subclusters", 
                             tumor_subcluster_partition_method = "random_trees")
    
    # Add results to seurat object
    infercnv_obj = readRDS(paste(paste("InferCNV/", samples[i], sep = ""), "run.final.infercnv_obj", sep=.Platform$file.sep))
    integrated.mapped.sample = integrated.mapped.sample[, colnames(integrated.mapped.sample@assays$RNA) %in% colnames(infercnv_obj@expr.data)]
    integrated.mapped.sample = infercnv::add_to_seurat(infercnv_output_path = paste("InferCNV/", samples[i], sep = ""), seurat_obj = integrated.mapped.sample, top_n = 10)
    save(integrated.mapped.sample, file = paste("InferCNV/", samples[i], "/Seurat.Rda", sep = ""))
  }
}

# STEP 3: Merge CNV output

knitr::opts_chunk$set(eval = FALSE)
# Import SRA run table to create vector of samples names for for loop
samples = read.csv("SraRunTable.csv")
samples = unique(samples$Sample.Name)
metadata.merged = data.frame()
i = 1
pb = txtProgressBar(min = 0, max = length(samples), initial = 0, style = 3) # Progress bar
for(i in 1:length(samples)){
  setTxtProgressBar(pb, i)
  if(is.na(file.size(paste("InferCNV/", samples[i], "/Seurat.Rda", sep = "")))){
    print(paste("Skipping", samples[i]))
  } else {
    load(paste("InferCNV/", samples[i], "/Seurat.Rda", sep = ""))
    metadata = integrated.mapped.sample@meta.data
    metadata.merged = merge(metadata, metadata.merged, all = TRUE)
  }
  close(pb)
}
save(metadata.merged, file = "Merged infercnv.Rda")

# STEP 4: Define clones from CNV output

# Load merged InferCNV
load("Merged infercnv.Rda")

# Load integrated and mapped seurat object
load("integrated.mapped.Rda")
metadata.merged = metadata.merged[metadata.merged$Organ == "Bone marrow",]

# Define clone ID according to CNVs
cnv = metadata.merged[,grep(names(metadata.merged), pattern = "has_loss|has_dupli", value = TRUE)]
rownames(cnv) = metadata.merged$barcodes
cnv = ifelse(cnv == TRUE, 1, 0)
cnv = t(cnv)

index = split(colnames(cnv), apply(cnv, 2, paste, collapse = ""))
names(index) = paste("Clone", seq(1:length(index)))
clones = stack(index)
names(clones) = c("barcodes", "CloneID")

# Create clone summary table
cnv.barcodes = metadata.merged[,c("barcodes", "Replicate_ID", "predicted.celltype.l1", grep(names(metadata.merged), pattern = "has_loss|has_dupli", value = TRUE))]
cnv.barcodes = merge(clones, cnv.barcodes)

# Create sample specific clone ID
samples = unique(cnv.barcodes$Replicate_ID)
i = 1
j = 1
sample.clones = data.frame()
for(i in 1:length(samples)) {
  s = cnv.barcodes[cnv.barcodes$Replicate_ID == samples[i],c("Replicate_ID", "CloneID")]
  c = unique(s$CloneID)
  c = c[!(is.na(c))]
  if(length(c) > 1) {
    for(j in 1:length(c)){
      s.c = s[s$CloneID == c[j] & !(is.na(s$CloneID)),]
      s.c$SampleCloneID = paste("Clone", j - 1)
      sample.clones = rbind(sample.clones, s.c)
    }
  }
}

# Summary table of clone count and frequency.  Primary clones are defined as the most frequent cells with the same CNV.  All others are subclones.
cnv.barcodes.summary = merge(unique(sample.clones), cnv.barcodes)
cnv.barcodes.summary = cnv.barcodes.summary[cnv.barcodes.summary$predicted.celltype.l1 == "B",]
cnv.barcodes.summary[cnv.barcodes.summary == TRUE] = 1
cnv.barcodes.summary[cnv.barcodes.summary == FALSE] = 0
cnv.barcodes.summary = cnv.barcodes.summary[rowSums(cnv.barcodes.summary[,c(6:ncol(cnv.barcodes.summary))]) != 0,]
clone.count = aggregate(data = cnv.barcodes.summary, barcodes~SampleCloneID+Replicate_ID, length)
names(clone.count)[3] = "CloneCount"
clone.max = slice_max(group_by(clone.count, Replicate_ID), CloneCount, n = 1)
clone.max$CloneAbundance = rep("Primary")
cell.count = data.frame(table(Replicate_ID = metadata.merged$Replicate_ID))
names(cell.count)[2] = "CellCount"
cnv.barcodes.summary = unique(cnv.barcodes.summary[,-c(2,4,5)])
cnv.barcodes.summary = merge(cnv.barcodes.summary, clone.count)
cnv.barcodes.summary = merge(clone.max, cnv.barcodes.summary, all = TRUE)
cnv.barcodes.summary = merge(cell.count, cnv.barcodes.summary)
cnv.barcodes.summary$CloneAbundance = ifelse(is.na(cnv.barcodes.summary$CloneAbundance), "Subclone", "Primary")
cnv.barcodes.summary$CNVCount = rowSums(cnv.barcodes.summary[,6:ncol(cnv.barcodes.summary)])
cnv.barcodes.summary$CloneFrequency = cnv.barcodes.summary$CloneCount/cnv.barcodes.summary$CellCount * 100
write.csv(cnv.barcodes.summary, file = "Table of clone characteristics.csv", row.names = FALSE)

# Save primary clone barcodes
primary.clone.barcodes = merge(cnv.barcodes, cnv.barcodes.summary[cnv.barcodes.summary$CloneAbundance == "Primary", c("Replicate_ID", "SampleCloneID")])
write.csv(primary.clone.barcodes[,c("Replicate_ID", "Barcode")], file = "Primary clone barcodes.csv", row.names = FALSE)

# STEP 5: Compute average expression in CNV clones

# Import barcodes that define malignant clones
primary.clone.barcodes = read.csv("Primary clone barcodes.csv")

# Subset seurat object with malignant clones
selected = integrated.mapped[,colnames(integrated.mapped)%in% primary.clone.barcodes$Barcode]

# Normalize counts
normalized = NormalizeData(object = selected)

# Compute averaged expression by sample
selected.average.exp = AverageExpression(normalized, group.by = "Replicate_ID", return.seurat = TRUE, slot = "data")
write.csv(selected.average.exp@assays$RNA@data, file = "Normalized averaged expression of malignant clones.csv")


