# Genomic and immune signatures predict clinical outcome in newly diagnosed multiple myeloma treated with immunotherapy regimens
###### Francesco Maura1*, Eileen M. Boyle2*, David Coffey1*, Kylee Maclachlan3*, Dylan Gagler2, Benjamin Diamond1, Hussein Ghamlouch2, Patrick Blaney2, Bachisio Ziccheddu1, Anthony Cirrincione1, Monika Chojnacka1, Yubao Wang2, Ariel Siegel2, James E. Hoffman1, Dickran Kazandjian1, Hani Hassoun3, Emily Guzman4, Sham Mailankody3, Urvi Shah3, Carlyn Tan3, Malin Hultcrantz3, Michael Scordo3, Gunjan Shah3, Heather Landau3, David J. Chung3, Sergio Giralt3, Yanming Zhang5, Arnaldo Arbini2, Ahmet Dogan6, Alexander M Lesokhin3, Faith E Davies2, Saad Usmani3, Neha Korde3, Gareth J Morgan2,7#, Ola Landgren1#

Contains code for defining malignant plasma cells using copy number variation in single cell RNA sequencing and performing the pathway enrichment analysis.

### Defining malignant plasma cells in single-cell RNA sequencing
This R scripts demonstrates how we performed a single-cell RNA sequencing analysis on CD138+/CD38+ FACS sorted plasma cells from a previously published study of 41 patients with multiple myeloma who failed to respond to a bortezomid-containing induction regimen and received daratumumab, carfilzomib, lenalidomide, and dexamethasone (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA675855).  To define the malignant plasma cell population, we used inferCNV to call copy number variants (CNV) in B cells using non-B cells as the germline reference (https://github.com/broadinstitute/inferCNV). For this analysis, we used the inferCNV “subclusters” mode and the “random_trees” method for partitioning the hierarchical clustering tree. The B cells were then grouped according to like CNVs and the dominate clone, which comprised the greatest number of B cells per sample, was chosen for downstream analyses. The average expression per sample was computed using the AverageExpression function in Seurat.

### RNAseq pathway analysis
This R code describes how we implement the RNAseq pathway analysis to validate the prognostic and biological impact of chr8 large chromosomal events. We implemeted the RNAseq analysis on the CoMMpass dataset (https://research.themmrf.org). The raw count data was filtered to remove genes with less than 10 reads in greater than 95% of samples. The trimmed mean of M-values (TMM) normalization was applied. This method estimates a scale factor used to reduce technical bias between samples with different library sizes. The voom transformation was performed to convert row counts in log2-counts per million (log2-CPM) and calculated the respective observation-level weights to be used in differential expression analysis. P values were corrected for multiple testing using Benjamini-Hochberg false discovery rate method. The Gen Set Enrichment Analysis (GSEA) was performed using the `fgsea` R package. The H Hallmark gene sets collection, retrieved from MSigDb database v7.4, was enriched with two INF signatures, ISG.RS and IFNG.GS, previously described to be associated with response to immunotherapy. Genes were ranked using the statistic derived from differential expression analysis with voom/limma pipeline.

