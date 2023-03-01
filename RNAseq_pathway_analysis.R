##############################################################
############
######## --- . RNAseq based GSEA pathway analysis --- ########
############
##############################################################
library(broom) #CRAN
library(vsn) #BIOCONDUCTOR
library("DESeq2") #BIOCONDUCTOR
library(limma) #BIOCONDUCTOR
library(lattice) #CRAN
library(org.Hs.eg.db) #BIOCONDUCTOR
library(MASS) #CRAN
library(RColorBrewer) #BIOCONDUCTOR
library(AnnotationDbi) #BIOCONDUCTOR
library(GenomicRanges) #BIOCONDUCTOR
library(GenomicFeatures) #BIOCONDUCTOR
library(rtracklayer) #BIOCONDUCTOR
library(biomaRt) #BIOCONDUCTOR
library(glmnet) #CRAN         *
library(Hmisc) #CRAN
library(ConsensusClusterPlus) #BIOCONDUCTOR    *
library(ggplot2)#CRAN
library("edgeR") #BIOCONDUCTOR
library(affy) #BIOCONDUCTOR
library(Biobase)#BIOCONDUCTOR
library(WGCNA) #CRAN
library(impute) #BIOCONDUCTOR
library(GO.db) #BIOCONDUCTOR
library(gplots) #CRAN
library(viridis) #CRAN       
library(stats) #CRAN
library(robustbase) #CRAN
library(fgsea) #BIOCONDUCTOR
library(ComplexHeatmap) #BIOCONDUCTOR
library(clusterProfiler)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(igraph)
library(tidyverse)
library(readxl)

################################################
###  
### Selection gains larger than 50% of the chr8
### Selection gains that involve the q arm using centromere position as reference
###
###############################################
pred<- read.delim("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Prediction_Model/AMoritz/prediction_may_2020.txt", sep="\t", stringsAsFactors = F)
cnv<- read.delim("~/Google Drive/My Drive/Share/2022_final/CNV/CoMMpass_CN_data_collapsed.txt")
head(cnv)
#           sample Chrom    start      end major minor
# 1 MMRF_1016_1_BM     1   571000   720500     1     0
# 2 MMRF_1016_1_BM     1   730700  9450900     2     1
# 3 MMRF_1016_1_BM     1  9456900  9820400     1     0
# 4 MMRF_1016_1_BM     1  9826900 16199000     2     1
# 5 MMRF_1016_1_BM     1 16204900 16299400     1     0
# 6 MMRF_1016_1_BM     1 16304400 74541800     2     1
gain8<- cnv[cnv$Chrom==8 & (cnv$end-cnv$start)>73182011 & cnv$major>2,] ### select gains larger than 50% of the chr8
gain8<- gain8[gain8$end > 45600001, ] ### select gains that invovle the q arm using centromere position as ref
pred<- pred[pred$study =="MMRF",]
pred$gain8<- 0
pred$gain8[pred$sample %in% gsub("_1_BM","", gain8$sample)]<- 1
pred_mmrf<- pred[pred$study=="MMRF",]### select only MMRF
pred_mmrf<- pred_mmrf[pred_mmrf$sample %in% gsub("_1_BM","",cnv$sample), ]


####################
###
### RNA seq analysis
###
####################

# . . . - Raw count file
exp<-read.table("~/Documents/Project/GEP_CoMMpass/Before/data/MMRF_CoMMpass_IA12a_E74GTF_HtSeq_Gene_Counts.txt",
                header=T,
                sep="\t",
                row.names=1,
                as.is=T)
exp<-exp[,grep("1_BM",colnames(exp))]
colnames(exp)<-gsub("_1_BM","",colnames(exp))
exp[1:5,1:5]
#                 MMRF_1446 MMRF_1755 MMRF_1906 MMRF_1683 MMRF_2475
# ENSG00000000003         3        58         0       622        22
# ENSG00000000005         0         1         0         0         0
# ENSG00000000419      1230      3126      2244      1249      1948
# ENSG00000000457       626      1175      1549      1272       753
# ENSG00000000460       210       407       431       361       135


# . . . - Pheno data
pdata<-pred_mmrf
rownames(pdata)<- pdata$sample
pdata[1:2,]
#              sample age gender   ecog    ISS SCT_first_line time_SCT SCT_line pfs_code pfs_time os_code os_time KAR chemo BORT LEN THAL DARA ELO duration
# MMRF_1016 MMRF_1016  56   MALE ecog<2 ISS1-2              0       NA       NA        0      673       1     692   1     1    0   1    0    0   0      673
# MMRF_1021 MMRF_1021  54 FEMALE ecog<2 ISS1-2              1      127        1        1      532       0    2135   0     0    1   1    0    0   0      127
#           continuos_treat combo study ABCF1 ACTG1 ARID2 BCL7A BHLHE41 BRAF BTG1 CCND1 CYLD DIS3 DTX1 DUSP2 EGR1 FAM46C FGFR3 FUBP1 HIST1H1B HIST1H1D HIST1H1E
# MMRF_1016               1     3  MMRF     0     1     0     0       0    0    0     0    0    0    0     0    0      0     1     0        0        0        0
# MMRF_1021               0     2  MMRF     0     1     0     0       0    0    0     0    0    0    0     0    0      0     0     1        0        0        0
#           IGLL5 IRF1 IRF4 KLHL6 KMT2B KRAS LCE1D LTB MAX NFKB2 NFKBIA NRAS PABPC1 PIM1 POT1 PRDM1 PRKD2 PTPN11 RASA2 RB1 RPL10 RPL5 RPRD1B RPS3A SAMHD1 SETD2
# MMRF_1016     0    0    0     0     0    0     0   0   0     0      0    0      0    0    0     0     0      0     0   0     0    0      0     0      0     0
# MMRF_1021     1    0    0     0     0    1     0   0   0     0      0    0      0    0    0     0     0      0     0   0     0    0      0     0      0     0
#           SP140 TBC1D29 TCL1A TGDS TP53 TRAF2 TRAF3 ZNF292 del4p15 del12p13 del12q21 gain1q21 del8p22 gain6p22 del20p12 ampMYC delRB1 del6q delCDKN2C delTP53
# MMRF_1016     0       0     0    0    0     0     0      0       0        0        0        0       0        1        0      0      0     0         1       0
# MMRF_1021     0       0     0    0    0     0     0      0       0        0        0        1       0        0        0      0      1     0         0       0
#           delFAM46c delCYLD delTRAF3 delTRAF2 del13q34 del14q23 HRD t_CCND1 t_MMSET t_MAFB t_MAF t_CCND3 t_MYC gain8
# MMRF_1016         0       0        0        0        0        0   1       0       1      0     0       0     0     0
# MMRF_1021         1       0        1        0        1        1   0       0       1      0     0       0     0     0


# . . . - Arrange expression, pheno and gene data to create ExpressionSet object
load("~/Documents/Project/GEP_CoMMpass/Before/analysis/annot.RData")
fdata<-aggregate(fdata[,2:3,],by=list(ensembl_gene_id=fdata$ensembl_gene_id),paste,collapse=";")
rownames(fdata)<-fdata$ensembl_gene_id
fdata<-fdata[rownames(exp),]
identical(rownames(exp),rownames(fdata))
exp2<-exp[,colnames(exp)%in%rownames(pdata)]
pdata<-pdata[colnames(exp2),]
identical(rownames(pdata),colnames(exp2))
summary(pdata)
hyper_8<-pdata[,c("HRD", "gain8")]
rownames(hyper_8)<- pdata$sample
hyper_8.1<-data.frame(rep(1,length(hyper_8$gain8)),hyper_8$gain8,row.names = rownames(hyper_8))
colnames(hyper_8.1)<-c("intercept","gain8")
exp2<-exp[,colnames(exp)%in%rownames(hyper_8.1)]
hyper_8.1<-hyper_8.1[colnames(exp2),]
identical(rownames(hyper_8.1),colnames(exp2))
dataset_8<-ExpressionSet(assayData=as.matrix(exp2),
                         featureData = new("AnnotatedDataFrame",fdata),
                         phenoData = new("AnnotatedDataFrame",hyper_8.1))
collrow_8<-collapseRows(exprs(dataset_8),
                        rowID=rownames(exprs(dataset_8)),
                        rowGroup=fData(dataset_8)$ensembl_gene_id,
                        method="function",
                        methodFunction = colSums)
exp2<-collrow_8$datETcollapsed


# . . . - Import gene annotation file
annot<-read.table("~/Documents/Project/GEP_CoMMpass/Before/data/homo_sapins_ensembl.txt",
                  header=T,
                  sep="\t",
                  as.is=T,
                  quote="",
                  fill=T)
head(annot)
#    GeneID chromosome      GeneID_ens   Symbol map_location type_of_gene description
# 1 6775087         MT            <NA> 12S rRNA            -         rRNA      s-rRNA
# 2 8923213         MT            <NA> 12S rRNA            -         rRNA      s-rRNA
# 3 8923219         MT            <NA> 16S rRNA            -         rRNA      l-rRNA
# 4 6775085         MT            <NA> 16S rRNA            -         rRNA      l-rRNA
# 5      NA          2 ENSG00000271924  5S_rRNA         <NA>         rRNA        <NA>
# 6      NA          X ENSG00000212595  5S_rRNA         <NA>         rRNA        <NA>
colnames(annot)<-c("entrez","chr","ensembl_gene_id","hgnc_symbol","cytoband","transcript_type","gene_name")
fdata<-merge(fdata,annot,by="ensembl_gene_id",all.x=T)
fdata<-fdata[!fdata$gene_name%in%"hypophosphatemic bone disease",]
fdata<-fdata[!fdata$gene_name%in%"Methylation modifier for class I HLA",]
fdata<-fdata[!fdata$gene_name%in%"Miyoshi muscular dystrophy 2",]
fdata<-fdata[!fdata$gene_name%in%"transient erythroblastopenia of childhood",]
rownames(fdata)<-fdata$ensembl_gene_id
fdata<-fdata[rownames(exp2),]


# . . . - Creation ExpressionSet object using expression, pheno and gene data.
dataset_8<-ExpressionSet(assayData=exp2,
                         featureData = new("AnnotatedDataFrame",fdata),
                         phenoData = new("AnnotatedDataFrame",pdata))
# . . . - Keeping the genes with at least 10 raw counts in at least 95% of the samples.
isexp_8<-rowSums(exprs(dataset_8)>=10)>=ncol(dataset_8)/100*5
dataset_8.f<-dataset_8[isexp_8,]
# . . . - Creation DGEList object
y<-DGEList(counts = exprs(dataset_8.f),
           genes = fData(dataset_8.f),
           samples=pData(dataset_8.f))
# . . . - Applying the library size normalization
y<-calcNormFactors(y)
# The y object consists of 3 sub-objects:
y$counts[1:5,1:5] # Expression matrix
#                 MMRF_1446 MMRF_1755 MMRF_1906 MMRF_1683 MMRF_2475
# ENSG00000000003         3        58         0       622        22
# ENSG00000000419      1230      3126      2244      1249      1948
# ENSG00000000457       626      1175      1549      1272       753
# ENSG00000000460       210       407       431       361       135
# ENSG00000000938        20       643        71        29      2898
y$samples[1:2,] # Pheno data
#           group lib.size norm.factors    sample age gender   ecog    ISS SCT_first_line time_SCT SCT_line pfs_code pfs_time os_code os_time KAR chemo BORT
# MMRF_1446     1 39487444        1.405 MMRF_1446  68   MALE ecog<2 ISS1-2              1      346        1        1      438       0    1542   1     1    0
# MMRF_1755     1 70698301        1.233 MMRF_1755  78   MALE ecog<2 ISS1-2              1      206        1        0      256       0     256   0     0    1
#           LEN THAL DARA ELO duration continuos_treat combo study ABCF1 ACTG1 ARID2 BCL7A BHLHE41 BRAF BTG1 CCND1 CYLD DIS3 DTX1 DUSP2 EGR1 FAM46C FGFR3
# MMRF_1446   1    0    0   0      346               0     3  MMRF     0     0     0     0       0    0    0     0    0    0    0     0    0      1     0
# MMRF_1755   1    0    0   0      206               0     2  MMRF     0     0     0     0       0    0    0     0    0    0    0     0    0      0     0
#           FUBP1 HIST1H1B HIST1H1D HIST1H1E IGLL5 IRF1 IRF4 KLHL6 KMT2B KRAS LCE1D LTB MAX NFKB2 NFKBIA NRAS PABPC1 PIM1 POT1 PRDM1 PRKD2 PTPN11 RASA2 RB1
# MMRF_1446     0        0        0        0     0    0    0     0     0    0     0   0   0     0      0    0      0    0    0     0     0      0     0   0
# MMRF_1755     0        0        0        0     0    0    0     0     0    0     0   0   0     0      0    0      0    0    0     0     0      0     0   0
#           RPL10 RPL5 RPRD1B RPS3A SAMHD1 SETD2 SP140 TBC1D29 TCL1A TGDS TP53 TRAF2 TRAF3 ZNF292 del4p15 del12p13 del12q21 gain1q21 del8p22 gain6p22
# MMRF_1446     0    0      0     0      0     0     0       0     0    0    0     0     0      0       0        0        0        0       0        0
# MMRF_1755     0    0      0     0      0     0     0       0     0    0    1     0     0      0       0        0        0        0       1        0
#           del20p12 ampMYC delRB1 del6q delCDKN2C delTP53 delFAM46c delCYLD delTRAF3 delTRAF2 del13q34 del14q23 HRD t_CCND1 t_MMSET t_MAFB t_MAF t_CCND3
# MMRF_1446        0      0      0     0         0       0         0       0        0        0        0        1   0       1       0      0     0       0
# MMRF_1755        1      0      0     1         0       1         0       1        0        0        0        0   1       1       0      0     0       0
#     t_MYC gain8
# MMRF_1446     0     0
# MMRF_1755     0     0
head(y$genes) # Gene annotation data
#                 ensembl_gene_id hgnc_symbol.x    hgnc_id entrez chr hgnc_symbol.y cytoband transcript_type
# ENSG00000000003 ENSG00000000003        TSPAN6 HGNC:11858   7105   X        TSPAN6   Xq22.1  protein_coding
# ENSG00000000419 ENSG00000000419          DPM1  HGNC:3005   8813  20          DPM1 20q13.13  protein_coding
# ENSG00000000457 ENSG00000000457         SCYL3 HGNC:19285  57147   1         SCYL3   1q24.2  protein_coding
# ENSG00000000460 ENSG00000000460      C1orf112 HGNC:25565  55732   1      C1orf112   1q24.2  protein_coding
# ENSG00000000938 ENSG00000000938           FGR  HGNC:3697   2268   1           FGR   1p35.3  protein_coding
# ENSG00000000971 ENSG00000000971           CFH  HGNC:4883   3075   1           CFH   1q31.3  protein_coding
#                                                                   gene_name
# ENSG00000000003                                               tetraspanin 6
# ENSG00000000419 dolichyl-phosphate mannosyltransferase subunit 1, catalytic
# ENSG00000000457                                    SCY1 like pseudokinase 3
# ENSG00000000460                         chromosome 1 open reading frame 112
# ENSG00000000938              FGR proto-oncogene, Src family tyrosine kinase
# ENSG00000000971                                         complement factor H


# . . . - Design matrix generation 
design_1<-model.matrix(~0+as.factor(gain8),
                       data = y$samples)
colnames(design_1)[c(1:2)]<-c("No_tri_chr8","tri_chr8")
head(design_1)
#           No_tri_chr8 tri_chr8
# MMRF_1446           1        0
# MMRF_1755           1        0
# MMRF_1906           1        0
# MMRF_1683           1        0
# MMRF_2475           1        0
# MMRF_1519           1        0

# . . . Voom/Limma pipeline
#       voom  -> transform RNAseq data ready for linear modelling
#       lmFit -> linear model for each gene given a series of arrays
v<-voom(y,design=design_1,plot=F)
fit<-lmFit(v,design_1)
#       Samples WITH tri_chr8  VS  Samples WITHOUT tri_chr8
contrast.matrix=makeContrasts(tri_chr8-No_tri_chr8,
                              levels=design_1)
fit=contrasts.fit(fit, contrast.matrix)
efit_8=eBayes(fit)
#       Extract a table of the top-ranked genes from a linear model fit
q1<-topTable(efit_8,coef=1,number=nrow(v),adjust="BH",p.value = 0.05)
dim(q1)
q2<-topTable(efit_8,coef=1,number=nrow(v),adjust="BH")
dim(q2)
q1<-q2


# . . . - Pathway Analysis - GSEA
setlist<-read.table("~/Desktop/Tools/Reference/gsea/h.all.v6.2.symbols.gmt_MOD_CART_paper.txt",sep="\t",fill=T,header=F,as.is=T)
gs<-vector("list",length=nrow(setlist))
names(gs)<-setlist[,1]
setlist<-setlist[,-c(1:2)]
for(i in 1:nrow(setlist)){
  temp<-as.vector(as.matrix(setlist[i,]))
  temp<-temp[temp!=""]
  gs[[i]]<-temp
}

res_8<-topTable(efit_8,coef=1,number=nrow(v),adjust="BH")
rnk_8<-res_8$t
names(rnk_8)<-res_8$hgnc_symbol.x
fgseaRes_8 <- fgseaMultilevel(pathways = gs, 
                              stats = rnk_8,
                              minSize=15,
                              maxSize=500)
fgseaRes_8<-fgseaRes_8[order(fgseaRes_8$pval),]
fwrite(fgseaRes_8,file="~/Google Drive/My Drive/Bachisio/Documents/Project/Team/MauraF/DARA-KRD_GeneExpression/trisomie/chr8/gsea_analysis_genesets_MODchr8_022023.txt",sep="\t",sep2=c("", ";", ""))



res_8<-read.table("~/Google Drive/My Drive/Bachisio/Documents/Project/Team/MauraF/DARA-KRD_GeneExpression/trisomie/chr8/gsea_analysis_genesets_MODchr8_022023.txt",sep="\t",fill=T,header=T,as.is=T)
res_8<-res_8[order(res_8$pathway),]

nes__8<-matrix(0,nrow=52,ncol=2)
fdr__8<-matrix(0,nrow=52,ncol=2)
nes__8[,1]<-res_8$NES
fdr__8[,1]<-res_8$padj
nes__8[,2]<-0
fdr__8[,2]<-1

rownames(nes__8)<-res_8$pathway
rownames(fdr__8)<-res_8$pathway
colnames(nes__8)<-c("Trisomie chr8","offset")
colnames(fdr__8)<-c("Trisomie chr8","offset")

heatmap.2(nes__8,
          trace="none",
          dendrogram="none",
          col=brewer.pal(11,"RdYlBu")[11:1],
          density.info="none",
          key.xlab = "NES",
          mar=c(2,15),
          labRow=gsub("HALLMARK_","",rownames(nes__8)),
          colsep=0:ncol(nes__8),
          rowsep=0:nrow(nes__8),
          sepcolor="grey30",
          sepwidth=c(0.001,0.001),
          cellnote=ifelse(fdr__8<0.05,"X",""),
          notecex=0.6,
          notecol="black",
          cexCol = 1)
legend("topright",
       legend="X: FDR >= 0.05",
       bty="n",
       inset=c(0.3,0))
