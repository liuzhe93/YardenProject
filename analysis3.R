
###Reference-based decomposition
#/Users/liuzhe/Desktop/UCSF/milk/scRNA/ref_based
setwd("/wynton/home/ahituv/liuzhe/projects/yarden_milk/scRNA_decom")
remove(list=ls())
library(Biobase)
library(BisqueRNA)
bulk_RNA<-read.table("Partek_Supply_RNAseq_Normalization_Normalized_counts.txt", header = T, sep = "\t")
bulk.matrix<-subset(bulk_RNA, select = c("gene_name", "ID1016.V1_S0", "ID1020.V3_S148", "ID1020.V2_S146", "ID1018.V2_S0",
                                         "ID1020.V5_S149", "ID1028.V1_S0", "ID1020.V1_S145", "ID1023.V2_S140", "ID1026.V1_S137",
                                         "VPL01.131_S150", "ID1016.V2_S0", "ID1024.V3_S139", "ID1025.V1_S138", "ID1029.V1_S0"))
bulk.matrix<-cbind(bulk.matrix$gene_name, bulk.matrix$ID1016.V1_S0, bulk.matrix$ID1020.V1_S145, bulk.matrix$ID1023.V2_S140,
                   bulk.matrix$ID1025.V1_S138, bulk.matrix$ID1026.V1_S137, bulk.matrix$ID1028.V1_S0,
                   bulk.matrix$ID1029.V1_S0)
colnames(bulk.matrix)<-c("gene_name", "1016-v1", "1020-v1", "1023-v2", "1025-v1", "1026-v1", "1028-v1", "1029-v1")
bulk.matrix<-as.data.frame(bulk.matrix)
rownames(bulk.matrix)<-bulk.matrix$gene_name
bulk.matrix$gene_name<-NULL
bulk.matrix<-as.matrix(bulk.matrix)
gene.list<-rownames(bulk.matrix)
bulk.matrix <- apply(bulk.matrix,2,as.numeric) 
rownames(bulk.matrix)<-gene.list
bulk.eset <- Biobase::ExpressionSet(assayData = bulk.matrix)
bulk.eset
#ExpressionSet (storageMode: lockedEnvironment)
#assayData: 21726 features, 7 samples 
#element names: exprs 
#protocolData: none
#phenoData: none
#featureData: none
#experimentData: use 'experimentData(object)'
#Annotation:

bulk.eset@assayData$exprs[1:4,1:7]
#1016-v1  1020-v1  1023-v2   1025-v1  1026-v1  1028-v1   1029-v1
#A1BG     1.519600 1.800380 1.026200 0.0600146 0.978740 0.116584  0.912329
#A1BG-AS1 0.226278 0.511430 0.377537 0.6591600 1.383000 0.192899  0.261494
#A2M      6.228290 0.115011 6.402980 0.1199290 1.823930 0.402569 18.405600
#A2M-AS1  0.384029 0.153315 0.205321 0.4794170 0.222518 0.121807  0.214742

library(Seurat)
library(Matrix)
raw_counts <- read.csv("supply_raw_counts.csv")
rownames(raw_counts)<-raw_counts$X
raw_counts$X<-NULL
raw_counts_t<-t(raw_counts)
raw_counts_t[1:4, 1:4]
#           1005_AAACCCAGTGATAGTA-1 1005_AAACGAAAGTTTCGAC-1 1005_AAAGAACGTTACGATC-1 1005_AAATGGACATCCAATG-1
#AL627309.1                       0                       0                       0                       0
#AL627309.5                       0                       0                       0                       0
#AL627309.4                       0                       0                       0                       0
#AP006222.2                       0                       0                       0                       0
celldata<-read.csv("supply_cell_labels.csv")
rownames(celldata)<-celldata$X
celldata$X<-NULL
cellmetadata <- celldata
head(cellmetadata)
library(dplyr)
library(tidyr)
cellmetadata$sampleid<-rownames(cellmetadata)
cellmetadata %>% separate(sampleid, c("sample", "barcode"), "_") ->cellmetadata
table(cellmetadata$sample)
#1005    1012 1016-v1 1020-v1 1020-v4 1023-v2 1024-v1 1025-v1 1026-v1 1028-v1 1029-v1 1030-v2 1031-v1 1033-v1 
# 491    1151    3256    2389     189    1812     352    5375    1309     892    1163     224    4295     719 

milk_scRNA <- CreateSeuratObject(counts = raw_counts_t, project = "milksupply")
milk_scRNA
#An object of class Seurat 
#22949 features across 23617 samples within 1 assay 
#Active assay: RNA (22949 features, 0 variable features)
milk_scRNA[["percent.mt"]]<-PercentageFeatureSet(milk_scRNA, pattern = "^MT-")
load_number<- dim(milk_scRNA)
library(dplyr)
milk_scRNA<-milk_scRNA %>% NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dim = 1:30)

milk_scRNA<-AddMetaData(milk_scRNA, metadata = cellmetadata)
head(cellmetadata)
#boxplot(as.numeric(cellmetadata$percent.mito))
Idents(milk_scRNA)<-milk_scRNA$general.clusters

levels(milk_scRNA)
#[1] "LC2"    "LC1"    "Mac/DC" "T cell" "NK"     "B cell" "Plasma"

milk_scRNA_subset<-subset(x = milk_scRNA, subset = sample %in% c("1016-v1", "1020-v1", "1023-v2", "1025-v1",
                                                                 "1026-v1", "1028-v1", "1029-v1"))
milk_scRNA_subset
#An object of class Seurat 
#22949 features across 16196 samples within 1 assay 
#Active assay: RNA (22949 features, 2000 variable features)
#2 dimensional reductions calculated: pca, umap

#milk_scRNA_subset<-milk_scRNA[,milk_scRNA@meta.data$sample %in% c("1016-v1", "1020-v1", "1023-v2", "1025-v1",
#                                                                  "1026-v1", "1028-v1", "1029-v1")]

#milk_scRNA_1016_v1<-subset(x = milk_scRNA, subset = sample == "1016-v1")
#milk_scRNA_1020_v1<-subset(x = milk_scRNA, subset = sample == "1020-v1")
#milk_scRNA_1023_v2<-subset(x = milk_scRNA, subset = sample == "1023-v2")
#milk_scRNA_1025_v1<-subset(x = milk_scRNA, subset = sample == "1025-v1")
#milk_scRNA_1026_v1<-subset(x = milk_scRNA, subset = sample == "1026-v1")
#milk_scRNA_1028_v1<-subset(x = milk_scRNA, subset = sample == "1028-v1")
#milk_scRNA_1029_v1<-subset(x = milk_scRNA, subset = sample == "1029-v1")

Idents(milk_scRNA_subset)<-milk_scRNA_subset$orig.ident
levels(milk_scRNA_subset)

Idents(milk_scRNA_subset)<-milk_scRNA_subset$general.clusters
levels(milk_scRNA_subset)

table(sc.eset[["SubjectName"]])
table(sc.eset[["SubjectName"]])
table(sc.eset[["cellType"]])
sum(is.na(sc.eset[["SubjectName"]]))
sum(is.na(sc.eset[["cellType"]]))


sc.eset <- BisqueRNA::SeuratToExpressionSet(milk_scRNA_subset, delimiter="_", position=1, version="v3")


#Split sample names by "_" and checked position 1. Found 7 individuals.
#Example: "1016-v1_AAACCCAAGTCCTACA-1" corresponds to individual "1016-v1".
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=F)
ref.based.estimates <- res$bulk.props
knitr::kable(ref.based.estimates, digits=2)
write.csv(ref.based.estimates, "percentage_individuals.csv", quote= F)


setwd("/Users/liuzhe/Desktop/UCSF/milk/scRNA/ref_based/")
remove(list=ls())

ref.based.estimates<-read.csv("percentage_individuals.csv", row.names = 1)

mydata<-ref.based.estimates
summary(mydata)
mydata<-as.data.frame(mydata)

mydata$sampleId<-rownames(mydata)
library("dplyr")
library(reshape2)
library(ggplot2)
data_long <- melt(mydata, id = c("sampleId"), measure = c("X1016.v1",  "X1020.v1",   "X1023.v2",
                                                          "X1025.v1",  "X1026.v1",   "X1028.v1",
                                                          "X1029.v1"))
pdf("barplot.pdf", width = 15, height = 6)
ggplot(data = data_long, mapping = aes(x = factor(variable), y = value, fill = sampleId)) + 
  geom_bar(stat = 'identity', position = 'dodge')
dev.off()
pdf("barplot_stacked.pdf", width = 15, height = 6)
ggplot(data = data_long, mapping = aes(x = factor(variable), y = value, fill = sampleId)) + geom_bar(stat = 'identity', position = 'stack')
dev.off()
#guides(fill = guide_legend(reverse = TRUE))


library(hdf5r)
library(Seurat) 
sample_1005<-Read10X_h5("Partek_scRNA_Feb_22-cellranger_new_exported (1)/filtered_feature_bc_matrix1005.h5")
seurat_1005<-CreateSeuratObject(sample_1005, project = "sample_1005")
seurat_1005
#An object of class Seurat 
#21331 features across 491 samples within 1 assay 
#Active assay: RNA (21331 features, 0 variable features)

sample_1012<-Read10X_h5("Partek_scRNA_Feb_22-cellranger_new_exported (1)/filtered_feature_bc_matrix1012.h5")
sample_1012<-CreateSeuratObject(sample_1012, project = "sample_1012")
sample_1012
#An object of class Seurat 
#21331 features across 1151 samples within 1 assay 
#Active assay: RNA (21331 features, 0 variable features)

sample_1016v1<-Read10X_h5("Partek_scRNA_Feb_22-cellranger_new_exported (1)/filtered_feature_bc_matrix1016v1.h5")
sample_1016v1<-CreateSeuratObject(sample_1016v1, project = "sample_1016v1")
sample_1016v1
#An object of class Seurat 
#21331 features across 3256 samples within 1 assay 
#Active assay: RNA (21331 features, 0 variable features)

sample_1020v1<-Read10X_h5("Partek_scRNA_Feb_22-cellranger_new_exported (1)/filtered_feature_bc_matrix1020v1.h5")
sample_1020v1<-CreateSeuratObject(sample_1005, project = "sample_1020v1")
sample_1020v1
#An object of class Seurat 
#21331 features across 491 samples within 1 assay 
#Active assay: RNA (21331 features, 0 variable features)

sample_1020v4<-Read10X_h5("Partek_scRNA_Feb_22-cellranger_new_exported (1)/filtered_feature_bc_matrix1020v4.h5")
sample_1020v4<-CreateSeuratObject(sample_1005, project = "sample_1020v4")
sample_1020v4
#An object of class Seurat 
#21331 features across 491 samples within 1 assay 
#Active assay: RNA (21331 features, 0 variable features)

sample_1023v2<-Read10X_h5("Partek_scRNA_Feb_22-cellranger_new_exported (1)/filtered_feature_bc_matrix1023v2.h5")
sample_1023v2<-CreateSeuratObject(sample_1005, project = "sample_1023v2")
sample_1023v2
#An object of class Seurat 
#21331 features across 491 samples within 1 assay 
#Active assay: RNA (21331 features, 0 variable features)

sample_1024v1<-Read10X_h5("Partek_scRNA_Feb_22-cellranger_new_exported (1)/filtered_feature_bc_matrix1024v1.h5")
sample_1024v1<-CreateSeuratObject(sample_1005, project = "sample_1024v1")
sample_1024v1
#An object of class Seurat 
#21331 features across 491 samples within 1 assay 
#Active assay: RNA (21331 features, 0 variable features)

sample_1025v1<-Read10X_h5("Partek_scRNA_Feb_22-cellranger_new_exported (1)/filtered_feature_bc_matrix1025v1.h5")
sample_1025v1<-CreateSeuratObject(sample_1005, project = "sample_1025v1")
sample_1025v1
#object of class Seurat 
#21331 features across 491 samples within 1 assay 
#Active assay: RNA (21331 features, 0 variable features)

sample_1026v1<-Read10X_h5("Partek_scRNA_Feb_22-cellranger_new_exported (1)/filtered_feature_bc_matrix1026v1.h5")
sample_1026v1<-CreateSeuratObject(sample_1005, project = "sample_1026v1")
sample_1026v1
#An object of class Seurat 
#21331 features across 491 samples within 1 assay 
#Active assay: RNA (21331 features, 0 variable features)

sample_1028v1<-Read10X_h5("Partek_scRNA_Feb_22-cellranger_new_exported (1)/filtered_feature_bc_matrix1028v1.h5")
sample_1028v1<-CreateSeuratObject(sample_1005, project = "sample_1028v1")
sample_1028v1
#An object of class Seurat 
#21331 features across 491 samples within 1 assay 
#Active assay: RNA (21331 features, 0 variable features)

sample_1029v1<-Read10X_h5("Partek_scRNA_Feb_22-cellranger_new_exported (1)/filtered_feature_bc_matrix1029v1.h5")
sample_1029v1<-CreateSeuratObject(sample_1005, project = "sample_1029v1")
sample_1029v1
#An object of class Seurat 
#21331 features across 491 samples within 1 assay 
#Active assay: RNA (21331 features, 0 variable features)

sample_1030v2<-Read10X_h5("Partek_scRNA_Feb_22-cellranger_new_exported (1)/filtered_feature_bc_matrix1030v2.h5")
sample_1030v2<-CreateSeuratObject(sample_1005, project = "sample_1030v2")
sample_1030v2
#An object of class Seurat 
#21331 features across 491 samples within 1 assay 
#Active assay: RNA (21331 features, 0 variable features)

sample_1031v1<-Read10X_h5("Partek_scRNA_Feb_22-cellranger_new_exported (1)/filtered_feature_bc_matrix1031v1.h5")
sample_1031v1<-CreateSeuratObject(sample_1005, project = "sample_1031v1")
sample_1031v1
#An object of class Seurat 
#21331 features across 491 samples within 1 assay 
#Active assay: RNA (21331 features, 0 variable features)

sample_1033v1<-Read10X_h5("Partek_scRNA_Feb_22-cellranger_new_exported (1)/filtered_feature_bc_matrix1033v1.h5")
sample_1033v1<-CreateSeuratObject(sample_1005, project = "sample_1033v1")
sample_1033v1
#An object of class Seurat 
#21331 features across 491 samples within 1 assay 
#Active assay: RNA (21331 features, 0 variable features)
