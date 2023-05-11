setwd("/Users/liuzhe/Desktop/UCSF/milk/analysis1/six_comparision_addIandJ/lastedversion/")
rm(list=ls())

#####normal breast sample sequencing data
rawdata <- read.table("/Users/liuzhe/Desktop/UCSF/milk/analysis1/datafile/Partek_cells_vs._fat_vs._breast_Quantify_to_annotation_model_(Partek_E_M)_Gene_counts.txt",sep="\t",header=T,row.names = 1)
metadata<-read.table("/Users/liuzhe/Desktop/UCSF/milk/analysis1/datafile/metadata.txt",header = T, sep = "\t")
dim(rawdata)
#[1]     9 30728
breast<-rawdata[7:9,4:30728]
breast_t<-t(breast)
breast_counts<-as.data.frame(breast_t)
breast_counts$ensemblID<-row.names(breast_counts)
#gene ID 转换
#reference： http://rvdsd.top/2018/06/24/BIoTools/%E5%90%84%E7%A7%8D%E5%9F%BA%E5%9B%A0%E5%90%8D%E8%BD%AC%E6%8D%A2/
library(biomaRt)
trans_gene_id <- breast_counts$ensemblID
ensembl <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
save(ensembl, file = "ensembl_human.RData")
result_gene_entrez <- getBM(values = trans_gene_id,
                            filters = 'ensembl_gene_id', 
                            mart = ensembl,
                            attributes = c('hgnc_symbol', 'ensembl_gene_id'))
dim(result_gene_entrez)
#[1] 30678     2
merged_breast_counts<-merge(breast_counts, result_gene_entrez, by.x = "ensemblID", by.y = "ensembl_gene_id")
counts_anno<-aggregate(merged_breast_counts,by=list(merged_breast_counts$hgnc_symbol),FUN=mean)
counts_anno<-counts_anno[,c(1,3,4,5)]
counts_anno<-counts_anno[-1,]
dim(counts_anno)
#[1] 23524     4
write.csv(counts_anno, "breast_readcounts_anno.csv", quote = F, row.names = F)


#####cells and fat layer samples sequencing data
sampleA<-read.table("/Users/liuzhe/Desktop/UCSF/milk/analysis1/addSample_IJ/A_featurecounts.txt", header = T)
sampleB<-read.table("/Users/liuzhe/Desktop/UCSF/milk/analysis1/addSample_IJ/B_featurecounts.txt", header = T)
sampleC<-read.table("/Users/liuzhe/Desktop/UCSF/milk/analysis1/addSample_IJ/C_featurecounts.txt", header = T)
sampleD<-read.table("/Users/liuzhe/Desktop/UCSF/milk/analysis1/addSample_IJ/D_featurecounts.txt", header = T)
sampleE<-read.table("/Users/liuzhe/Desktop/UCSF/milk/analysis1/addSample_IJ/E_featurecounts.txt", header = T)
sampleF<-read.table("/Users/liuzhe/Desktop/UCSF/milk/analysis1/addSample_IJ/F_featurecounts.txt", header = T)
sampleI<-read.table("/Users/liuzhe/Desktop/UCSF/milk/analysis1/addSample_IJ/I_featurecounts.txt", header = T)
sampleJ<-read.table("/Users/liuzhe/Desktop/UCSF/milk/analysis1/addSample_IJ/J_featurecounts.txt", header = T)
sampleinfo<-read.csv("/Users/liuzhe/Desktop/UCSF/milk/analysis1/addSample_IJ/Metadata.csv")
dim(sampleA)
#[1] 60642     7
dim(sampleB)
#[1] 60642     7
dim(sampleC)
#[1] 60642     7
dim(sampleD)
#[1] 60642     7
dim(sampleE)
#[1] 60642     7
dim(sampleF)
#[1] 60642     7
dim(sampleI)
#[1] 60642     7
dim(sampleJ)
#[1] 60642     7
readcounts<-cbind(sampleA, sampleB$result.BAligned.sortedByCoord.out.bam.sorted, 
                  sampleC$result.CAligned.sortedByCoord.out.bam.sorted, sampleD$result.DAligned.sortedByCoord.out.bam.sorted,
                  sampleD$result.DAligned.sortedByCoord.out.bam.sorted, sampleE$result.EAligned.sortedByCoord.out.bam.sorted,
                  sampleF$result.FAligned.sortedByCoord.out.bam.sorted, sampleI$result.IAligned.sortedByCoord.out.bam.sorted,
                  sampleJ$result.JAligned.sortedByCoord.out.bam.sorted)
readcounts<-subset(readcounts, select = c("Geneid", "Length", "result.AAligned.sortedByCoord.out.bam.sorted",
                                          "sampleB$result.BAligned.sortedByCoord.out.bam.sorted",
                                          "sampleC$result.CAligned.sortedByCoord.out.bam.sorted",
                                          "sampleD$result.DAligned.sortedByCoord.out.bam.sorted",
                                          "sampleE$result.EAligned.sortedByCoord.out.bam.sorted",
                                          "sampleF$result.FAligned.sortedByCoord.out.bam.sorted",
                                          "sampleI$result.IAligned.sortedByCoord.out.bam.sorted",
                                          "sampleJ$result.JAligned.sortedByCoord.out.bam.sorted"))
colnames(readcounts)<-c("Geneid", "Length", "A", "B", "C", "D", "E", "F", "I", "J")
allsamples_counts<-merge(readcounts, counts_anno, by.x = "Geneid", by.y = "Group.1")
write.csv(allsamples_counts, "allsamples_readcounts_anno.csv", quote = F, row.names = F)
row.names(allsamples_counts)<-allsamples_counts$Geneid
allsamples_counts$Geneid<-NULL
head(allsamples_counts)

library("DESeq2")
library("limma")
library("edgeR")
library("dplyr") 
library("biomaRt")
library("tidyverse")
library("stringr")
allsamples_counts<-read.csv("/Users/liuzhe/Desktop/UCSF/milk/analysis1/six_comparision_addIandJ/allsamples_readcounts_anno.csv", header = T, row.names = 1)
kb<-allsamples_counts$Length/1000
rpk <- allsamples_counts[2:12] / kb
rpk
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
write.csv(tpm,file="tpm.csv",quote=F)

####################compare1: identification of DEGs between fat layers vs cells#########################################
setwd("/Users/liuzhe/Desktop/UCSF/milk/analysis1/six_comparision_addIandJ/lastedversion/Fat_Cell")
rawdata<-subset(allsamples_counts, select = c("B", "D", "F", "J", "A", "C", "E", "I"))
dim(rawdata)
#[1] 21993     8
gset <- rawdata[rowMeans(rawdata)>0,] # remove low expressed genes
dim(gset)
#[1] 20664     8
coldata <- data.frame(
  condition = factor(c("milk_fat_layer", "milk_fat_layer", "milk_fat_layer", "milk_fat_layer",
                       "milk_cells", "milk_cells", "milk_cells", "milk_cells")))
str(coldata)
condition = factor(c("milk_fat_layer", "milk_fat_layer", "milk_fat_layer", "milk_fat_layer",
                     "milk_cells", "milk_cells", "milk_cells", "milk_cells"))
subject <- factor(c("ID1012","ID1013","ID1014","ID1012","ID1012","ID1013","ID1014","ID1012"))
rownames(coldata)<-colnames(gset)
coldata$subject<-subject
# add paired information into design
dds <- DESeqDataSetFromMatrix(countData = gset,
                              colData = coldata,
                              design = ~ subject + condition) 
dds$condition<- relevel(dds$condition, ref = "milk_cells") 
dds <- DESeq(dds)
nrDEG_DESeq2 <- as.data.frame(results(dds))
rld <- rlog(dds)
# we exported the normalized expression matrix for further analysis
normal_gset <- assay(rld) 
nrDEG_DESeq2 = nrDEG_DESeq2[order(nrDEG_DESeq2$log2FoldChange),] 
## define DEGs
nrDEG <- nrDEG_DESeq2
nrDEG$Group = "notsignificant"
logFC_cutoff <- 2
nrDEG$Group[which( (nrDEG$padj < 0.01) & (nrDEG$log2FoldChange > logFC_cutoff) )] = "upregulated"
nrDEG$Group[which( (nrDEG$padj < 0.01) & (nrDEG$log2FoldChange < -logFC_cutoff) )] = "downregulated"
table(nrDEG$Group)
#downregulated notsignificant    upregulated 
#         780          19876            8
write.csv(nrDEG,"DEGs_for_Fat-cell.csv",quote = F)
library("ggplot2")
deg<-nrDEG
head(deg)
deg$color<-ifelse(deg$padj<0.01 & abs(deg$log2FoldChange)>logFC_cutoff, ifelse(deg$log2FoldChange< -logFC_cutoff, "blue", "red"),"gray")
color<-c(red = "red", gray = "gray", blue = "blue")
deg$symbol<-row.names(deg)
p <- ggplot(data = deg, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=Group)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  theme_bw() +
  xlim(-25,25)

riskgenes = c("SCARNA21", "FCGR3B", "ANXA1", "EMP1", "KRT23", "CD53", "CORO1A", "CLDN1", "KLF6", "PMAIP1", 
              "PGLYRP1", "SELL", "PROK2", "AIF1", "KRT23", "CXCR2", "CD69", "GPR65", "FCGR3B", "HLA-DOA",
              "UBE4A", "FGFBP1", "HERC2", "TACC2", "SSNA1", "LNPEP", "PRKDC", "VPS13D", "MRPL41", "EMILIN3",
              "PAX7", "DNAI1", "LINC01785", "PPP1R12A-AS1", "GGTLC1", "FGF14-AS2", "AGAP2-AS1","ZNF572",
              "GLP1R", "CEROX1", "PAX7", "DNAI1", "LINC01785","PPP1R12A-AS1", "GGTLC1", "FGF14-AS2",
              "AGAP2-AS1", "ZNF572", "GLP1R", "CEROX1", "HOXA10", "HOXB8", "UBE4A", "RPL17P25", "GGACT",
              "FBXL22", "GEMIN8P4", "LRRC37A7P", "LRRC26", "GAMT", "HOXB7", "GNAQP1", "SMIM22", "CALML5",
              "TACC2", "TMEM160", "EMILIN3", "PRKDC", "CRLS1", "MYORG","GAMT","LRRC26", "LRRC37A7P","UBE4A",
              "HOXB8","HOXA10","CEROX1", "LINC01785")
riskgenes <- unique(riskgenes)
for_label<-deg[riskgenes,]
#for_label <- deg %>% 
#  filter(abs(log2FoldChange) >1& padj < 0.05)
pdf("VP_DEGs_Fat-cell_add.pdf")
p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )+ theme(legend.position="none")
dev.off()
#riskgenes = c("gene1", "gene2", "gene3", "gene4")
#for_label<-deg[riskgenes,]
normal_gset<-as.data.frame(normal_gset)
deg_down<-subset(deg, Group == "downregulated")
deg_up<-subset(deg, Group == "upregulated")
deg_select<-rbind(deg_down,deg_up)
normalized_exp<-normal_gset[rownames(deg_select),]
exp_deg<-cbind(deg_select,normalized_exp)
exp_deg<-subset(exp_deg, select = c("B", "D", "F", "J", "A", "C", "E", "I", "log2FoldChange", "pvalue", "padj", "Group"))
colnames(exp_deg)<-c("FatLayer_B", "FatLayer_D", "FatLayer_F", "FatLayer_J", 
                     "Cells_A", "Cells_C", "Cells_E", "Cells_I",
                     "log2FoldChange", "raw_Pvalue","adjustedPvalue","Regulation")
#exp_deg<-cbind(round(exp_deg[1:10],3),exp_deg[,11])
colnames(exp_deg)[12]<-"regulation"
write.csv(exp_deg,"Table1_Fat-cell.csv",quote = F)
library("pheatmap")
rt<-exp_deg[,1:8]
metadata_fat_cell<-sampleinfo[c(1:6,9:10),]
annotation_col = data.frame(Type = factor(metadata_fat_cell$Type))
rownames(annotation_col) = c("Cells_A", "FatLayer_B", "Cells_C", "FatLayer_D", "Cells_E", "FatLayer_F", "Cells_I", "FatLayer_J")
annotation_row = matrix(exp_deg$regulation,nrow=length(exp_deg$regulation), ncol=1,byrow=TRUE)
rownames(annotation_row) = rownames(exp_deg)
colnames(annotation_row)<-"gene_type"
annotation_row<-as.data.frame(annotation_row)
ann_colors = list(Type = c(Cells = "blue", Fat = "red"),
                  gene_type = c(downregulated = "blue", upregulated  = "red"))
pdf(file="DESeq2_median_of_ratios/heatmap_Fat-cell.pdf",width = 6,height = 8)
pheatmap(rt,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 780, gaps_col = 4, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
rt_top50<-rbind(rt[1:25,],rt[781:788,])
pdf(file="DESeq2_median_of_ratios/heatmap_Fat-cell_selected_top50.pdf",width = 6,height = 8)
pheatmap(rt_top50,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 25, gaps_col = 4, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
tpm_selected<-subset(tpm, select = c("A", "B", "C", "D", "E", "F", "J", "I"))
colnames(tpm_selected)<-c("Cells_A","FatLayer_B","Cells_C","FatLayer_D","Cells_E","FatLayer_F","FatLayer_J","Cells_I")
tpm_selected<-subset(tpm_selected, select=c("FatLayer_B","FatLayer_D","FatLayer_F","FatLayer_J","Cells_A","Cells_C","Cells_E","Cells_I"))
selected_genes<-rownames(rt)
tpm_data<-tpm_selected[selected_genes,]
pdf(file="TPM/heatmap_Fat-cell.pdf",width = 6,height = 8)
pheatmap(tpm_data,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 780, gaps_col = 4, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
tpm_data_top50<-rbind(tpm_data[1:25,],tpm_data[781:788,])
pdf(file="TPM/heatmap_Fat-cell_selected_top50.pdf",width = 6,height = 8)
pheatmap(tpm_data_top50,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 25, gaps_col = 4, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()

library(ggpubr)
library("DOSE")
library("org.Hs.eg.db")
library("topGO")
library("clusterProfiler")
library("pathview")
###up DEG functional enrichment analysis
deg_up<-subset(exp_deg, regulation == "upregulated")
up_genes <- rownames(deg_up)
test = bitr(up_genes, #dataset
            fromType="SYMBOL",
            toType="ENTREZID", 
            OrgDb="org.Hs.eg.db") 
head(test,2)
go_BP <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1,
                  keyType = 'ENTREZID')
y <- setReadable(go_BP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
### The geneID column is translated to symbol
write.csv(y,"GOBP_enrich_up_Fat-cell_readable.csv",row.names =FALSE)
#visualization
GOBP_up<-summary(go_BP)
GOBP_up<-GOBP_up[1:20,]
GOBP_up$yvalue<-GOBP_up$Count/5
GOBP_up$xvalue<-GOBP_up$Description
library(ggpubr)
pdf("GOBP_up_Fat-cell.barplot.pdf",height=6,width = 12)
ggbarplot(GOBP_up, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("red","pink"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

#GOBP_up$minuslog10Pval=-log10(GOBP_up$pvalue)
#pdf("GOBP_up_Fat-cell.barplot_log_transformed_pvalue.pdf",height=6,width = 12)
#ggbarplot(GOBP_up, x = "xvalue", y = "yvalue", orientation = "horiz",
#          xlab = "", ylab="Gene Ratio",
#          fill = "minuslog10Pval",
#          color = "white",
#          width = 0.6,
#          position = position_dodge())+gradient_fill(c("pink","red"))+
#  theme(legend.key.size = unit(0.3, "cm"))+
#  theme(legend.position = "right")
#dev.off()

###down DEG functional enrichment analysis
deg_down<-subset(exp_deg, regulation == "downregulated")
down_genes <- rownames(deg_down)
test = bitr(down_genes, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
head(test,2)
go_BP <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1,
                  keyType = 'ENTREZID')
#visualization
y <- setReadable(go_BP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
### The geneID column is translated to symbol
write.csv(y,"GOBP_enrich_down_Fat-cell_readable.csv",row.names =FALSE)
GOBP_down<-summary(go_BP)
GOBP_down<-GOBP_down[1:20,]
GOBP_down$yvalue<-GOBP_down$Count/687
GOBP_down$xvalue<-GOBP_down$Description
pdf("GOBP_down.barplot_Fat-cell.pdf",height=6,width = 12)
ggbarplot(GOBP_down, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("blue","purple"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

GOBP_down$minuslog10Pval=-log10(GOBP_down$pvalue)
pdf("GOBP_down_Fat-cell.barplot_log_transformed_pvalue.pdf",height=6,width = 12)
ggbarplot(GOBP_down, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "minuslog10Pval",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("purple","blue"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()




########compare2: identification of DEGs between fat layers vs normal breast tissue#########################################
setwd("/Users/liuzhe/Desktop/UCSF/milk/analysis1/six_comparision_addIandJ/lastedversion/Fat_Breast")
rawdata<-subset(allsamples_counts, select = c("B", "D", "F", "J", "SRR2148259.sra", "SRR2148260.sra", "SRR2148261"))
dim(rawdata)
#[1] 21993     7
gset <- rawdata[rowMeans(rawdata)>0,] # 剔除表达量低的基因
dim(gset)
#[1] 21936     7
coldata <- data.frame(
  condition = factor(c("milk_fat_layer", "milk_fat_layer", "milk_fat_layer", "milk_fat_layer",
                       "normal_breast", "normal_breast", "normal_breast")))
str(coldata)
condition = factor(c("milk_fat_layer", "milk_fat_layer", "milk_fat_layer", "milk_fat_layer",
                     "normal_breast", "normal_breast", "normal_breast"))
subject <- factor(c("ID1012","ID1013","ID1014","ID1012","B1","B2","B3"))
rownames(coldata)<-colnames(gset)
coldata$subject<-subject
# 注意在design中加上配对信息
data = apply(gset, 2, as.integer) ## DESeq2分析需要是整数
row.names(data)<-row.names(gset)
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design =  ~ condition) 
dds$condition<- relevel(dds$condition, ref = "normal_breast") 
dds <- DESeq(dds)
nrDEG_DESeq2 <- as.data.frame(results(dds))
rld <- rlog(dds)
# 这里我还提取了标准化后的表达矩阵，可以用于后续的热图绘制等等
normal_gset <- assay(rld) 
nrDEG_DESeq2 = nrDEG_DESeq2[order(nrDEG_DESeq2$log2FoldChange),] 
## 4.4定义差异基因
nrDEG <- nrDEG_DESeq2
nrDEG$Group = "notsignificant"
logFC_cutoff <- 2
nrDEG$Group[which( (nrDEG$padj < 0.01) & (nrDEG$log2FoldChange > logFC_cutoff) )] = "upregulated"
nrDEG$Group[which( (nrDEG$padj < 0.01) & (nrDEG$log2FoldChange < -logFC_cutoff) )] = "downregulated"
table(nrDEG$Group)
#downregulated notsignificant    upregulated 
#         3884          16729           1323
write.csv(nrDEG,"DEGs_for_Fat-breast.csv",quote = F)
library("ggplot2")
deg<-nrDEG
head(deg)
deg$color<-ifelse(deg$padj<0.01 & abs(deg$log2FoldChange)>logFC_cutoff, ifelse(deg$log2FoldChange< -logFC_cutoff, "blue", "red"),"gray")
color<-c(red = "red", gray = "gray", blue = "blue")
deg$symbol<-row.names(deg)
p <- ggplot(data = deg, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=Group)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  theme_bw() +
  xlim(-25,25)
for_label <- deg %>% 
  filter(abs(log2FoldChange) >2& padj < 0.01)
pdf("VP_DEGs_Fat-breast.pdf")
p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )+ theme(legend.position="none")
dev.off()
#riskgenes = c("gene1", "gene2", "gene3", "gene4")
#for_label<-deg[riskgenes,]
normal_gset<-as.data.frame(normal_gset)
deg_down<-subset(deg, Group == "downregulated")
deg_up<-subset(deg, Group == "upregulated")
deg_select<-rbind(deg_down,deg_up)
normalized_exp<-normal_gset[rownames(deg_select),]
exp_deg<-cbind(deg_select,normalized_exp)
exp_deg<-subset(exp_deg, select = c("B", "D", "F", "J", "SRR2148259.sra", "SRR2148260.sra", "SRR2148261",
                                    "log2FoldChange", "pvalue", "padj", "Group"))
colnames(exp_deg)<-c("FatLayer_B", "FatLayer_D", "FatLayer_F", "FatLayer_J", 
                     "NormalBreast_SRR2148259", "NormalBreast_SRR2148260", "NormalBreast_SRR2148261",
                     "log2FoldChange", "raw_Pvalue","adjustedPvalue","Regulation")
#exp_deg<-cbind(round(exp_deg[1:10],3),exp_deg[,11])
colnames(exp_deg)[11]<-"regulation"
write.csv(exp_deg,"Table2_Fat-breast.csv",quote = F)
library("pheatmap")
rt<-exp_deg[,1:7]
annotation_col = data.frame(Type = c("Fat", "Fat", "Fat", "Fat", "NormalBreast", "NormalBreast", "NormalBreast"))
annotation_col$Type<-as.factor(annotation_col$Type)
rownames(annotation_col) = c("FatLayer_B", "FatLayer_D", "FatLayer_F", "FatLayer_J", 
                             "NormalBreast_SRR2148259", "NormalBreast_SRR2148260", "NormalBreast_SRR2148261")
annotation_row = matrix(exp_deg$regulation,nrow=length(exp_deg$regulation), ncol=1,byrow=TRUE)
rownames(annotation_row) = rownames(exp_deg)
colnames(annotation_row)<-"gene_type"
annotation_row<-as.data.frame(annotation_row)
ann_colors = list(Type = c(NormalBreast = "blue", Fat = "red"),
                  gene_type = c(downregulated = "blue", upregulated  = "red"))

pdf(file="DESeq2_median_of_ratios/heatmap_Fat-breast.pdf",width = 6,height = 8)
pheatmap(rt,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 3884, gaps_col = 4, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
rt_top50<-rbind(rt[1:25,],rt[3885:3909,])
pdf(file="DESeq2_median_of_ratios/heatmap_Fat-breast_selected_top50.pdf",width = 6,height = 8)
pheatmap(rt_top50,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 25, gaps_col = 4, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()

tpm_selected<-subset(tpm, select = c("B", "D", "F", "J", "SRR2148259.sra", "SRR2148260.sra", "SRR2148261"))
colnames(tpm_selected)<-colnames(rt)
tpm_selected<-subset(tpm_selected, select=c("FatLayer_B","FatLayer_D","FatLayer_F","FatLayer_J","NormalBreast_SRR2148259","NormalBreast_SRR2148260","NormalBreast_SRR2148261"))
selected_genes<-rownames(rt)
tpm_data<-tpm_selected[selected_genes,]
pdf(file="TPM/heatmap_Fat-breast.pdf",width = 6,height = 8)
pheatmap(tpm_data,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 3884, gaps_col = 4, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
tpm_data_top50<-rbind(tpm_data[1:25,],tpm_data[3885:3909,])
pdf(file="TPM/heatmap_Fat-breast_selected_top50.pdf",width = 6,height = 8)
pheatmap(tpm_data_top50,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 25, gaps_col = 4, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()

library(ggpubr)
library("DOSE")
library("org.Hs.eg.db")
library("topGO")
library("clusterProfiler")
library("pathview")
###up DEG functional enrichment analysis
deg_up<-subset(exp_deg, regulation == "upregulated")
up_genes <- rownames(deg_up)
test = bitr(up_genes, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
head(test,2)
go_BP <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1,
                  keyType = 'ENTREZID')
y <- setReadable(go_BP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
### The geneID column is translated to symbol
write.csv(y,"GOBP_enrich_up_Fat-breast_readable.csv",row.names =FALSE)
GOBP_up<-summary(go_BP)
GOBP_up<-GOBP_up[1:20,]
GOBP_up$yvalue<-GOBP_up$Count/1030
GOBP_up$xvalue<-GOBP_up$Description
library(ggpubr)
pdf("GOBP_up_Fat-breast.barplot.pdf",height=6,width = 12)
ggbarplot(GOBP_up, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("red","pink"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

GOBP_up$minuslog10Pval=-log10(GOBP_up$pvalue)
pdf("GOBP_up_Fat-breast.barplot_log_transformed_pvalue.pdf",height=6,width = 12)
ggbarplot(GOBP_up, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "minuslog10Pval",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("pink","red"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

###down DEG functional enrichment analysis
deg_down<-subset(exp_deg, regulation == "downregulated")
down_genes <- rownames(deg_down)
test = bitr(down_genes, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
head(test,2)
go_BP <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1,
                  keyType = 'ENTREZID')
y <- setReadable(go_BP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
### The geneID column is translated to symbol
write.csv(y,"GOBP_enrich_down_Fat-breast_readable.csv",row.names =FALSE)
GOBP_down<-summary(go_BP)
GOBP_down<-GOBP_down[1:20,]
GOBP_down$yvalue<-GOBP_down$Count/2773
GOBP_down$xvalue<-GOBP_down$Description
pdf("GOBP_down.barplot_Fat-breast.pdf",height=6,width = 12)
ggbarplot(GOBP_down, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("blue","purple"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

GOBP_down$minuslog10Pval=-log10(GOBP_down$pvalue)
pdf("GOBP_down_Fat-breast.barplot_log_transformed_pvalue.pdf",height=6,width = 12)
ggbarplot(GOBP_down, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "minuslog10Pval",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("purple","blue"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()




##############compare3:identification of DEGs between cells vs normal breast tissue#########################################
setwd("/Users/liuzhe/Desktop/UCSF/milk/analysis1/six_comparision_addIandJ/lastedversion/Cell_Breast")
rawdata<-subset(allsamples_counts, select = c("A", "C", "E", "I", "SRR2148259.sra", "SRR2148260.sra", "SRR2148261"))
dim(rawdata)
#[1] 21993     7
gset <- rawdata[rowMeans(rawdata)>0,] # 剔除表达量低的基因
dim(gset)
#[1] 21960     7
coldata <- data.frame(
  condition = factor(c("milk_cell", "milk_cell", "milk_cell", "milk_cell",
                       "normal_breast", "normal_breast", "normal_breast")))
str(coldata)
condition = factor(c("milk_cell", "milk_cell", "milk_cell", "milk_cell",
                     "normal_breast", "normal_breast", "normal_breast"))
subject <- factor(c("ID1012","ID1013","ID1014","ID1012","B1","B2","B3"))
rownames(coldata)<-colnames(gset)
coldata$subject<-subject
# 注意在design中加上配对信息
data = apply(gset, 2, as.integer) ## DESeq2分析需要是整数
row.names(data)<-row.names(gset)
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design =  ~ condition) 
dds$condition<- relevel(dds$condition, ref = "normal_breast") 
dds <- DESeq(dds)
nrDEG_DESeq2 <- as.data.frame(results(dds))
rld <- rlog(dds)
# 这里我还提取了标准化后的表达矩阵，可以用于后续的热图绘制等等
normal_gset <- assay(rld) 
nrDEG_DESeq2 = nrDEG_DESeq2[order(nrDEG_DESeq2$log2FoldChange),] 
## 4.4定义差异基因
nrDEG <- nrDEG_DESeq2
nrDEG$Group = "notsignificant"
logFC_cutoff <- 2
nrDEG$Group[which( (nrDEG$padj < 0.01) & (nrDEG$log2FoldChange > logFC_cutoff) )] = "upregulated"
nrDEG$Group[which( (nrDEG$padj < 0.01) & (nrDEG$log2FoldChange < -logFC_cutoff) )] = "downregulated"
table(nrDEG$Group)
#downregulated notsignificant    upregulated 
#         3769          16633           1558
write.csv(nrDEG,"DEGs_for_Cell-breast.csv",quote = F)
library("ggplot2")
deg<-nrDEG
head(deg)
deg$color<-ifelse(deg$padj<0.01 & abs(deg$log2FoldChange)>logFC_cutoff, ifelse(deg$log2FoldChange< -logFC_cutoff, "blue", "red"),"gray")
color<-c(red = "red", gray = "gray", blue = "blue")
deg$symbol<-row.names(deg)
p <- ggplot(data = deg, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=Group)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  theme_bw() +
  xlim(-25,25)
for_label <- deg %>% 
  filter(abs(log2FoldChange) >2& padj < 0.01)
pdf("VP_DEGs_Cell-breast.pdf")
p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )+ theme(legend.position="none")
dev.off()
#riskgenes = c("gene1", "gene2", "gene3", "gene4")
#for_label<-deg[riskgenes,]
normal_gset<-as.data.frame(normal_gset)
deg_down<-subset(deg, Group == "downregulated")
deg_up<-subset(deg, Group == "upregulated")
deg_select<-rbind(deg_down,deg_up)
normalized_exp<-normal_gset[rownames(deg_select),]
exp_deg<-cbind(deg_select,normalized_exp)
exp_deg<-subset(exp_deg, select = c("A", "C", "E", "I", "SRR2148259.sra", "SRR2148260.sra", "SRR2148261",
                                    "log2FoldChange", "pvalue", "padj", "Group"))
colnames(exp_deg)<-c("Cell_A", "Cell_C", "Cell_E", "Cell_I", 
                     "NormalBreast_SRR2148259", "NormalBreast_SRR2148260", "NormalBreast_SRR2148261",
                     "log2FoldChange", "raw_Pvalue","adjustedPvalue","Regulation")
#exp_deg<-cbind(round(exp_deg[1:10],3),exp_deg[,11])
colnames(exp_deg)[11]<-"regulation"
write.csv(exp_deg,"Table3_Cell-breast.csv",quote = F)
library("pheatmap")
rt<-exp_deg[,1:7]

annotation_col = data.frame(Type = c("Cell", "Cell", "Cell", "Cell", "NormalBreast", "NormalBreast", "NormalBreast"))
annotation_col$Type<-as.factor(annotation_col$Type)
rownames(annotation_col) = c("Cell_A", "Cell_C", "Cell_E", "Cell_I", 
                             "NormalBreast_SRR2148259", "NormalBreast_SRR2148260", "NormalBreast_SRR2148261")
annotation_row = matrix(exp_deg$regulation,nrow=length(exp_deg$regulation), ncol=1,byrow=TRUE)
rownames(annotation_row) = rownames(exp_deg)
colnames(annotation_row)<-"gene_type"
annotation_row<-as.data.frame(annotation_row)
ann_colors = list(Type = c(NormalBreast = "blue", Cell = "red"),
                  gene_type = c(downregulated = "blue", upregulated  = "red"))

pdf(file="DESeq2_median_of_ratios/heatmap_Cell-breast.pdf",width = 6,height = 8)
pheatmap(rt,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 3769, gaps_col = 4, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
rt_top50<-rbind(rt[1:25,],rt[3770:3794,])
pdf(file="DESeq2_median_of_ratios/heatmap_Cell-breast_selected_top50.pdf",width = 6,height = 8)
pheatmap(rt_top50,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 25, gaps_col = 4, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()

tpm_selected<-subset(tpm, select = c("A", "C", "E", "I", "SRR2148259.sra", "SRR2148260.sra", "SRR2148261"))
colnames(tpm_selected)<-colnames(rt)
selected_genes<-rownames(rt)
tpm_data<-tpm_selected[selected_genes,]
pdf(file="TPM/heatmap_Cell-breast.pdf",width = 6,height = 8)
pheatmap(tpm_data,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 3769, gaps_col = 4, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
tpm_data_top50<-rbind(tpm_data[1:25,],tpm_data[3770:3794,])
pdf(file="TPM/heatmap_Cell-breast_selected_top50.pdf",width = 6,height = 8)
pheatmap(tpm_data_top50,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 25, gaps_col = 4, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()

library(ggpubr)
library("DOSE")
library("org.Hs.eg.db")
library("topGO")
library("clusterProfiler")
library("pathview")
###up DEG functional enrichment analysis
deg_up<-subset(exp_deg, regulation == "upregulated")
up_genes <- rownames(deg_up)
test = bitr(up_genes, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
head(test,2)
go_BP <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1,
                  keyType = 'ENTREZID')
y <- setReadable(go_BP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
### The geneID column is translated to symbol
write.csv(y,"GOBP_enrich_up_Cell-breast_readable.csv",row.names =FALSE)
#结果可视化
GOBP_up<-summary(go_BP)
GOBP_up<-GOBP_up[1:20,]
GOBP_up$yvalue<-GOBP_up$Count/1177
GOBP_up$xvalue<-GOBP_up$Description
library(ggpubr)
pdf("GOBP_up_Cell-breast.barplot.pdf",height=6,width = 12)
ggbarplot(GOBP_up, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("red","pink"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

GOBP_up$minuslog10Pval=-log10(GOBP_up$pvalue)
pdf("GOBP_up_Cell-breast.barplot_log_transformed_pvalue.pdf",height=6,width = 12)
ggbarplot(GOBP_up, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "minuslog10Pval",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("pink","red"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

###down DEG functional enrichment analysis
deg_down<-subset(exp_deg, regulation == "downregulated")
down_genes <- rownames(deg_down)
test = bitr(down_genes, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
head(test,2)
go_BP <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1,
                  keyType = 'ENTREZID')
y <- setReadable(go_BP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
### The geneID column is translated to symbol
write.csv(y,"GOBP_enrich_down_Cell-breast_readable.csv",row.names =FALSE)
GOBP_down<-summary(go_BP)
GOBP_down<-GOBP_down[1:20,]
GOBP_down$yvalue<-GOBP_down$Count/2491
GOBP_down$xvalue<-GOBP_down$Description
pdf("GOBP_down.barplot_Cell-breast.pdf",height=6,width = 12)
ggbarplot(GOBP_down, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("blue","purple"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

GOBP_down$minuslog10Pval=-log10(GOBP_down$pvalue)
pdf("GOBP_down_Cell-breast.barplot_log_transformed_pvalue.pdf",height=6,width = 12)
ggbarplot(GOBP_down, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "minuslog10Pval",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("purple","blue"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()




##############compare4:identification of DEGs between fatcell vs normal breast tissue#########################################
setwd("/Users/liuzhe/Desktop/UCSF/milk/analysis1/six_comparision_addIandJ/lastedversion/FatCell_Breast")
rawdata<-subset(allsamples_counts, select = c("A","C","E","I","B","D","F","J","SRR2148259.sra","SRR2148260.sra","SRR2148261"))
dim(rawdata)
#[1] 21993     11
gset <- rawdata[rowMeans(rawdata)>0,] # 剔除表达量低的基因
dim(gset)
#[1] 21967     11
coldata <- data.frame(
  condition = factor(c("milk", "milk", "milk", "milk",
                       "milk", "milk", "milk", "milk",
                       "normal_breast", "normal_breast", "normal_breast")))
str(coldata)
condition = factor(c("milk", "milk", "milk", "milk",
                     "milk", "milk", "milk", "milk",
                     "normal_breast", "normal_breast", "normal_breast"))
subject <- factor(c("ID1012","ID1013","ID1014","ID1012",
                    "ID1012","ID1013","ID1014","ID1012",
                    "B1","B2","B3"))
rownames(coldata)<-colnames(gset)
coldata$subject<-subject
# 注意在design中加上配对信息
data = apply(gset, 2, as.integer) ## DESeq2分析需要是整数
row.names(data)<-row.names(gset)
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design =  ~ condition) 
dds$condition<- relevel(dds$condition, ref = "normal_breast") 
dds <- DESeq(dds)
nrDEG_DESeq2 <- as.data.frame(results(dds))
rld <- rlog(dds)
# 这里我还提取了标准化后的表达矩阵，可以用于后续的热图绘制等等
normal_gset <- assay(rld) 
nrDEG_DESeq2 = nrDEG_DESeq2[order(nrDEG_DESeq2$log2FoldChange),] 
## 4.4定义差异基因
nrDEG <- nrDEG_DESeq2
nrDEG$Group = "notsignificant"
logFC_cutoff <- 2
nrDEG$Group[which( (nrDEG$padj < 0.01) & (nrDEG$log2FoldChange > logFC_cutoff) )] = "upregulated"
nrDEG$Group[which( (nrDEG$padj < 0.01) & (nrDEG$log2FoldChange < -logFC_cutoff) )] = "downregulated"
table(nrDEG$Group)
#downregulated notsignificant    upregulated 
#         2819          17921           1227 
write.csv(nrDEG,"DEGs_for_FatCell-breast.csv",quote = F)
library("ggplot2")
deg<-nrDEG
head(deg)
deg$color<-ifelse(deg$padj<0.01 & abs(deg$log2FoldChange)>logFC_cutoff, ifelse(deg$log2FoldChange< -logFC_cutoff, "blue", "red"),"gray")
color<-c(red = "red", gray = "gray", blue = "blue")
deg$symbol<-row.names(deg)
p <- ggplot(data = deg, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=Group)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  theme_bw() +
  xlim(-25,25)
for_label <- deg %>% 
  filter(abs(log2FoldChange) >2& padj < 0.01)
pdf("VP_DEGs_FatCell-breast.pdf")
p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )+ theme(legend.position="none")
dev.off()
#riskgenes = c("gene1", "gene2", "gene3", "gene4")
#for_label<-deg[riskgenes,]
normal_gset<-as.data.frame(normal_gset)
deg_down<-subset(deg, Group == "downregulated")
deg_up<-subset(deg, Group == "upregulated")
deg_select<-rbind(deg_down,deg_up)
normalized_exp<-normal_gset[rownames(deg_select),]
exp_deg<-cbind(deg_select,normalized_exp)
exp_deg<-subset(exp_deg, select = c("A", "C", "E", "I", "B", "D", "F", "J", "SRR2148259.sra", "SRR2148260.sra", "SRR2148261",
                                    "log2FoldChange", "pvalue", "padj", "Group"))
colnames(exp_deg)<-c("Cell_A", "Cell_C", "Cell_E", "Cell_I", 
                     "FatLayer_B",	"FatLayer_D",	"FatLayer_F",	"FatLayer_J",
                     "NormalBreast_SRR2148259", "NormalBreast_SRR2148260", "NormalBreast_SRR2148261",
                     "log2FoldChange", "raw_Pvalue","adjustedPvalue","Regulation")
#exp_deg<-cbind(round(exp_deg[1:10],3),exp_deg[,11])
colnames(exp_deg)[15]<-"regulation"
write.csv(exp_deg,"Table4_FatCell-breast.csv",quote = F)
library("pheatmap")
rt<-exp_deg[,1:11]
annotation_col = data.frame(Type = c("Cell", "Cell", "Cell", "Cell",
                                     "Fat", "Fat", "Fat", "Fat", 
                                     "NormalBreast", "NormalBreast", "NormalBreast"))
annotation_col$Type<-as.factor(annotation_col$Type)
rownames(annotation_col) = c("Cell_A", "Cell_C", "Cell_E", "Cell_I", 
                             "FatLayer_B", "FatLayer_D", "FatLayer_F", "FatLayer_J", 
                             "NormalBreast_SRR2148259", "NormalBreast_SRR2148260", "NormalBreast_SRR2148261")
annotation_row = matrix(exp_deg$regulation,nrow=length(exp_deg$regulation), ncol=1,byrow=TRUE)
rownames(annotation_row) = rownames(exp_deg)
colnames(annotation_row)<-"gene_type"
annotation_row<-as.data.frame(annotation_row)
ann_colors = list(Type = c(NormalBreast = "blue", Cell = "red", Fat = "red"),
                  gene_type = c(downregulated = "blue", upregulated  = "red"))

pdf(file="DESeq2_median_of_ratios/heatmap_FatCell-breast.pdf",width = 6,height = 8)
pheatmap(rt,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 2819, gaps_col = 8, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
rt_top50<-rbind(rt[1:25,],rt[2820:2854,])
pdf(file="DESeq2_median_of_ratios/heatmap_FatCell-breast_selected_top50.pdf",width = 6,height = 8)
pheatmap(rt_top50,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 25, gaps_col = 8, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()

tpm_selected<-subset(tpm, select = c("A", "C", "E", "I", "B", "D", "F", "J", "SRR2148259.sra", "SRR2148260.sra", "SRR2148261"))
colnames(tpm_selected)<-colnames(rt)
selected_genes<-rownames(rt)
tpm_data<-tpm_selected[selected_genes,]
pdf(file="TPM/heatmap_FatCell-breast.pdf",width = 6,height = 8)
pheatmap(tpm_data,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 2819, gaps_col = 8, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
tpm_data_top50<-rbind(tpm_data[1:25,],tpm_data[2820:2854,])
pdf(file="TPM/heatmap_FatCell-breast_selected_top50.pdf",width = 6,height = 8)
pheatmap(tpm_data_top50,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 25, gaps_col = 8, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()

library(ggpubr)
library("DOSE")
library("org.Hs.eg.db")
library("topGO")
library("clusterProfiler")
library("pathview")
###up DEG functional enrichment analysis
deg_up<-subset(exp_deg, regulation == "upregulated")
up_genes <- rownames(deg_up)
test = bitr(up_genes, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
head(test,2)
go_BP <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1,
                  keyType = 'ENTREZID')
y <- setReadable(go_BP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
### The geneID column is translated to symbol
write.csv(y,"GOBP_enrich_up_FatCell-breast_readable.csv",row.names =FALSE)
#结果可视化
GOBP_up<-summary(go_BP)
GOBP_up<-GOBP_up[1:20,]
GOBP_up$yvalue<-GOBP_up$Count/907
GOBP_up$xvalue<-GOBP_up$Description
library(ggpubr)
pdf("GOBP_up_FatCell-breast.barplot.pdf",height=6,width = 12)
ggbarplot(GOBP_up, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("red","pink"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

GOBP_up$minuslog10Pval=-log10(GOBP_up$pvalue)
pdf("GOBP_up_FatCell-breast.barplot_log_transformed_pvalue.pdf",height=6,width = 12)
ggbarplot(GOBP_up, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "minuslog10Pval",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("pink","red"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

###down DEG functional enrichment analysis
deg_down<-subset(exp_deg, regulation == "downregulated")
down_genes <- rownames(deg_down)
test = bitr(down_genes, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
head(test,2)
go_BP <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1,
                  keyType = 'ENTREZID')
y <- setReadable(go_BP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
### The geneID column is translated to symbol
write.csv(y,"GOBP_enrich_down_FatCell-breast_readable.csv",row.names =FALSE)
GOBP_down<-summary(go_BP)
GOBP_down<-GOBP_down[1:20,]
GOBP_down$yvalue<-GOBP_down$Count/2009
GOBP_down$xvalue<-GOBP_down$Description
pdf("GOBP_down.barplot_FatCell-breast.pdf",height=6,width = 12)
ggbarplot(GOBP_down, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("blue","purple"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

GOBP_down$minuslog10Pval=-log10(GOBP_down$pvalue)
pdf("GOBP_down_FatCell-breast.barplot_log_transformed_pvalue.pdf",height=6,width = 12)
ggbarplot(GOBP_down, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "minuslog10Pval",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("purple","blue"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()


##############compare5:identification of DEGs between fatbreast vs cell tissue#########################################
setwd("/Users/liuzhe/Desktop/UCSF/milk/analysis1/six_comparision_addIandJ/lastedversion/FatBreast_Cell")
rawdata<-subset(allsamples_counts, select = c("A","C","E","I","B","D","F","J","SRR2148259.sra","SRR2148260.sra","SRR2148261"))
dim(rawdata)
#[1] 21993     11
gset <- rawdata[rowMeans(rawdata)>0,] # 剔除表达量低的基因
dim(gset)
#[1] 21967     11
coldata <- data.frame(
  condition = factor(c("Cell", "Cell", "Cell", "Cell",
                       "FatBreast", "FatBreast", "FatBreast", "FatBreast",
                       "FatBreast", "FatBreast", "FatBreast")))
str(coldata)
condition = factor(c("Cell", "Cell", "Cell", "Cell",
                     "FatBreast", "FatBreast", "FatBreast", "FatBreast",
                     "FatBreast", "FatBreast", "FatBreast"))
subject <- factor(c("ID1012","ID1013","ID1014","ID1012",
                    "ID1012","ID1013","ID1014","ID1012",
                    "B1","B2","B3"))
rownames(coldata)<-colnames(gset)
coldata$subject<-subject
# 注意在design中加上配对信息
data = apply(gset, 2, as.integer) ## DESeq2分析需要是整数
row.names(data)<-row.names(gset)
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design =  ~ condition) 
dds$condition<- relevel(dds$condition, ref = "Cell") 
dds <- DESeq(dds)
nrDEG_DESeq2 <- as.data.frame(results(dds))
rld <- rlog(dds)
# 这里我还提取了标准化后的表达矩阵，可以用于后续的热图绘制等等
normal_gset <- assay(rld) 
nrDEG_DESeq2 = nrDEG_DESeq2[order(nrDEG_DESeq2$log2FoldChange),] 
## 4.4定义差异基因
nrDEG <- nrDEG_DESeq2
nrDEG$Group = "notsignificant"
logFC_cutoff <- 2
nrDEG$Group[which( (nrDEG$padj < 0.01) & (nrDEG$log2FoldChange > logFC_cutoff) )] = "upregulated"
nrDEG$Group[which( (nrDEG$padj < 0.01) & (nrDEG$log2FoldChange < -logFC_cutoff) )] = "downregulated"
table(nrDEG$Group)
#downregulated notsignificant    upregulated 
#         187          21493            287 
write.csv(nrDEG,"DEGs_for_FatBreast-cell.csv",quote = F)
library("ggplot2")
deg<-nrDEG
head(deg)
deg$color<-ifelse(deg$padj<0.01 & abs(deg$log2FoldChange)>logFC_cutoff, ifelse(deg$log2FoldChange< -logFC_cutoff, "blue", "red"),"gray")
color<-c(red = "red", gray = "gray", blue = "blue")
deg$symbol<-row.names(deg)
p <- ggplot(data = deg, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=Group)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  theme_bw() +
  xlim(-30,30)
for_label <- deg %>% 
  filter(abs(log2FoldChange) >2& padj < 0.01)
pdf("VP_DEGs_FatBreast-cell.pdf")
p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )+ theme(legend.position="none")
dev.off()
#riskgenes = c("gene1", "gene2", "gene3", "gene4")
#for_label<-deg[riskgenes,]
normal_gset<-as.data.frame(normal_gset)
deg_down<-subset(deg, Group == "downregulated")
deg_up<-subset(deg, Group == "upregulated")
deg_select<-rbind(deg_down,deg_up)
normalized_exp<-normal_gset[rownames(deg_select),]
exp_deg<-cbind(deg_select,normalized_exp)
exp_deg<-subset(exp_deg, select = c("B", "D", "F", "J", "SRR2148259.sra", "SRR2148260.sra", "SRR2148261",
                                    "A", "C", "E", "I", "log2FoldChange", "pvalue", "padj", "Group"))
colnames(exp_deg)<-c("FatLayer_B",	"FatLayer_D",	"FatLayer_F",	"FatLayer_J",
                     "NormalBreast_SRR2148259", "NormalBreast_SRR2148260", "NormalBreast_SRR2148261",
                     "Cell_A", "Cell_C", "Cell_E", "Cell_I", 
                     "log2FoldChange", "raw_Pvalue","adjustedPvalue","Regulation")
#exp_deg<-cbind(round(exp_deg[1:10],3),exp_deg[,11])
colnames(exp_deg)[15]<-"regulation"
write.csv(exp_deg,"Table5_FatBreast-cell.csv",quote = F)
library("pheatmap")
rt<-exp_deg[,1:11]
annotation_col = data.frame(Type = c("Fat", "Fat", "Fat", "Fat", "NormalBreast", "NormalBreast", "NormalBreast",
                                     "Cell", "Cell", "Cell", "Cell"))
annotation_col$Type<-as.factor(annotation_col$Type)
rownames(annotation_col) = c("FatLayer_B", "FatLayer_D", "FatLayer_F", "FatLayer_J", 
                             "NormalBreast_SRR2148259", "NormalBreast_SRR2148260", "NormalBreast_SRR2148261",
                             "Cell_A", "Cell_C", "Cell_E", "Cell_I")
annotation_row = matrix(exp_deg$regulation,nrow=length(exp_deg$regulation), ncol=1,byrow=TRUE)
rownames(annotation_row) = rownames(exp_deg)
colnames(annotation_row)<-"gene_type"
annotation_row<-as.data.frame(annotation_row)
ann_colors = list(Type = c(NormalBreast = "red", Cell = "blue", Fat = "red"),
                  gene_type = c(downregulated = "blue", upregulated  = "red"))

pdf(file="DESeq2_median_of_ratios/heatmap_FatBreast-cell.pdf",width = 6,height = 8)
pheatmap(rt,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 187, gaps_col = 7, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
rt_top50<-rbind(rt[1:25,],rt[188:212,])
pdf(file="DESeq2_median_of_ratios/heatmap_FatBreast-cell_selected_top50.pdf",width = 6,height = 8)
pheatmap(rt_top50,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 25, gaps_col = 7, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()

tpm_selected<-subset(tpm, select = c("B", "D", "F", "J", "SRR2148259.sra", "SRR2148260.sra", "SRR2148261","A", "C", "E", "I"))
colnames(tpm_selected)<-colnames(rt)
selected_genes<-rownames(rt)
tpm_data<-tpm_selected[selected_genes,]
pdf(file="TPM/heatmap_FatBreast-cell.pdf",width = 6,height = 8)
pheatmap(tpm_data,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 187, gaps_col = 7, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
tpm_data_top50<-rbind(tpm_data[1:25,],tpm_data[188:212,])
pdf(file="TPM/heatmap_FatBreast-cell_selected_top50.pdf",width = 6,height = 8)
pheatmap(tpm_data_top50,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 25, gaps_col = 7, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()

library(ggpubr)
library("DOSE")
library("org.Hs.eg.db")
library("topGO")
library("clusterProfiler")
library("pathview")
###up DEG functional enrichment analysis
deg_up<-subset(exp_deg, regulation == "upregulated")
up_genes <- rownames(deg_up)
test = bitr(up_genes, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
head(test,2)
go_BP <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1,
                  keyType = 'ENTREZID')
y <- setReadable(go_BP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
### The geneID column is translated to symbol
write.csv(y,"GOBP_enrich_up_FatBreast-cell_readable.csv",row.names =FALSE)
#结果可视化
GOBP_up<-summary(go_BP)
GOBP_up<-GOBP_up[1:20,]
GOBP_up$yvalue<-GOBP_up$Count/217
GOBP_up$xvalue<-GOBP_up$Description
library(ggpubr)
pdf("GOBP_up_FatBreast-cell.barplot.pdf",height=6,width = 12)
ggbarplot(GOBP_up, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("red","pink"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

GOBP_up$minuslog10Pval=-log10(GOBP_up$pvalue)
pdf("GOBP_up_FatBreast-cell.barplot_log_transformed_pvalue.pdf",height=6,width = 12)
ggbarplot(GOBP_up, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "minuslog10Pval",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("pink","red"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

###down DEG functional enrichment analysis
deg_down<-subset(exp_deg, regulation == "downregulated")
down_genes <- rownames(deg_down)
test = bitr(down_genes, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
head(test,2)
go_BP <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1,
                  keyType = 'ENTREZID')
y <- setReadable(go_BP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
### The geneID column is translated to symbol
write.csv(y,"GOBP_enrich_down_FatBreast-cell_readable.csv",row.names =FALSE)
GOBP_down<-summary(go_BP)
GOBP_down<-GOBP_down[1:20,]
GOBP_down$yvalue<-GOBP_down$Count/166
GOBP_down$xvalue<-GOBP_down$Description
pdf("GOBP_down.barplot_FatBreast-cell.pdf",height=6,width = 12)
ggbarplot(GOBP_down, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("blue","purple"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

GOBP_down$minuslog10Pval=-log10(GOBP_down$pvalue)
pdf("GOBP_down_FatBreast-cell.barplot_log_transformed_pvalue.pdf",height=6,width = 12)
ggbarplot(GOBP_down, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "minuslog10Pval",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("purple","blue"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()



##############compare6:identification of DEGs between cellbreast vs fat breast tissue#########################################
setwd("/Users/liuzhe/Desktop/UCSF/milk/analysis1/six_comparision_addIandJ/lastedversion/CellBreast_Fat")
rawdata<-subset(allsamples_counts, select = c("A","C","E","I","B","D","F","J","SRR2148259.sra","SRR2148260.sra","SRR2148261"))
dim(rawdata)
#[1] 21993     11
gset <- rawdata[rowMeans(rawdata)>0,] # 剔除表达量低的基因
dim(gset)
#[1] 21967     11
coldata <- data.frame(
  condition = factor(c("CellBreast", "CellBreast", "CellBreast", "CellBreast",
                       "Fat", "Fat", "Fat", "Fat",
                       "CellBreast", "CellBreast", "CellBreast")))
str(coldata)
condition = factor(c("CellBreast", "CellBreast", "CellBreast", "CellBreast",
                     "Fat", "Fat", "Fat", "Fat",
                     "CellBreast", "CellBreast", "CellBreast"))
subject <- factor(c("ID1012","ID1013","ID1014","ID1012",
                    "ID1012","ID1013","ID1014","ID1012",
                    "B1","B2","B3"))
rownames(coldata)<-colnames(gset)
coldata$subject<-subject
# 注意在design中加上配对信息
data = apply(gset, 2, as.integer) ## DESeq2分析需要是整数
row.names(data)<-row.names(gset)
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design =  ~ condition) 
dds$condition<- relevel(dds$condition, ref = "Fat") 
dds <- DESeq(dds)
nrDEG_DESeq2 <- as.data.frame(results(dds))
rld <- rlog(dds)
# 这里我还提取了标准化后的表达矩阵，可以用于后续的热图绘制等等
normal_gset <- assay(rld) 
nrDEG_DESeq2 = nrDEG_DESeq2[order(nrDEG_DESeq2$log2FoldChange),] 
## 4.4定义差异基因
nrDEG <- nrDEG_DESeq2
nrDEG$Group = "notsignificant"
logFC_cutoff <- 2
nrDEG$Group[which( (nrDEG$padj < 0.01) & (nrDEG$log2FoldChange > logFC_cutoff) )] = "upregulated"
nrDEG$Group[which( (nrDEG$padj < 0.01) & (nrDEG$log2FoldChange < -logFC_cutoff) )] = "downregulated"
table(nrDEG$Group)
#downregulated notsignificant    upregulated 
#           25          20220           1722 
write.csv(nrDEG,"DEGs_for_CellBreast-fat.csv",quote = F)
library("ggplot2")
deg<-nrDEG
head(deg)
deg$color<-ifelse(deg$padj<0.01 & abs(deg$log2FoldChange)>logFC_cutoff, ifelse(deg$log2FoldChange< -logFC_cutoff, "blue", "red"),"gray")
color<-c(red = "red", gray = "gray", blue = "blue")
deg$symbol<-row.names(deg)
p <- ggplot(data = deg, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=Group)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  theme_bw() +
  xlim(-30,30)
for_label <- deg %>% 
  filter(abs(log2FoldChange) >1& padj < 0.05)
pdf("VP_DEGs_CellBreast-fat.pdf")
p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )+ theme(legend.position="none")
dev.off()
#riskgenes = c("gene1", "gene2", "gene3", "gene4")
#for_label<-deg[riskgenes,]
normal_gset<-as.data.frame(normal_gset)
deg_down<-subset(deg, Group == "downregulated")
deg_up<-subset(deg, Group == "upregulated")
deg_select<-rbind(deg_down,deg_up)
normalized_exp<-normal_gset[rownames(deg_select),]
exp_deg<-cbind(deg_select,normalized_exp)
exp_deg<-subset(exp_deg, select = c("A", "C", "E", "I","SRR2148259.sra", "SRR2148260.sra", "SRR2148261",
                                    "B", "D", "F", "J", "log2FoldChange", "pvalue", "padj", "Group"))
colnames(exp_deg)<-c("Cell_A", "Cell_C", "Cell_E", "Cell_I",
                     "NormalBreast_SRR2148259", "NormalBreast_SRR2148260", "NormalBreast_SRR2148261",
                     "FatLayer_B",	"FatLayer_D",	"FatLayer_F",	"FatLayer_J",
                     "log2FoldChange", "raw_Pvalue","adjustedPvalue","Regulation")
#exp_deg<-cbind(round(exp_deg[1:10],3),exp_deg[,11])
colnames(exp_deg)[15]<-"regulation"
write.csv(exp_deg,"Table6_CellBreast-fat.csv",quote = F)
library("pheatmap")
rt<-exp_deg[,1:11]
annotation_col = data.frame(Type = c("Cell", "Cell", "Cell", "Cell","NormalBreast", "NormalBreast", "NormalBreast",
                                     "Fat", "Fat", "Fat", "Fat"))
annotation_col$Type<-as.factor(annotation_col$Type)
rownames(annotation_col) = c("Cell_A", "Cell_C", "Cell_E", "Cell_I",
                             "NormalBreast_SRR2148259", "NormalBreast_SRR2148260", "NormalBreast_SRR2148261",
                             "FatLayer_B", "FatLayer_D", "FatLayer_F", "FatLayer_J")
annotation_row = matrix(exp_deg$regulation,nrow=length(exp_deg$regulation), ncol=1,byrow=TRUE)
rownames(annotation_row) = rownames(exp_deg)
colnames(annotation_row)<-"gene_type"
annotation_row<-as.data.frame(annotation_row)
ann_colors = list(Type = c(NormalBreast = "red", Cell = "red", Fat = "blue"),
                  gene_type = c(downregulated = "blue", upregulated  = "red"))

pdf(file="DESeq2_median_of_ratios/heatmap_CellBreast-fat.pdf",width = 6,height = 8)
pheatmap(rt,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 25, gaps_col = 7, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
rt_top50<-rbind(rt[1:25,],rt[26:50,])
pdf(file="DESeq2_median_of_ratios/heatmap_CellBreast-fat_selected_top50.pdf",width = 6,height = 8)
pheatmap(rt_top50,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 25, gaps_col = 7, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()

tpm_selected<-subset(tpm, select = c("A", "C", "E", "I", "SRR2148259.sra", "SRR2148260.sra", "SRR2148261", "B", "D", "F", "J"))
colnames(tpm_selected)<-colnames(rt)
selected_genes<-rownames(rt)
tpm_data<-tpm_selected[selected_genes,]
pdf(file="TPM/heatmap_CellBreast-fat.pdf",width = 6,height = 8)
pheatmap(tpm_data,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 25, gaps_col = 7, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
tpm_data_top50<-rbind(tpm_data[1:25,],tpm_data[26:50,])
pdf(file="TPM/heatmap_CellBreast-fat_selected_top50.pdf",width = 6,height = 8)
pheatmap(tpm_data_top50,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 25, gaps_col = 7, angle_col = "45",
         annotation_colors = ann_colors,  show_rownames = T, show_colnames = T, main = "Title", fontsize = 6)
#fontsize = 1
dev.off()

library(ggpubr)
library("DOSE")
library("org.Hs.eg.db")
library("topGO")
library("clusterProfiler")
library("pathview")
###up DEG functional enrichment analysis
deg_up<-subset(exp_deg, regulation == "upregulated")
up_genes <- rownames(deg_up)
test = bitr(up_genes, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
head(test,2)
go_BP <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1,
                  keyType = 'ENTREZID')
y <- setReadable(go_BP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
### The geneID column is translated to symbol
write.csv(y,"GOBP_enrich_up_CellBreast-fat_readable.csv",row.names =FALSE)
#结果可视化
GOBP_up<-summary(go_BP)
GOBP_up<-GOBP_up[1:20,]
GOBP_up$yvalue<-GOBP_up$Count/1300
GOBP_up$xvalue<-GOBP_up$Description
library(ggpubr)
pdf("GOBP_up_CellBreast-fat.barplot.pdf",height=6,width = 12)
ggbarplot(GOBP_up, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("red","pink"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

GOBP_up$minuslog10Pval=-log10(GOBP_up$pvalue)
pdf("GOBP_up_CellBreast-fat.barplot_log_transformed_pvalue.pdf",height=6,width = 12)
ggbarplot(GOBP_up, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "minuslog10Pval",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("pink","red"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

###down DEG functional enrichment analysis
deg_down<-subset(exp_deg, regulation == "downregulated")
down_genes <- rownames(deg_down)
test = bitr(down_genes, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
head(test,2)
go_BP <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1,
                  keyType = 'ENTREZID')
y <- setReadable(go_BP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
### The geneID column is translated to symbol
write.csv(y,"GOBP_enrich_down_CellBreast-fat_readable.csv",row.names =FALSE)
GOBP_down<-summary(go_BP)
GOBP_down<-GOBP_down[1:20,]
GOBP_down$yvalue<-GOBP_down$Count/15
GOBP_down$xvalue<-GOBP_down$Description
pdf("GOBP_down.barplot_CellBreast-fat.pdf",height=6,width = 12)
ggbarplot(GOBP_down, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("blue","purple"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

GOBP_down$minuslog10Pval=-log10(GOBP_down$pvalue)
pdf("GOBP_down_CellBreast-fat.barplot_log_transformed_pvalue.pdf",height=6,width = 12)
ggbarplot(GOBP_down, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "minuslog10Pval",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("purple","blue"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

