setwd("/Users/liuzhe/Desktop/UCSF/milk/analysis2/six_comparsion")
rm(list=ls())

#一 载入R包，数据
# 加载包
library("DESeq2")
library("limma")
library("edgeR")
library("dplyr") 
library("biomaRt")
library("tidyverse")
library("stringr")


setwd("/Users/liuzhe/Desktop/UCSF/milk/analysis2/six_comparsion")
rawdata <- read.table("/Users/liuzhe/Desktop/UCSF/milk/analysis2/datafile/Partek_Supply_RNAseq_Quantify_to_annotation_model_(Partek_E_M)_Gene_counts.txt",sep="\t",row.names = 1)
dim(rawdata)
#[1]    15 21736
rawdata<-rawdata[1:15,4:21736]
rawdata<-t(rawdata)
class(rawdata)
#[1] "matrix" "array" 
dim(rawdata)
#[1] 21733    15
rawdata<-data.frame(rawdata)
sampleA<-read.table("/Users/liuzhe/Desktop/UCSF/milk/analysis1/addSample_IJ/A_featurecounts.txt", header = T)
mydict<-subset(sampleA, select = c("Geneid", "Length"))
rawdata$gene_name<-rawdata$Sample.name
merged_data<-merge(rawdata, mydict, by.x = "gene_name", by.y = "Geneid")
merged_data<-subset(merged_data, select = c("gene_name", "Length", "ID1016.V1_S0", "ID1020.V3_S148", "ID1020.V2_S146",
                  "ID1018.V2_S0", "ID1020.V5_S149", "ID1028.V1_S0", "ID1020.V1_S145", "ID1023.V2_S140",
                  "ID1026.V1_S137", "VPL01.131_S150", "ID1016.V2_S0", "ID1024.V3_S139", "ID1025.V1_S138", "ID1029.V1_S0"))
kb<-merged_data$Length/1000
rpk <- apply(merged_data[,3:16], 2, as.numeric)/ kb
rownames(rpk)<-merged_data$gene_name
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
write.csv(tpm,file="/Users/liuzhe/Desktop/UCSF/milk/analysis2/six_comparsion/tpm.csv",quote=F)
write.csv(rawdata, "/Users/liuzhe/Desktop/UCSF/milk/analysis2/six_comparsion/allsamples_readcounts_anno.csv",quote=F,row.names = F)


#1.读入原始数据及分组信息
allsamples_counts<-read.csv("/Users/liuzhe/Desktop/UCSF/milk/analysis2/datafile/readcount_all.csv", row.name = 1)
metadata<-read.table("/Users/liuzhe/Desktop/UCSF/milk/analysis2/datafile/metadata.txt", sep = "\t", header = T)

rawdata<-allsamples_counts
dim(rawdata)
#[1] 21733    14
gset <- rawdata[rowMeans(rawdata)>0,] # 剔除表达量低的基因
dim(gset)
#[1] 21733     14

metadata_high<-subset(metadata, Supply == "High")
metadata_low<-subset(metadata, Supply == "Low")
metadata_normal<-subset(metadata, Supply == "Normal")

######################compare1: high -- low##########################################################
setwd("/Users/liuzhe/Desktop/UCSF/milk/analysis2/six_comparsion/High_Low")
sample_selected<-c(metadata_high$Sample.name, metadata_low$Sample.name)
sample_selected<-gsub("-", "\\.", sample_selected)
counts_anno_high_low<-subset(gset, select = sample_selected)
group_list<-c("High", "High", "High", "Low", "Low", "Low", "Low")
group_list <- factor(group_list,levels = c("Low", "High"))
## 表达矩阵
data = apply(counts_anno_high_low, 2, as.integer) ## DESeq2分析需要是整数
row.names(data) <- row.names(counts_anno_high_low)
## 4.2分组矩阵，配对分析与常规分析最大的区别就在分组矩阵
condition = group_list
# 配对分析要加上这段代码，知道谁和谁是一对，比如1,1是一对，5,5是一对
coldata <- data.frame(row.names = colnames(data), condition)
# 注意在design中加上配对信息
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design = ~condition) 
dds$condition<- relevel(dds$condition, ref = "Low") 
## 4.3差异表达矩阵，还是和常规分析一样
dds <- DESeq(dds)
nrDEG_DESeq2 <- as.data.frame(results(dds))
rld <- rlog(dds)
# 这里我还提取了标准化后的表达矩阵，可以用于后续的热图绘制等等
normal_gset <- assay(rld) 
nrDEG_DESeq2 = nrDEG_DESeq2[order(nrDEG_DESeq2$log2FoldChange),] 
## 4.4定义差异基因
nrDEG <- nrDEG_DESeq2
nrDEG$Group = "notsignificant"
logFC_cutoff <- 1
nrDEG$Group[which( (nrDEG$padj < 0.05) & (nrDEG$log2FoldChange > logFC_cutoff) )] = "upregulated"
nrDEG$Group[which( (nrDEG$padj < 0.05) & (nrDEG$log2FoldChange < -logFC_cutoff) )] = "downregulated"
table(nrDEG$Group)
#downregulated notsignificant    upregulated 
#           24          21694             15
write.csv(nrDEG,"DEGs_for_High-Low.csv",quote = F)
library("ggplot2")
deg<-nrDEG
head(deg)
deg$color<-ifelse(deg$padj<0.05 & abs(deg$log2FoldChange)>logFC_cutoff, ifelse(deg$log2FoldChange< -logFC_cutoff, "blue", "red"),"gray")
color<-c(red = "red", gray = "gray", blue = "blue")
deg$symbol<-row.names(deg)
p <- ggplot(data = deg, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=Group)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw() +
  xlim(-10,10)
for_label <- deg %>% 
  filter(abs(log2FoldChange) >1& padj < 0.05)
pdf("VP_DEGs_High-Low.pdf")
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
exp_deg<-subset(exp_deg, select = c("ID1024.V3_S139", "ID1025.V1_S138", "ID1026.V1_S137", "ID1018.V2_S0", 
                                    "ID1023.V2_S140", "ID1028.V1_S0", "ID1029.V1_S0", "log2FoldChange",
                                    "pvalue","padj","Group"))
colnames(exp_deg)<-c("High_ID1024.V3_S139", "High_ID1025.V1_S138", "High_ID1026.V1_S137", "Low_ID1018.V2_S0", 
                     "Low_ID1023.V2_S140", "Low_ID1028.V1_S0", "Low_ID1029.V1_S0","log2FoldChange",
                     "raw_Pvalue","adjustedPvalue","Regulation" )
exp_deg<-cbind(round(exp_deg[1:10],3),exp_deg[,11])
colnames(exp_deg)[11]<-"regulation"
write.csv(exp_deg,"Table1_High-Low.csv",quote = F)
library("pheatmap")
rt<-exp_deg[,1:7]
metadata_high_low<-rbind(metadata_high, metadata_low)
annotation_col = data.frame(Type = factor(metadata_high_low$Supply))
rownames(annotation_col) = colnames(exp_deg)[1:7]
annotation_row = matrix(exp_deg$regulation,nrow=length(exp_deg$regulation), ncol=1,byrow=TRUE)
rownames(annotation_row) = rownames(exp_deg)
colnames(annotation_row)<-"gene_type"
ann_colors = list(Type = c("Low" = "blue", "High" = "red"))
pdf(file="DESeq/heatmap_High-Low.pdf",width = 6,height = 8)
pheatmap(rt,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = T,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         gaps_row = 24,gaps_col = 3,angle_col = "45", annotation_colors = ann_colors, 
         main = "Title", fontsize = 6)
#fontsize = 1
dev.off()

tpm_selected<-subset(tpm, select = c("ID1024.V3_S139", "ID1025.V1_S138", "ID1026.V1_S137", "ID1018.V2_S0", 
                                     "ID1023.V2_S140", "ID1028.V1_S0", "ID1029.V1_S0"))
colnames(tpm_selected)<-colnames(rt)
selected_genes<-rownames(rt)
tpm_data<-tpm_selected[selected_genes,]
pdf(file="TPM/heatmap_High-Low.pdf",width = 6,height = 8)
pheatmap(tpm_data,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = T,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         gaps_row = 24,gaps_col = 3,angle_col = "45", annotation_colors = ann_colors, 
         main = "Title", fontsize = 6)
dev.off()
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
write.csv(y,"GOBP_enrich_up_High-Low.csv",row.names =FALSE)
#结果可视化
GOBP_up<-summary(go_BP)
GOBP_up<-GOBP_up[1:20,]
GOBP_up$yvalue<-GOBP_up$Count/14
GOBP_up$xvalue<-GOBP_up$Description
library(ggpubr)
pdf("GOBP_up_High-Low.barplot.pdf",height=6,width = 12)
ggbarplot(GOBP_up, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("red","pink"))+
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
write.csv(y,"GOBP_enrich_down_High-Low.csv",row.names =FALSE)
#结果可视化
GOBP_down<-summary(go_BP)
GOBP_down<-GOBP_down[1:20,]
GOBP_down$yvalue<-GOBP_down$Count/19
GOBP_down$xvalue<-GOBP_down$Description
pdf("GOBP_down.barplot_High-Low.pdf",height=6,width = 12)
ggbarplot(GOBP_down, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("blue","purple"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()


######################compare2: high -- normal##########################################################
setwd("/Users/liuzhe/Desktop/UCSF/milk/analysis2/six_comparsion/High_Normal/")
sample_selected<-c(metadata_high$Sample.name, metadata_normal$Sample.name)
sample_selected<-gsub("-", "\\.", sample_selected)
counts_anno_high_normal<-subset(gset, select = sample_selected)
group_list<-c("High", "High", "High", "Normal", "Normal", "Normal", "Normal", "Normal", "Normal", "Normal")
group_list <- factor(group_list,levels = c("Normal", "High"))
## 表达矩阵
data = apply(counts_anno_high_low, 2, as.integer) ## DESeq2分析需要是整数
row.names(data) <- row.names(counts_anno_high_normal)
## 4.2分组矩阵，配对分析与常规分析最大的区别就在分组矩阵
condition = group_list
# 配对分析要加上这段代码，知道谁和谁是一对，比如1,1是一对，5,5是一对
coldata <- data.frame(row.names = colnames(data), condition)
# 注意在design中加上配对信息
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design = ~condition) 
dds$condition<- relevel(dds$condition, ref = "Normal") 
## 4.3差异表达矩阵，还是和常规分析一样
dds <- DESeq(dds)
nrDEG_DESeq2 <- as.data.frame(results(dds))
rld <- rlog(dds)
# 这里我还提取了标准化后的表达矩阵，可以用于后续的热图绘制等等
normal_gset <- assay(rld) 
nrDEG_DESeq2 = nrDEG_DESeq2[order(nrDEG_DESeq2$log2FoldChange),] 
## 4.4定义差异基因
nrDEG <- nrDEG_DESeq2
nrDEG$Group = "notsignificant"
logFC_cutoff <- 1
nrDEG$Group[which( (nrDEG$padj < 0.05) & (nrDEG$log2FoldChange > logFC_cutoff) )] = "upregulated"
nrDEG$Group[which( (nrDEG$padj < 0.05) & (nrDEG$log2FoldChange < -logFC_cutoff) )] = "downregulated"
table(nrDEG$Group)
#downregulated notsignificant    upregulated 
#           78          21610             45 
write.csv(nrDEG,"DEGs_for_High-Normal.csv",quote = F)
library("ggplot2")
deg<-nrDEG
head(deg)
deg$color<-ifelse(deg$padj<0.05 & abs(deg$log2FoldChange)>logFC_cutoff, ifelse(deg$log2FoldChange< -logFC_cutoff, "blue", "red"),"gray")
color<-c(red = "red", gray = "gray", blue = "blue")
deg$symbol<-row.names(deg)
p <- ggplot(data = deg, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=Group)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw() +
  xlim(-10,10)
for_label <- deg %>% 
  filter(abs(log2FoldChange) >1& padj < 0.05)
pdf("VP_DEGs_High-Normal.pdf")
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
exp_deg<-subset(exp_deg, select = c("ID1024.V3_S139", "ID1025.V1_S138", "ID1026.V1_S137", "ID1016.V1_S0", "ID1016.V2_S0",
                                    "ID1020.V1_S145", "ID1020.V2_S146", "ID1020.V3_S148", "ID1020.V5_S149", "VPL01.131_S150",
                                    "log2FoldChange", "pvalue","padj","Group"))
colnames(exp_deg)<-c("High_ID1024.V3_S139", "High_ID1025.V1_S138", "High_ID1026.V1_S137", 
                     "Normal_ID1016.V1_S0", "Normal_ID1016.V2_S0", "Normal_ID1020.V1_S145", "Normal_ID1020.V2_S146",
                     "Normal_ID1020.V3_S148", "Normal_ID1020.V5_S149", "Normal_VPL01.131_S150", "log2FoldChange",
                     "raw_Pvalue","adjustedPvalue","Regulation" )
exp_deg<-cbind(round(exp_deg[1:13],3),exp_deg[,14])
colnames(exp_deg)[14]<-"regulation"
write.csv(exp_deg,"Table2_High-Normal.csv",quote = F)
library("pheatmap")
rt<-exp_deg[,1:10]
metadata_high_normal<-rbind(metadata_high, metadata_normal)
annotation_col = data.frame(Type = factor(metadata_high_normal$Supply))
rownames(annotation_col) = colnames(exp_deg)[1:10]
annotation_row = matrix(exp_deg$regulation,nrow=length(exp_deg$regulation), ncol=1,byrow=TRUE)
rownames(annotation_row) = rownames(exp_deg)
colnames(annotation_row)<-"gene_type"
ann_colors = list(Type = c("Normal" = "blue", "High" = "red"))
pdf(file="DESeq2/heatmap_High-Normal.pdf",width = 6,height = 8)
pheatmap(rt,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = T,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         gaps_row = 78,gaps_col = 3,angle_col = "45", annotation_colors = ann_colors, 
         main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
tpm$geneid<-rownames(tpm)
tpm$geneid<-gsub("-","\\.",tpm$geneid)
rownames(tpm)<-tpm$geneid
tpm$geneid<-NULL
tpm_selected<-subset(tpm, select = c("ID1024.V3_S139", "ID1025.V1_S138", "ID1026.V1_S137", "ID1016.V1_S0",
                                     "ID1016.V2_S0", "ID1020.V1_S145", "ID1020.V2_S146",
                                     "ID1020.V3_S148", "ID1020.V5_S149", "VPL01.131_S150"))
colnames(tpm_selected)<-colnames(rt)
selected_genes<-rownames(rt)
tpm_data<-tpm_selected[selected_genes,]
pdf(file="TPM/heatmap_High-Normal.pdf",width = 6,height = 8)
pheatmap(tpm_data,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         gaps_row = 78,gaps_col = 3,angle_col = "45", annotation_colors = ann_colors, 
         main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
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
write.csv(y,"GOBP_enrich_up_High-Normal.csv",row.names =FALSE)
#结果可视化
GOBP_up<-summary(go_BP)
GOBP_up<-GOBP_up[1:20,]
GOBP_up$yvalue<-GOBP_up$Count/41
GOBP_up$xvalue<-GOBP_up$Description
library(ggpubr)
pdf("GOBP_up_High-Normal.barplot.pdf",height=6,width = 12)
ggbarplot(GOBP_up, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("red","pink"))+
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
write.csv(y,"GOBP_enrich_down_High-Normal.csv",row.names =FALSE)
#结果可视化
GOBP_down<-summary(go_BP)
GOBP_down<-GOBP_down[1:20,]
GOBP_down$yvalue<-GOBP_down$Count/55
GOBP_down$xvalue<-GOBP_down$Description
pdf("GOBP_down.barplot_High-Normal.pdf",height=6,width = 12)
ggbarplot(GOBP_down, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("blue","purple"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

######################compare3: low -- normal##########################################################
setwd("/Users/liuzhe/Desktop/UCSF/milk/analysis2/six_comparsion/Low_Normal/")
sample_selected<-c(metadata_low$Sample.name, metadata_normal$Sample.name)
sample_selected<-gsub("-", "\\.", sample_selected)
counts_anno_low_normal<-subset(gset, select = sample_selected)
group_list<-c("Low", "Low", "Low", "Low", "Normal", "Normal", "Normal", "Normal", "Normal", "Normal", "Normal")
group_list <- factor(group_list,levels = c("Normal", "Low"))
## 表达矩阵
data = apply(counts_anno_low_normal, 2, as.integer) ## DESeq2分析需要是整数
row.names(data) <- row.names(counts_anno_low_normal)
## 4.2分组矩阵，配对分析与常规分析最大的区别就在分组矩阵
condition = group_list
# 配对分析要加上这段代码，知道谁和谁是一对，比如1,1是一对，5,5是一对
coldata <- data.frame(row.names = colnames(data), condition)
# 注意在design中加上配对信息
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design = ~condition) 
dds$condition<- relevel(dds$condition, ref = "Normal") 
## 4.3差异表达矩阵，还是和常规分析一样
dds <- DESeq(dds)
nrDEG_DESeq2 <- as.data.frame(results(dds))
rld <- rlog(dds)
# 这里我还提取了标准化后的表达矩阵，可以用于后续的热图绘制等等
normal_gset <- assay(rld) 
nrDEG_DESeq2 = nrDEG_DESeq2[order(nrDEG_DESeq2$log2FoldChange),] 
## 4.4定义差异基因
nrDEG <- nrDEG_DESeq2
nrDEG$Group = "notsignificant"
logFC_cutoff <- 1
nrDEG$Group[which( (nrDEG$padj < 0.05) & (nrDEG$log2FoldChange > logFC_cutoff) )] = "upregulated"
nrDEG$Group[which( (nrDEG$padj < 0.05) & (nrDEG$log2FoldChange < -logFC_cutoff) )] = "downregulated"
table(nrDEG$Group)
#downregulated notsignificant    upregulated 
#        26          21684             23
write.csv(nrDEG,"DEGs_for_Low-Normal.csv",quote = F)
library("ggplot2")
deg<-nrDEG
head(deg)
deg$color<-ifelse(deg$padj<0.05 & abs(deg$log2FoldChange)>logFC_cutoff, ifelse(deg$log2FoldChange< -logFC_cutoff, "blue", "red"),"gray")
color<-c(red = "red", gray = "gray", blue = "blue")
deg$symbol<-row.names(deg)
p <- ggplot(data = deg, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=Group)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw() +
  xlim(-10,10)
for_label <- deg %>% 
  filter(abs(log2FoldChange) >1& padj < 0.05)
pdf("VP_DEGs_Low-Normal.pdf")
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
exp_deg<-subset(exp_deg, select = c("ID1018.V2_S0", "ID1023.V2_S140", "ID1028.V1_S0", "ID1029.V1_S0",
                                    "ID1016.V1_S0", "ID1016.V2_S0", "ID1020.V1_S145", "ID1020.V2_S146",
                                    "ID1020.V3_S148", "ID1020.V5_S149", "VPL01.131_S150",
                                    "log2FoldChange", "pvalue","padj","Group"))
colnames(exp_deg)<-c("Low_ID1018.V2_S0", "Low_ID1023.V2_S140", "Low_ID1028.V1_S0", "Low_ID1029.V1_S0",
                     "Normal_ID1016.V1_S0", "Normal_ID1016.V2_S0", "Normal_ID1020.V1_S145", "Normal_ID1020.V2_S146",
                     "Normal_ID1020.V3_S148", "Normal_ID1020.V5_S149", "Normal_VPL01.131_S150", "log2FoldChange",
                     "raw_Pvalue","adjustedPvalue","Regulation" )
exp_deg<-cbind(round(exp_deg[1:14],3),exp_deg[,15])
colnames(exp_deg)[15]<-"regulation"
write.csv(exp_deg,"Table3_Low-Normal.csv",quote = F)
library("pheatmap")
rt<-exp_deg[,1:11]
metadata_low_normal<-rbind(metadata_low, metadata_normal)
annotation_col = data.frame(Type = factor(metadata_low_normal$Supply))
rownames(annotation_col) = colnames(exp_deg)[1:11]
annotation_row = matrix(exp_deg$regulation,nrow=length(exp_deg$regulation), ncol=1,byrow=TRUE)
rownames(annotation_row) = rownames(exp_deg)
colnames(annotation_row)<-"gene_type"
ann_colors = list(Type = c("Normal" = "blue", "Low" = "red"))
pdf(file="DESeq2/heatmap_Low-Normal.pdf",width = 6,height = 8)
pheatmap(rt,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = T,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         gaps_row = 26,gaps_col = 4,angle_col = "45", annotation_colors = ann_colors, 
         main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
tpm$geneid<-rownames(tpm)
tpm$geneid<-gsub("-","\\.",tpm$geneid)
rownames(tpm)<-tpm$geneid
tpm$geneid<-NULL
tpm_selected<-subset(tpm, select = c("ID1018.V2_S0", "ID1023.V2_S140", "ID1028.V1_S0", "ID1029.V1_S0",
                                     "ID1016.V1_S0", "ID1016.V2_S0", "ID1020.V1_S145", "ID1020.V2_S146",
                                     "ID1020.V3_S148", "ID1020.V5_S149", "VPL01.131_S150"))
colnames(tpm_selected)<-colnames(rt)
selected_genes<-rownames(rt)
tpm_data<-tpm_selected[selected_genes,]
pdf(file="TPM/heatmap_Low-Normal.pdf",width = 6,height = 8)
pheatmap(tpm_data,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         gaps_row = 26,gaps_col = 4,angle_col = "45", annotation_colors = ann_colors, 
         main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
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
write.csv(y,"GOBP_enrich_up_Low-Normal.csv",row.names =FALSE)
#结果可视化
GOBP_up<-summary(go_BP)
GOBP_up<-GOBP_up[1:20,]
GOBP_up$yvalue<-GOBP_up$Count/18
GOBP_up$xvalue<-GOBP_up$Description
library(ggpubr)
pdf("GOBP_up_Low-Normal.barplot.pdf",height=6,width = 12)
ggbarplot(GOBP_up, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("red","pink"))+
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
write.csv(y,"GOBP_enrich_down_Low-Normal.csv",row.names =FALSE)
#结果可视化
GOBP_down<-summary(go_BP)
GOBP_down<-GOBP_down[1:20,]
GOBP_down$yvalue<-GOBP_down$Count/8
GOBP_down$xvalue<-GOBP_down$Description
pdf("GOBP_down.barplot_Low-Normal.pdf",height=6,width = 12)
ggbarplot(GOBP_down, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("blue","purple"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()

######################compare4: highlow -- normal##########################################################
setwd("/Users/liuzhe/Desktop/UCSF/milk/analysis2/six_comparsion/HighLow_Normal/")
sample_selected<-c(metadata_high$Sample.name, metadata_low$Sample.name, metadata_normal$Sample.name)
sample_selected<-gsub("-", "\\.", sample_selected)
counts_anno_highlow_normal<-subset(gset, select = sample_selected)
group_list<-c("HighLow", "HighLow", "HighLow",
              "HighLow", "HighLow", "HighLow", "HighLow",
              "Normal", "Normal", "Normal", "Normal", "Normal", "Normal", "Normal")
group_list <- factor(group_list,levels = c("HighLow", "Normal"))
## 表达矩阵
data = apply(counts_anno_highlow_normal, 2, as.integer) ## DESeq2分析需要是整数
row.names(data) <- row.names(counts_anno_highlow_normal)
## 4.2分组矩阵，配对分析与常规分析最大的区别就在分组矩阵
condition = group_list
# 配对分析要加上这段代码，知道谁和谁是一对，比如1,1是一对，5,5是一对
coldata <- data.frame(row.names = colnames(data), condition)
# 注意在design中加上配对信息
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design = ~condition) 
dds$condition<- relevel(dds$condition, ref = "Normal") 
## 4.3差异表达矩阵，还是和常规分析一样
dds <- DESeq(dds)
nrDEG_DESeq2 <- as.data.frame(results(dds))
rld <- rlog(dds)
# 这里我还提取了标准化后的表达矩阵，可以用于后续的热图绘制等等
normal_gset <- assay(rld) 
nrDEG_DESeq2 = nrDEG_DESeq2[order(nrDEG_DESeq2$log2FoldChange),] 
## 4.4定义差异基因
nrDEG <- nrDEG_DESeq2
nrDEG$Group = "notsignificant"
logFC_cutoff <- 1
nrDEG$Group[which( (nrDEG$padj < 0.05) & (nrDEG$log2FoldChange > logFC_cutoff) )] = "upregulated"
nrDEG$Group[which( (nrDEG$padj < 0.05) & (nrDEG$log2FoldChange < -logFC_cutoff) )] = "downregulated"
table(nrDEG$Group)
#downregulated notsignificant    upregulated 
#          14          21706             13 
write.csv(nrDEG,"DEGs_for_HighLow-Normal.csv",quote = F)
library("ggplot2")
deg<-nrDEG
head(deg)
deg$color<-ifelse(deg$padj<0.05 & abs(deg$log2FoldChange)>logFC_cutoff, ifelse(deg$log2FoldChange< -logFC_cutoff, "blue", "red"),"gray")
color<-c(red = "red", gray = "gray", blue = "blue")
deg$symbol<-row.names(deg)
p <- ggplot(data = deg, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=Group)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw() +
  xlim(-10,10)
for_label <- deg %>% 
  filter(abs(log2FoldChange) >1& padj < 0.05)
pdf("VP_DEGs_HighLow-Normal.pdf")
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
exp_deg<-subset(exp_deg, select = c("ID1024.V3_S139", "ID1025.V1_S138", "ID1026.V1_S137", 
                                    "ID1018.V2_S0", "ID1023.V2_S140", "ID1028.V1_S0", "ID1029.V1_S0",
                                    "ID1016.V1_S0", "ID1016.V2_S0", "ID1020.V1_S145", "ID1020.V2_S146",
                                    "ID1020.V3_S148", "ID1020.V5_S149", "VPL01.131_S150",
                                    "log2FoldChange", "pvalue","padj","Group"))
colnames(exp_deg)<-c("High_ID1024.V3_S139", "High_ID1025.V1_S138", "High_ID1026.V1_S137", 
                    "Low_ID1018.V2_S0", "Low_ID1023.V2_S140", "Low_ID1028.V1_S0", "Low_ID1029.V1_S0",
                     "Normal_ID1016.V1_S0", "Normal_ID1016.V2_S0", "Normal_ID1020.V1_S145", "Normal_ID1020.V2_S146",
                     "Normal_ID1020.V3_S148", "Normal_ID1020.V5_S149", "Normal_VPL01.131_S150", "log2FoldChange",
                     "raw_Pvalue","adjustedPvalue","Regulation" )
exp_deg<-cbind(round(exp_deg[1:17],3),exp_deg[,18])
colnames(exp_deg)[18]<-"regulation"
write.csv(exp_deg,"Table4_HighLow-Normal.csv",quote = F)
library("pheatmap")
rt<-exp_deg[,1:14]
metadata_highlow_normal<-rbind(metadata_high,metadata_low, metadata_normal)
annotation_col = data.frame(Type = factor(metadata_highlow_normal$Supply))
rownames(annotation_col) = colnames(exp_deg)[1:14]
annotation_row = matrix(exp_deg$regulation,nrow=length(exp_deg$regulation), ncol=1,byrow=TRUE)
rownames(annotation_row) = rownames(exp_deg)
colnames(annotation_row)<-"gene_type"
ann_colors = list(Type = c("Normal" = "blue", "Low" = "orange", "High" = "red"))
pdf(file="DESeq2/heatmap_HighLow-Normal.pdf",width = 6,height = 8)
pheatmap(rt,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         gaps_row = 14,gaps_col = 7,angle_col = "45", annotation_colors = ann_colors, 
         main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
tpm$geneid<-rownames(tpm)
tpm$geneid<-gsub("-","\\.",tpm$geneid)
rownames(tpm)<-tpm$geneid
tpm$geneid<-NULL
tpm_selected<-subset(tpm, select = c("ID1016.V1_S0", "ID1020.V3_S148", "ID1020.V2_S146",
                                      "ID1018.V2_S0", "ID1023.V2_S140", "ID1028.V1_S0", "ID1029.V1_S0",
                                     "ID1016.V1_S0", "ID1016.V2_S0", "ID1020.V1_S145", "ID1020.V2_S146",
                                     "ID1020.V3_S148", "ID1020.V5_S149", "VPL01.131_S150"))
colnames(tpm_selected)<-colnames(rt)
selected_genes<-rownames(rt)
tpm_data<-tpm_selected[selected_genes,]
pdf(file="TPM/heatmap_HighLow-Normal.pdf",width = 6,height = 8)
pheatmap(tpm_data,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         gaps_row = 14,gaps_col = 7,angle_col = "45", annotation_colors = ann_colors, 
         main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
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
write.csv(y,"GOBP_enrich_up_HighLow-Normal.csv",row.names =FALSE)
#结果可视化
GOBP_up<-summary(go_BP)
GOBP_up<-GOBP_up[1:20,]
GOBP_up$yvalue<-GOBP_up$Count/10
GOBP_up$xvalue<-GOBP_up$Description
library(ggpubr)
pdf("GOBP_up_HighLow-Normal.barplot.pdf",height=6,width = 12)
ggbarplot(GOBP_up, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("red","pink"))+
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
write.csv(y,"GOBP_enrich_down_HighLow-Normal.csv",row.names =FALSE)
#结果可视化
GOBP_down<-summary(go_BP)
GOBP_down<-GOBP_down[1:20,]
GOBP_down$yvalue<-GOBP_down$Count/6
GOBP_down$xvalue<-GOBP_down$Description
pdf("GOBP_down.barplot_HighLow-Normal.pdf",height=6,width = 12)
ggbarplot(GOBP_down, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("blue","purple"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()



######################compare5: highnormal -- low##########################################################
setwd("/Users/liuzhe/Desktop/UCSF/milk/analysis2/six_comparsion/HighNormal_Low/")
sample_selected<-c(metadata_high$Sample.name, metadata_normal$Sample.name, metadata_low$Sample.name)
sample_selected<-gsub("-", "\\.", sample_selected)
counts_anno_highnormal_low<-subset(gset, select = sample_selected)
group_list<-c("HighNormal", "HighNormal", "HighNormal",
              "HighNormal", "HighNormal", "HighNormal", "HighNormal","HighNormal", "HighNormal", "HighNormal",
              "Low", "Low", "Low", "Low")
group_list <- factor(group_list,levels = c("HighNormal", "Low"))
## 表达矩阵
data = apply(counts_anno_highnormal_low, 2, as.integer) ## DESeq2分析需要是整数
row.names(data) <- row.names(counts_anno_highnormal_low)
## 4.2分组矩阵，配对分析与常规分析最大的区别就在分组矩阵
condition = group_list
# 配对分析要加上这段代码，知道谁和谁是一对，比如1,1是一对，5,5是一对
coldata <- data.frame(row.names = colnames(data), condition)
# 注意在design中加上配对信息
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design = ~condition) 
dds$condition<- relevel(dds$condition, ref = "Low") 
## 4.3差异表达矩阵，还是和常规分析一样
dds <- DESeq(dds)
nrDEG_DESeq2 <- as.data.frame(results(dds))
rld <- rlog(dds)
# 这里我还提取了标准化后的表达矩阵，可以用于后续的热图绘制等等
normal_gset <- assay(rld) 
nrDEG_DESeq2 = nrDEG_DESeq2[order(nrDEG_DESeq2$log2FoldChange),] 
## 4.4定义差异基因
nrDEG <- nrDEG_DESeq2
nrDEG$Group = "notsignificant"
logFC_cutoff <- 1
nrDEG$Group[which( (nrDEG$padj < 0.05) & (nrDEG$log2FoldChange > logFC_cutoff) )] = "upregulated"
nrDEG$Group[which( (nrDEG$padj < 0.05) & (nrDEG$log2FoldChange < -logFC_cutoff) )] = "downregulated"
table(nrDEG$Group)
#downregulated notsignificant    upregulated 
#            9          21717              7 
write.csv(nrDEG,"DEGs_for_HighNormal-Low.csv",quote = F)
library("ggplot2")
deg<-nrDEG
head(deg)
deg$color<-ifelse(deg$padj<0.05 & abs(deg$log2FoldChange)>logFC_cutoff, ifelse(deg$log2FoldChange< -logFC_cutoff, "blue", "red"),"gray")
color<-c(red = "red", gray = "gray", blue = "blue")
deg$symbol<-row.names(deg)
p <- ggplot(data = deg, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=Group)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw() +
  xlim(-10,10)
for_label <- deg %>% 
  filter(abs(log2FoldChange) >1& padj < 0.05)
pdf("VP_DEGs_HighNormal-Low.pdf")
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
exp_deg<-subset(exp_deg, select = c("ID1024.V3_S139", "ID1025.V1_S138", "ID1026.V1_S137", 
                                    "ID1016.V1_S0", "ID1016.V2_S0", "ID1020.V1_S145", "ID1020.V2_S146",
                                    "ID1020.V3_S148", "ID1020.V5_S149", "VPL01.131_S150",
                                    "ID1018.V2_S0", "ID1023.V2_S140", "ID1028.V1_S0", "ID1029.V1_S0",
                                    "log2FoldChange", "pvalue","padj","Group"))
colnames(exp_deg)<-c("High_ID1024.V3_S139", "High_ID1025.V1_S138", "High_ID1026.V1_S137", 
                     "Normal_ID1016.V1_S0", "Normal_ID1016.V2_S0", "Normal_ID1020.V1_S145", "Normal_ID1020.V2_S146",
                     "Normal_ID1020.V3_S148", "Normal_ID1020.V5_S149", "Normal_VPL01.131_S150",
                     "Low_ID1018.V2_S0", "Low_ID1023.V2_S140", "Low_ID1028.V1_S0", "Low_ID1029.V1_S0", "log2FoldChange",
                     "raw_Pvalue","adjustedPvalue","Regulation" )
exp_deg<-cbind(round(exp_deg[1:17],3),exp_deg[,18])
colnames(exp_deg)[18]<-"regulation"
write.csv(exp_deg,"Table5_HighNormal-Low.csv",quote = F)
library("pheatmap")
rt<-exp_deg[,1:14]
metadata_highnormal_low<-rbind(metadata_high,metadata_normal, metadata_low)
annotation_col = data.frame(Type = factor(metadata_highnormal_low$Supply))
rownames(annotation_col) = colnames(exp_deg)[1:14]
annotation_row = matrix(exp_deg$regulation,nrow=length(exp_deg$regulation), ncol=1,byrow=TRUE)
rownames(annotation_row) = rownames(exp_deg)
colnames(annotation_row)<-"gene_type"
ann_colors = list(Type = c("Normal" = "blue", "Low" = "orange", "High" = "red"))
pdf(file="DESeq2/heatmap_HighNormal-Low.pdf",width = 6,height = 8)
pheatmap(rt,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         gaps_row = 9,gaps_col = 10,angle_col = "45", annotation_colors = ann_colors, 
         main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
tpm$geneid<-rownames(tpm)
tpm$geneid<-gsub("-","\\.",tpm$geneid)
rownames(tpm)<-tpm$geneid
tpm$geneid<-NULL
tpm_selected<-subset(tpm, select = c("ID1016.V1_S0", "ID1020.V3_S148", "ID1020.V2_S146",
                                     "ID1016.V1_S0", "ID1016.V2_S0", "ID1020.V1_S145", "ID1020.V2_S146",
                                     "ID1020.V3_S148", "ID1020.V5_S149", "VPL01.131_S150",
                                     "ID1018.V2_S0", "ID1023.V2_S140", "ID1028.V1_S0", "ID1029.V1_S0"))
colnames(tpm_selected)<-colnames(rt)
selected_genes<-rownames(rt)
tpm_data<-tpm_selected[selected_genes,]
pdf(file="TPM/heatmap_HighNormal-Low.pdf",width = 6,height = 8)
pheatmap(tpm_data,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         gaps_row = 9,gaps_col = 10,angle_col = "45", annotation_colors = ann_colors, 
         main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
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
write.csv(y,"GOBP_enrich_up_HighNormal-Low.csv",row.names =FALSE)
#结果可视化
GOBP_up<-summary(go_BP)
GOBP_up<-GOBP_up[1:20,]
GOBP_up$yvalue<-GOBP_up$Count/3
GOBP_up$xvalue<-GOBP_up$Description
library(ggpubr)
pdf("GOBP_up_HighNormal-Low.barplot.pdf",height=6,width = 12)
ggbarplot(GOBP_up, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("red","pink"))+
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
write.csv(y,"GOBP_enrich_down_HighNormal-Low.csv",row.names =FALSE)
#结果可视化
GOBP_down<-summary(go_BP)
GOBP_down<-GOBP_down[1:20,]
GOBP_down$yvalue<-GOBP_down$Count/6
GOBP_down$xvalue<-GOBP_down$Description
pdf("GOBP_down.barplot_HighNormal-Low.pdf",height=6,width = 12)
ggbarplot(GOBP_down, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("blue","purple"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()



######################compare6: lownormal -- high##########################################################
setwd("/Users/liuzhe/Desktop/UCSF/milk/analysis2/six_comparsion/LowNormal_High/")
sample_selected<-c(metadata_low$Sample.name, metadata_normal$Sample.name, metadata_high$Sample.name)
sample_selected<-gsub("-", "\\.", sample_selected)
counts_anno_lownormal_high<-subset(gset, select = sample_selected)
group_list<-c("LowNormal", "LowNormal", "LowNormal", "LowNormal",
              "LowNormal", "LowNormal", "LowNormal", "LowNormal","LowNormal", "LowNormal", "LowNormal",
              "High", "High", "High")
group_list <- factor(group_list,levels = c("LowNormal", "High"))
## 表达矩阵
data = apply(counts_anno_lownormal_high, 2, as.integer) ## DESeq2分析需要是整数
row.names(data) <- row.names(counts_anno_lownormal_high)
## 4.2分组矩阵，配对分析与常规分析最大的区别就在分组矩阵
condition = group_list
# 配对分析要加上这段代码，知道谁和谁是一对，比如1,1是一对，5,5是一对
coldata <- data.frame(row.names = colnames(data), condition)
# 注意在design中加上配对信息
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design = ~condition) 
dds$condition<- relevel(dds$condition, ref = "High") 
## 4.3差异表达矩阵，还是和常规分析一样
dds <- DESeq(dds)
nrDEG_DESeq2 <- as.data.frame(results(dds))
rld <- rlog(dds)
# 这里我还提取了标准化后的表达矩阵，可以用于后续的热图绘制等等
normal_gset <- assay(rld) 
nrDEG_DESeq2 = nrDEG_DESeq2[order(nrDEG_DESeq2$log2FoldChange),] 
## 4.4定义差异基因
nrDEG <- nrDEG_DESeq2
nrDEG$Group = "notsignificant"
logFC_cutoff <- 1
nrDEG$Group[which( (nrDEG$padj < 0.05) & (nrDEG$log2FoldChange > logFC_cutoff) )] = "upregulated"
nrDEG$Group[which( (nrDEG$padj < 0.05) & (nrDEG$log2FoldChange < -logFC_cutoff) )] = "downregulated"
table(nrDEG$Group)
#downregulated notsignificant    upregulated 
#         37          21644             52 
write.csv(nrDEG,"DEGs_for_LowNormal-High.csv",quote = F)
library("ggplot2")
deg<-nrDEG
head(deg)
deg$color<-ifelse(deg$padj<0.05 & abs(deg$log2FoldChange)>logFC_cutoff, ifelse(deg$log2FoldChange< -logFC_cutoff, "blue", "red"),"gray")
color<-c(red = "red", gray = "gray", blue = "blue")
deg$symbol<-row.names(deg)
p <- ggplot(data = deg, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=Group)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw() +
  xlim(-10,10)
for_label <- deg %>% 
  filter(abs(log2FoldChange) >1& padj < 0.05)
pdf("VP_DEGs_LowNormal-High.pdf")
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
exp_deg<-subset(exp_deg, select = c("ID1018.V2_S0", "ID1023.V2_S140", "ID1028.V1_S0", "ID1029.V1_S0",
                                    "ID1016.V1_S0", "ID1016.V2_S0", "ID1020.V1_S145", "ID1020.V2_S146",
                                    "ID1020.V3_S148", "ID1020.V5_S149", "VPL01.131_S150",
                                    "ID1024.V3_S139", "ID1025.V1_S138", "ID1026.V1_S137", 
                                    "log2FoldChange", "pvalue","padj","Group"))
colnames(exp_deg)<-c("Low_ID1018.V2_S0", "Low_ID1023.V2_S140", "Low_ID1028.V1_S0", "Low_ID1029.V1_S0",
                     "Normal_ID1016.V1_S0", "Normal_ID1016.V2_S0", "Normal_ID1020.V1_S145", "Normal_ID1020.V2_S146",
                     "Normal_ID1020.V3_S148", "Normal_ID1020.V5_S149", "Normal_VPL01.131_S150",
                     "High_ID1024.V3_S139", "High_ID1025.V1_S138", "High_ID1026.V1_S137", "log2FoldChange",
                     "raw_Pvalue","adjustedPvalue","Regulation" )
exp_deg<-cbind(round(exp_deg[1:17],3),exp_deg[,18])
colnames(exp_deg)[18]<-"regulation"
write.csv(exp_deg,"Table6_LowNormal-High.csv",quote = F)
library("pheatmap")
rt<-exp_deg[,1:14]
metadata_lownormal_high<-rbind(metadata_low,metadata_normal, metadata_high)
annotation_col = data.frame(Type = factor(metadata_lownormal_high$Supply))
rownames(annotation_col) = colnames(exp_deg)[1:14]
annotation_row = matrix(exp_deg$regulation,nrow=length(exp_deg$regulation), ncol=1,byrow=TRUE)
rownames(annotation_row) = rownames(exp_deg)
colnames(annotation_row)<-"gene_type"
ann_colors = list(Type = c("Normal" = "blue", "Low" = "orange", "High" = "red"))
pdf(file="DESeq2/heatmap_LowNormal-High.pdf",width = 6,height = 8)
pheatmap(rt,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         gaps_row = 37,gaps_col = 11,angle_col = "45", annotation_colors = ann_colors, 
         main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
tpm$geneid<-rownames(tpm)
tpm$geneid<-gsub("-","\\.",tpm$geneid)
rownames(tpm)<-tpm$geneid
tpm$geneid<-NULL
tpm_selected<-subset(tpm, select = c("ID1018.V2_S0", "ID1023.V2_S140", "ID1028.V1_S0", "ID1029.V1_S0",
                                     "ID1016.V1_S0", "ID1016.V2_S0", "ID1020.V1_S145", "ID1020.V2_S146",
                                     "ID1020.V3_S148", "ID1020.V5_S149", "VPL01.131_S150",
                                     "ID1016.V1_S0", "ID1020.V3_S148", "ID1020.V2_S146"))
colnames(tpm_selected)<-colnames(rt)
selected_genes<-rownames(rt)
tpm_data<-tpm_selected[selected_genes,]
pdf(file="TPM/heatmap_LowNormal-High.pdf",width = 6,height = 8)
pheatmap(tpm_data,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         gaps_row = 37,gaps_col = 11,angle_col = "45", annotation_colors = ann_colors, 
         main = "Title", fontsize = 6)
#fontsize = 1
dev.off()
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
write.csv(y,"GOBP_enrich_up_LowNormal-High.csv",row.names =FALSE)
#结果可视化
GOBP_up<-summary(go_BP)
GOBP_up<-GOBP_up[1:20,]
GOBP_up$yvalue<-GOBP_up$Count/43
GOBP_up$xvalue<-GOBP_up$Description
library(ggpubr)
pdf("GOBP_up_LowNormal-High.barplot.pdf",height=6,width = 12)
ggbarplot(GOBP_up, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("red","pink"))+
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
write.csv(y,"GOBP_enrich_down_LowNormal-High.csv",row.names =FALSE)
#结果可视化
GOBP_down<-summary(go_BP)
GOBP_down<-GOBP_down[1:20,]
GOBP_down$yvalue<-GOBP_down$Count/34
GOBP_down$xvalue<-GOBP_down$Description
pdf("GOBP_down.barplot_LowNormal-High.pdf",height=6,width = 12)
ggbarplot(GOBP_down, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("blue","purple"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()


