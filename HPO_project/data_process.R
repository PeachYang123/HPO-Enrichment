library(ape)
library(systemPipeR)
library(systemPipeRdata)
library(DESeq2)
getwd()
setwd("D:/自然语言处理/项目_work")
genWorkenvir(workflow = "rnaseq")
setwd("rnaseq")
gene_counts <- read.table(file.choose(), header=TRUE)
sampleNames <- names(gene_counts)
countData <- as.matrix(gene_counts)
rownames(countData) <- rownames(gene_counts)
database <- data.frame(name=sampleNames, condition=c("h", "s", "s", "h", "s", "h", "h", "s", "h", "s"))
rownames(database) <- sampleNames
dds <- DESeqDataSetFromMatrix(countData, colData=database, design= ~condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, "D:/自然语言处理/项目_work/res_des_output.xlsx")
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata, "D:/自然语言处理/项目_work/all_des_output.csv", row.names=FALSE)

plotMA(res, main="DESeq2", ylim=c(-2, 2))

library(ggplot2) 
resdata$change <- as.factor( ifelse( 
                                     resdata$padj<0.01 & abs(resdata$log2FoldChange)>1, 
                                     ifelse(resdata$log2FoldChange>1, "Up", "Down"), 
                                     "NoDiff" )
                             )
valcano <- ggplot(data=resdata, 
                  aes(x=log2FoldChange, y=-log10(padj), 
                      color=change)) +
          geom_point(alpha=0.8, size=1) + 
          theme_bw(base_size=15) +
          theme(panel.grid.minor=element_blank(), 
                panel.grid.major=element_blank() ) + 
          ggtitle("DESeq2 Valcano") + 
          scale_color_manual(name="", 
                             values=c("red", "green", "black"), 
                             limits=c("Up", "Down", "NoDiff")) + 
          geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
          geom_hline(yintercept=-log10(0.01), lty=2, col="gray", lwd=0.5) 
valcano

#绘制热图
library(pheatmap)
nt <- normTransform(dds)
exp <- log2(gene_counts+0.01)
df <- as.data.frame(colData(dds)[, c("name","condition")])
pheatmap(exp, cluster_rows=T, show_rownames=F, cluster_cols=T,annotation_col=df, fontsize=6)



#筛选差异表达基因
summary(res)
table(res$padj<0.05)
diff_gene_deseq2 <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene_deseq2)
head(diff_gene_deseq2)


#基因id转换
exp <- read.delim("D:/自然语言处理/项目_work/GSE174704_gene_counts.txt/gene_counts.txt",stringsAsFactors = FALSE)
data = data.frame(exp)
library(stringi)
data$Ensembl_ID=stri_sub(data$Ensembl_ID,1,15)#保留前15位
ensembl_id <- stri_sub(rownames(diff_gene_deseq2),1,15)
library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
gene_symbol <- bitr(ensembl_id, fromType="ENSEMBL", toType=c("SYMBOL", "ENTREZID"), OrgDb="org.Hs.eg.db")
head(gene_symbol)
#匹配到表达列表
data=data.frame(gene_symbol,diff_gene_deseq2[match(gene_symbol$ENSEMBL,rownames(diff_gene_deseq2)),])
df <- merge(gene_symbol,df)
df=df[,-4]#去除重复的Ensembl_ID列
write.csv(data, "D:/自然语言处理/项目_work/data.csv")
write.csv(gene_symbol, "D:/自然语言处理/项目_work/gene_symbol_data.csv")

