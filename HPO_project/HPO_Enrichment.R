#载入相关包
library(fdrtool)
library(ggplot2)

#数据准备
setwd("E:/大三下/自然语言处理_资料/HPO_project")
#导入处理过的HPO数据，保留entrez-gene-id，HPO-Term-ID，HPO-Term-Name
HPO_Data_Set <- read.table("./genes_to_phenotype_edit.txt", header = T, sep = "\t", quote = "")
#导入去重后的HPO-Term-ID，HPO-Term-Name
HPOID_to_TermName <- read.table("./HPOID_to_TermName.txt", header = T, sep = "\t", quote = "")
#导入差异基因集合，已转化为enrezid
Gene_Set <- unique(read.table('df_gene.csv',header=T))

#计算每个entrezid对应的HPOID,生成list
HPO_GENE <- data.frame(HPO_Data_Set$entrez.gene.id)
HPO_GENE_UN <- unique(HPO_GENE)
GENE_to_HPOID <- sapply(HPO_GENE_UN$HPO_Data_Set.entrez.gene.id, function(x) HPO_Data_Set$HPO.Term.ID[HPO_Data_Set$entrez.gene.id==x])
for (i in 1:length(GENE_to_HPOID)){
  names(GENE_to_HPOID)[i] <- HPO_GENE_UN$HPO_Data_Set.entrez.gene.id[i]
}
Gene_Set <- intersect(Gene_Set$ENTREZID, names(GENE_to_HPOID))
Gene_Set <- data.frame(ENTREZID=Gene_Set)

#计算每个HPOID对应的entrezid,生成list
HPOID <- data.frame(HPO_Data_Set$HPO.Term.ID)
HPOID_UN <- unique(HPOID)
HPOID_to_GENE <- sapply(HPOID_UN$HPO_Data_Set.HPO.Term.ID, function(x) HPO_Data_Set$entrez.gene.id[HPO_Data_Set$HPO.Term.ID==x])
for (i in 1:length(HPOID_to_GENE)){
  names(HPOID_to_GENE)[i] <- HPOID_UN$HPO_Data_Set.HPO.Term.ID[i]
}

#p值计算
#参数解释，在某个HPO term中：N:HPO中的所有gene数目；K:导入的所有差异基因数目；M:该term在HPO中对应的基因数目；Q：该term在差异基因中对应的基因数目
GetPValue <- function(Q,M,N,K)
  1 - sum(sapply(0:(Q - 1), function(i) choose(M, i) * choose(N - M, K - i)/ choose(N, K)))

#计算N,K值
N <- length(HPO_GENE_UN$HPO_Data_Set.entrez.gene.id)
K <- length(Gene_Set$ENTREZID)

#计算M
vM <- vector()
for (i in 1:length(HPOID_UN$HPO_Data_Set.HPO.Term.ID)){
  vM[i] <- length(HPOID_to_GENE[[i]])
}

#计算Q
vQ <- vector()
#计算每个差异基因对应的所有HPOID
Gene_Set_HPOID <- list()
num <- 1
for (i in 1:length(GENE_to_HPOID)){
  for (j in 1:length(Gene_Set$ENTREZID)){
    if(names(GENE_to_HPOID)[i]==as.character(Gene_Set$ENTREZID[j])) {
      Gene_Set_HPOID[num] <- GENE_to_HPOID[i]
      names(Gene_Set_HPOID)[num] <- Gene_Set$ENTREZID[j]
      num <- num+1
      }
  }
}
num2 <- 0
Gene_Set_HPOID_Unlist <- unlist(Gene_Set_HPOID)
for (i in 1:length(HPOID_UN$HPO_Data_Set.HPO.Term.ID)){
  for(j in 1:length(Gene_Set_HPOID_Unlist)){
    if(HPOID_UN$HPO_Data_Set.HPO.Term.ID[i]==Gene_Set_HPOID_Unlist[j]){
      num2 <- num2+1
    }
  }
  vQ[i] <- num2
  num2 <- 0
}

#p值计算和FDR矫正
pvalues <- sapply(1:length(HPOID_UN$HPO_Data_Set.HPO.Term.ID), function(i) GetPValue(vQ[i], vM[i], N, K))
fdr <- fdrtool(pvalues,statistic = "pvalue")
FDR <- fdr$qval

res <- data.frame(
  "HPO_ID" = HPOID_UN$HPO_Data_Set.HPO.Term.ID,
  "Term_Name" = HPOID_to_TermName$HPO.Term.Name,
  "GeneRatio" = paste(vQ, K, sep = '/'),
  "pvalues" = pvalues,
  "qvalues" = FDR,
  "Count" = vQ
)
res <- subset(res, Count>0)
res <- subset(res, pvalues<0.01)
res1 <- res[order(res$pvalues),]
res2 <- res[order(res$qvalues),]
write.csv(res1, "./result1.csv")
write.csv(res2, "./result2.csv")

#绘制图
# 条形图绘制
data <- read.csv("result2.csv",header = TRUE)
ggplot(data = data,aes(x =Term_Name ,y = Count,fill = -log10(qvalues))) + geom_bar(stat="identity") + 
  scale_x_discrete(limits=data$Term_Name) + coord_flip() + labs(title = "HPO_Enrichment") + 
  theme(plot.title = element_text(size = 20,face = "bold"),axis.text = element_text(size = 12,face = "bold"),
        axis.title.x =element_text(size=14), axis.title.y=element_text(size=16),panel.background 
        = element_rect(fill="white", colour='gray')) 
  