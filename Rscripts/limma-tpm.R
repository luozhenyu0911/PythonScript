#limma的核心函数是lmFit和eBayes， 前者是用于线性拟合，后者根据前者的拟合结果进行统计推断。
#lmFit RNAseq需要:表达矩阵，分组对象, 差异比较矩阵
#与芯片数据不同，limma使用voom方法来对RNA-seq数据raw count进行normalization.
#limma包接受多种数据类型：counts（需要用voom进行normalization）,log2(rpkm/FPKM+1)。
#edgeR 和 DESeq2 主要接受raw counts矩阵数据。

#limma为高斯分布(正态分布)，可接受标准化后数据，如TPM.最好把FPKM转为TPM进行分析。  
#edgeR,DESeq2为负二项分布，接受read counts。

rm(list = ls()) 
options(stringsAsFactors = F)
library(ggplot2)
library(limma)

if(T){
countdata<-read.csv('GSE115796_All_gene_FPKM.csv',row.names = 1)
countdata <-na.omit(countdata)

fpkmToTpm <- function(fpkm) {exp(log(fpkm) - log(sum(fpkm)) + log(1e6))}
tpms <- apply(countdata,2,fpkmToTpm)   #使用TPM进行后续分析。

colSums(tpms)

exprSet <-log2(tpms+1)  #对tpm取log

}

#exprSet <- exprSet[,c(1,2,7,8)]
#exprSet <- exprSet[,c(1,2,13,14)]

#exprSet <- exprSet[,c(3,4,9,10)]
#exprSet <- exprSet[,c(3,4,15,16)]

#exprSet <- exprSet[,c(5,6,11,12)]
#exprSet <- exprSet[,c(5,6,17,18)]

exprSet <- exprSet[,c(15,16,17,18)]
#exprSet <- exprSet[,c(9,10,11,12)]
#exprSet <- exprSet[,c(3,4,5,6)]

head(exprSet)
exprSet<-exprSet[rowSums(exprSet)>0,]
boxplot(log2(exprSet+1),outline=T,las=2) 

if(T){
exprSet <-na.omit(exprSet)

#第一步：构建分组矩阵
sample=colnames(exprSet)
group=factor(c(rep('H2',ncol(exprSet)/2),rep('H3',ncol(exprSet)/2)),levels=c('H3','H2'))#control在前
phno=data.frame(sample,group)

design <- model.matrix(~0 + factor(phno$group))
colnames(design)=levels(factor(phno$group))
rownames(design)=colnames(exprSet)
design

cont.matrix=makeContrasts(contrasts=c('H2-H3'),levels = design)  # contrasts:trt-untrt
cont.matrix

#第三步：做差异分析，提取差异分析结果
fit <- lmFit(exprSet, design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)

tempOutput = topTable(fit2, coef='H2-H3', n=Inf)
DEG_limma_voom = na.omit(tempOutput)
#write.csv(DEG_limma_voom,"DEG_limma_voom.csv",quote = F) #自行筛选DEG
DEgenes <- topTable(fit2, coef="H2-H3", number=Inf,p.value=0.05, adjust="BH", lfc=1)
res <- decideTests(fit2, p.value=0.05, lfc=1)
summary(res)

resdata<-merge(as.data.frame(DEG_limma_voom), as.data.frame(exprSet), by='row.names',sort=F)                     
head(resdata)
resdata$significant <- 'unchanged'  
resdata$significant[resdata$adj.P.Val <= 0.05 & resdata$logFC>=1] <- 'upregulated'     
resdata$significant[resdata$adj.P.Val <= 0.05 & resdata$logFC<=-1] <- 'downregulated'  
head(resdata)
table(resdata$significant)

DEG=as.data.frame(DEG_limma_voom)
DEG=na.omit(DEG)  
nrDEG=DEG   
nrDEG$change = as.factor(ifelse(nrDEG$adj.P.Val < 0.05 & abs(nrDEG$logFC) >1, 
                                ifelse(nrDEG$logFC > 1 ,'Up','Down'),'NoDiff')) 
this_tile <- paste0("H2 vs H3", 
                    '\nup gene : ',nrow(nrDEG[nrDEG$change =='Up',]) , 
                    '\ndown gene : ',nrow(nrDEG[nrDEG$change =='Down',])) 
#this_tile <- 'blue'
valcano <- ggplot(data=nrDEG, aes(x=logFC, y=-log10(adj.P.Val), color=change)) + 
  geom_point(alpha=0.4, size=4) + 
  xlim(-8,8)+ylim(0,2)+
  theme_set(theme_set(theme_bw(base_size=20)))+ 
  xlab("log2FoldChange") + ylab("-log10(padj)") + 
  ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+  
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
  geom_hline(yintercept=-log10(0.05), lty=2, col="gray", lwd=0.5)+
  scale_colour_manual(values = c('red','blue','green'), limits=c("Up", "Down", "NoDiff")) 
valcano
}

summary(res)
table(resdata$significant)
ggsave(valcano,filename = 'volcano-limma-L2vsL3.png') 
write.csv(resdata,'resdata-limma-L2vsL3.csv')
