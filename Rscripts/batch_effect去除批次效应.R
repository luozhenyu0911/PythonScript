
#=============================================================================About CPM and TPM
#Counts
countsdata<-read.table(file = paste0(getwd(), '/expr-matrix.txt'),
                       header = T,sep = '\t',quote = '',check.names = F)  
exprSet<-countsdata[,7:ncol(countsdata)]
rownames(exprSet)<-countsdata[,1]
prefix<-"counts"

#CPM
cpm <- t(t(exprSet)/colSums(exprSet) * 1000000)
avg_cpm <- data.frame(avg_cpm=rowMeans(cpm))

#TPM 
kb <- countsdata$Length / 1000
rpk <- exprSet/ kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
avg_tpm <- data.frame(avg_tpm=rowMeans(tpm))

write.csv(avg_tpm,paste0(prefix,"_avg_tpm.csv"))
write.csv(avg_cpm,paste0(prefix,"_avg_cpm.csv"))
write.csv(tpm,paste0(prefix,"_tpm.csv"))
write.csv(cpm,paste0(prefix,"_cpm.csv"))


#======================================================================================================
#用boxplot()或者聚类查看是否需要去除批次。如果没有，就不去除。
rm(list=ls())
options(stringsAsFactors = F)

# library(bladderbatch)
# data(bladderdata)  
# class(bladderEset)    # bladder的属性是EsetExpressionSet，所以可用pData和exprs方法.
# pheno <- pData(bladderEset)
# edata <- exprs(bladderEset) 

FPKM<-read.table(file = paste0(getwd(), '/fpkmmerge.txt'),
                 header = T,sep = '\t',quote = '',check.names = F)  
rownames(FPKM)<-FPKM[,1]
exprSet <-FPKM[,-1]
exprSet<-as.matrix(exprSet)     # 表达矩阵 matrix
exprSet[is.na(exprSet)] <- 0    # NA值转为0
exprSet <-exprSet[rowSums(exprSet)>0,]       #过滤0值，必要，否则报错.
#exprSet <-exprSet[rowSums(exprSet==0)==0,]  #过滤含有0值的行
dim(exprSet)

batch1=exprSet[,c(1,2,5,6,11,12,15,16)]
batch2=exprSet[,c(9,10,19,20)]
batch3=exprSet[,c(3,4,7,8,13,14,17,18)]

par(mfrow = c(2, 1))
boxplot(log2(batch1+1),las=2,col=rainbow(7))
boxplot(log2(batch2+1),las=2,col=rainbow(7))
boxplot(log2(batch3+1),las=2,col=rainbow(7))
dev.off()
boxplot(log2(exprSet+1),las=2)

#boxplot(log2(exp1+1),las=2) #批次1
#boxplot(log2(exp2+1),las=2) #批次2
#若批次2 中有某一样本boxplot的中位数偏离其他样本，
#则使用normalizeBetweenArrays()拉此样本回到正常水平，该函数只能是在同一个数据集里面使用。
#library(limma)
#exp2_norm=normalizeBetweenArrays(log2(exp2+1))
#boxplot(exp_norm,las=2)

#x_merge=cbind(exp1,exp2)#两个批次数据合并
#boxplot(log2(x_merge+1))


coldata <-read.csv('phenodata-batch.csv',row.names = 1)

#--------------------------------------------------------------------------------------ComBat,效果最好。
#Sva包使用中在没有对数据集进行0值过滤或者有NA值时会错误，这时候需要去除0值和NA值。
# use FPKM  
library(sva)

#使用Hcluster查看聚类情况
dist_matrix <- dist(t(exprSet)) 
cluster <- hclust(dist_matrix, method = "complete") 
par(mfrow = c(2, 1))
plot(cluster, labels = coldata$batch) 
#plot(cluster, labels = coldata$Condition)
plot(cluster, labels = coldata$sample)
dev.off()

#校正批次效应,model可以有,也可以没有。如有，即告诉combat，有些分组本来就有差别，不要矫枉过正。
model = model.matrix(~as.factor(Condition), data=coldata)  #matrix
combat_edata <- ComBat(dat = exprSet, batch = coldata$batch, mod = model) #matrix
#combat_edata_nomod <- ComBat(dat = exprSet, batch = coldata$batch)
dim(combat_edata)
#dim(combat_edata_nomod)

#table(combat_edata_nomod<0)
#combat_edata_nomod[combat_edata_nomod<0]<- 0   # 负值转为0
#combat_edata_nomod <-combat_edata_nomod[rowSums(ccombat_edata_nomod)> 0,] 

table(combat_edata<0)
combat_edata[combat_edata<0]<- 0   # 负值转为0
combat_edata <-combat_edata[rowSums(combat_edata)> 0,] 
boxplot(log2(combat_edata+1),las=2)

write.csv(combat_edata,"combat_data_batch.csv") #作为差异分析新的表达矩阵。
#write.csv(combat_edata_nomod,"combat_data_nomod.csv")

dist_mat_combat <- dist(t(combat_edata)) 
clustering_combat <- hclust(dist_mat_combat, method = "complete") 
par(mfrow = c(3, 1))
plot(clustering_combat, labels = coldata$batch) 
plot(clustering_combat, labels = coldata$Condition)
plot(clustering_combat, labels = coldata$sample)
dev.off()

#pca
library(ggfortify)
library(ggplot2)

pheno<- read.csv('phenodata-batch.csv',row.names = 1)
screeplot(prcomp(t(exprSet)), type='lines')  #碎石图，决定选择几个主成分。
#df <- cbind(data,pheno)

before_noframe <- autoplot(prcomp(t(exprSet)), data=pheno, colour='Condition',size=5,
              shape='batch', frame=F, frame.type = 'norm') +
  geom_text(aes(label=rownames(t(exprSet)), vjust=1, colour = Condition),size=4)+
  #geom_vline(xintercept = 0, colour="#990000",linetype="dashed") +
  #geom_hline(yintercept = 0, colour="#990000",linetype="dashed") +
  theme_bw()
before_noframe
ggsave(before_noframe,filename = 'PCA-noSCALE-before-noframe.png')
dev.off()

before<- autoplot(prcomp(t(exprSet)), data=pheno, colour='Condition',size=5,
                           shape='batch', frame=T, frame.type = 'norm') +
  #geom_text(aes(label=rownames(t(exprSet)), vjust=1, colour = Source),size=4)+
  #geom_vline(xintercept = 0, colour="#990000",linetype="dashed") +
  #geom_hline(yintercept = 0, colour="#990000",linetype="dashed") +
  theme_bw()
before
ggsave(before,filename = 'PCA-noSCALE-before.png')
dev.off()

after <- autoplot(prcomp(t(combat_edata)), data=pheno, colour='Condition',size=5,
                  shape='Condition', frame=T, frame.type = 'norm') +
  geom_text(aes(label=rownames(t(combat_edata)), vjust=1, colour = Condition),size=4)+
  #geom_vline(xintercept = 0, colour="#990000",linetype="dashed") +
  #geom_hline(yintercept = 0, colour="#990000",linetype="dashed") +
  theme_bw()
after
ggsave(after,filename = 'PCA-noSCALE-after.png')
dev.off()

after_noframe <- autoplot(prcomp(t(combat_edata)), data=pheno, colour='Condition',size=5,
                  shape='Condition', frame=F, frame.type = 'norm') +
  geom_text(aes(label=rownames(t(combat_edata)), vjust=1, colour = Condition),size=4)+
  #geom_vline(xintercept = 0, colour="#990000",linetype="dashed") +
  #geom_hline(yintercept = 0, colour="#990000",linetype="dashed") +
  theme_bw()
after_noframe
ggsave(after_noframe,filename = 'PCA-noSCALE-after-noframe.png')
dev.off()

save(before,before_noframe,after,after_noframe,file = 'pca.Rdata')
load('pca.Rdata')

#----------------------------------------------------------------------limma::removeBatchEffect
#1.理解boxplot的重要性，来看数据集是否需要log，以便后面才能用limma包进行差异分析.
#2.normalizeBetweenArrays只能是在同一个数据集里面用来去除样本的差异，
#不同数据集需要用limma 的 removeBatchEffect函数去除批次效应。
#关于log，因为做差异分析的limma包要求表达矩阵exprSet中的数据是经过log的。
rm(list = ls())
options(stringsAsFactors = F)
library(limma)

coldata <-read.csv('phenodata-batch.csv',row.names = 1)
coldata$Condition <-as.factor(coldata$Condition)
batch <-as.factor(coldata$batch)
group_list <- coldata$Condition

design <- model.matrix(~0+group_list)
colnames(design)=levels(group_list)
rownames(design)=rownames(coldata)
design

cont.matrix=makeContrasts(contrasts=c('LN-CK'),levels = design)  # contrasts:trt-untrt
cont.matrix

fpkm<-read.table(file = paste0(getwd(), '/fpkmmerge.txt'),
                       header = T,sep = '\t',quote = '',check.names = F,row.names = 1)  
head(fpkm)
fpkm <-fpkm[rowSums(fpkm)>0,]  
exprSet<-fpkm
exprSet <-exprSet[rowSums(exprSet)>0,]  

expr_limma <- removeBatchEffect(exprSet, batch = coldata$batch, design = design)

boxplot(log2(exprSet+1),las=2)
boxplot(log2(expr_limma+1),las=2)

#差异分析
fit=lmFit(expr_limma,design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)

options(digits = 4)
allGeneSets <- topTable(fit2, coef='LN-CK', number=Inf) 
DEgeneSets <- topTable(fit2, coef='LN-CK', number=Inf, p.value=0.05, adjust="BH")
res <- decideTests(fit2, p.value=0.05)
summary(res)
df=allGeneSets[allGeneSets$adj.P.Val<0.05 ,]
df

dist_mat <- dist(t(expr_limma)) 
clustering <- hclust(dist_mat, method = "complete") 
plot(clustering, labels = coldata$batch) 
plot(clustering, labels = coldata$Condition)

#pca

library(ggfortify)
library(ggplot2) 

pheno<- read.csv('phenodata-source.csv',row.names = 1)

data <- plotPCA(vst, "batch", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
p<- ggplot(data, aes(PC1, PC2, color=pheno$Source,shape=pheno$Condition)) +
  geom_point(size=4) +
  #geom_text(aes(label=rownames(data),vjust=1))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()
p
ggsave('PCA-limma-after-noframe.png')

after <- autoplot(prcomp(t(vstMat)), data=pheno, colour='Source',size=4,
                  shape='Condition') +
  #geom_text(aes(label=rownames(t(combat_edata)), vjust=1, colour = Source,size=4))+
  #geom_vline(xintercept = 0, colour="#990000",linetype="dashed") +
  #geom_hline(yintercept = 0, colour="#990000",linetype="dashed") +
  theme_bw()
dev.off()
after
ggsave('PCA-gglimma-after-noframe.png')

#-----------------------------------------------------------------------------preprocessCore  
rm(list=ls())
options(stringsAsFactors = F)

library(preprocessCore)
countsdata<-read.table(file = paste0(getwd(), '/fpkmmerge.txt'),
                       header = T,sep = '\t',quote = '',check.names = F)
exprSetm<-countsdata[,-1]
rownames(exprSetm)<-countsdata[,1]
exprSet<-as.matrix(exprSetm)
exprSet[is.na(exprSet)] <- 0    # NA值转为0
exprSet <-exprSet[rowSums(exprSet)>0,]  
dataMat <- exprSet

coldata <-read.csv('phenodata-batch.csv',row.names = 1)

dataMatNorm <- normalize.quantiles(dataMat) #将数据按分位数法标准化去除批次效应
dataMatNorm <- as.data.frame(dataMatNorm)
names(dataMatNorm) <-names(exprSetm)
rownames(dataMatNorm) <-rownames(dataMat)
dim(dataMatNorm)
table(dataMatNorm<0)
dataMatNorm[dataMatNorm<0]<- 0   # 负值转为0
dataMatNorm <-dataMatNorm[rowSums(dataMatNorm)> 0,] 

write.csv(dataMatNorm,"dataMatNorm.csv")

whichbatch <- coldata$batch
data1_removeBE_quantile=dataMatNorm[, whichbatch=="1"]
data2_removeBE_quantile=dataMatNorm[, whichbatch=="2"]
data3_removeBE_quantile=dataMatNorm[, whichbatch=="3"]
data4_removeBE_quantile=dataMatNorm[, whichbatch=="4"]

dist_mat <- dist(t(dataMatNorm)) 
clustering <- hclust(dist_mat, method = "complete") 
plot(clustering, labels = coldata$batch) 
plot(clustering, labels = coldata$Condition)
plot(clustering, labels = coldata$sample.1)
#pca
library(ggfortify)
library(ggplot2)

pheno<- read.csv('phenodata-source.csv',row.names = 1)

before <- autoplot(prcomp(t(exprSet)), data=pheno, colour='Source',size=4,
                   shape='Condition', frame=T, frame.type = 'norm') +
  #geom_text(aes(label=rownames(t(exprSet)), vjust=1, colour = Source),size=4)+
  #geom_vline(xintercept = 0, colour="#990000",linetype="dashed") +
  #geom_hline(yintercept = 0, colour="#990000",linetype="dashed") +
  theme_bw()
dev.off()
before
#ggsave('PCA-noSCALE-before-noframe.png')

after <- autoplot(prcomp(t(dataMatNorm)), data=pheno, colour='Condition',size=4,
                  shape='Condition') +
  geom_text(aes(label=rownames(t(combat_edata)), vjust=1, colour = Condition),size=4)+
  #geom_vline(xintercept = 0, colour="#990000",linetype="dashed") +
  #geom_hline(yintercept = 0, colour="#990000",linetype="dashed") +
  theme_bw()
dev.off()
after
ggsave('PCA-preprocessCore-after-noframe.png')

