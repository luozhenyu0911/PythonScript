library(ggpubr)
library(ggthemes)
library(Rtsne)
setwd('D:\\ARyuyan\\strand-specificnew\\antisenseDEG')
data<-couts_FPKM[,2:ncol(couts_FPKM)]
colnames(couts_FPKM)
rownames(data)=couts_FPKM$X1
head(data)
#计算tsne
set.seed(1)
tsne.info<-Rtsne(t(data),perplexity=3)
#显示tsne计算结果前6行
head(tsne.info$Y)
colnames(tsne.info$Y)<-c("tSNE_1","tSNE_2")

tSNE.data<-data.frame(sample=colnames(data),
                      Type=c(rep("d4",3),rep("l4",3),
                             rep("d12",3),rep("l12",3)),
                      tsne.info$Y)
#绘制tsne散点图
ggscatter(tSNE.data, x="tSNE_1",y="tSNE_2",color = "Type",
          ellipse = T ,size = 4,
          shape = "Type",
          # palette = c("","","",""),设置颜色
          main="tSNE plot",
          label = "sample")+ theme_base ()
#------PCA----------
#### RNA-seq样本PCA分析####
#加载的R包
# install.packages(c("ggpubr","ggthemes","gmodels"))
library(ggpubr)
#加载差异基因表达矩阵
library(gmodels)
library(ggpubr)
library(ggplot2)
library(ggthemes)
head(data)#每一列为一个样本，每一行为一个基因
#计算PCA
pca.info <- prcomp(data,scale=TRUE)
#显示PCA计算结果
head(pca.info$rotation)
pca.data <- data.frame(sample = rownames(pca.info$rotation),
                       Type = c(rep("d4",3),rep("l4",3),
                                rep("d12",3),rep("l12",3)),
                       pca.info$rotation)
#绘图
ggscatter(pca.data,x = "PC1",y = "PC2",
          label = "sample",color = "Type") + theme_base()
###其它图形修饰参数自己摸索吧，哈哈~
