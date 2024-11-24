#Mfuzz为模糊聚类算法的一种
#模糊聚类又叫做软聚类，相比“硬聚类”，其最大的特点是允许同一数据属于多个不同的类
#即一个基因可能被聚到多类里面。
#用0-1的数字表示与所有分类的隶属度，即为membership。

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", update=F)
if (!requireNamespace("Mfuzz", quietly = TRUE))
  BiocManager::install("Mfuzz", update=F)


rm(list=ls())

library(Mfuzz)
library(RColorBrewer)

setwd('C:/Users/ZQQ/Desktop/植物生信SHOP网课资料/cluster/')

#加载FPKM数据，不要Z-score或log。
data <- read.delim(header=T,file='data/mfuzz.exp.txt',row.names = 1,
                   sep="\t",na.strings = "-",check.names=F)

n_clust=9  #人为指定分类数

#将矩阵转换为ExpressionSet类型，用于Mfuzz分类
eset <- new("ExpressionSet",exprs = as.matrix(data))
class(eset)

eset <- filter.std(eset,visu = F,min.std=0.3) # min.std，可以自行调节，根据标准差去除样本间差异太小的基因
eset <- standardise(eset) #标准化，使每个基因的均值为0，标准差为1
cl <- mfuzz(eset, c = n_clust, m = mestimate(eset) )  #m为模糊指数

#输出
write.table(cl$membership,file='mfuzz.membership.txt',
            col.names = NA,quote=F,sep="\t")#membership越大表明该基因属于某一类可能性越大。
write.table(cl$cluster,file='mfuzz.cluster.txt',
            col.names = F,quote=F,sep="\t") 


color=colorRampPalette(rev(brewer.pal(n=11,'Spectral')))(100)
color = rev(rainbow(10))
# 画图，每个cluster单独画,每图一页PDF
pdf('mfuzz.sep.pdf',h=8,w=10)
par(mar=c(3,5,3,2))
mfuzz.plot2(eset,cl, xlab='',ylab='Expression', #x、y轴名称
  centre=T,    #绘制中心线
  time.labels=colnames(data), x11=F,centre.col="blue",centre.lwd=2,
  #colo='fancy', #默认颜色
  colo=color,#也可以自定义颜色，提示警告可以忽略
  col.axis = 'black',col.lab = 'black',col.main = 'black',
  col.sub =  'black',col =  'black',
  cex=1.5,cex.lab=1.8,cex.axis=1.5,cex.main=2,cex.sub=1.8) 
dev.off()


# 画图，所有cluster画到一张图上
pdf('mfuzz.onepage.pdf',h=8,w=14)
par(mar=c(3,5,3,2))
mfuzz.plot2(eset,cl,
  xlab='',ylab='Expression',
  mfrow=c(3,3), #同一页显示，根据聚类数目分
  centre.col="black",centre.lwd=2,
  centre=T, time.labels=colnames(data),
  colo='fancy', # colo=color, 
  x11=F,cex=1.5,cex.lab=1.8,cex.axis=1.5,cex.main=2,cex.sub=1.8) 
dev.off()
