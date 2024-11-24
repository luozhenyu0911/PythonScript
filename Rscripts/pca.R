#利用ggarrange将多个图像排列在一起：
library(ggpubr)
ggarrange(p1,p2,p3,p4,p5,p6, ncol = 2,nrow =3,widths = c(1,2),heights = c(1,1,2))
#六个图，3行2列排列。

#===================================================================================================
#主成分分析是一种降维的数据处理方法。从生物学上说，
#目的是判断组内样本重复性是否够好（图上距离较近），组间样本的差异是否足够大（图上距离较远）。
#若PCA分析结果不好，则后续差异分析结果不可靠；
#若存在离群样本，可剔除该样本再进行后续分析，以确保后续结果有意义。


#画PCA图时要求：行是样本，列是gene。
#library("edgeR")  

#PCA注意事项:
#1.一般说来，PCA之前原始数据需要中心化（center=T）。
#除了中心化以外，定标(Scale=T) 也是数据处理前需要考虑的一点。
#如果数据没有定标，则原始数据中方差大的变量对主成分的贡献会很大。
#2.但定标(scale)可能会有一些负面效果:
#定标后，如果变量中有噪音,就在无形中把噪音和信息的权重变得相近，但PCA本身无法区分信号和噪音。
#在这样的情形下，做定标反而影响了我们对数据异常情况的检查和核对。
#如果变量之间的数据的处于不同数量级或者变量之间的均值/方差相差很大时，建议进行标准化scale.

rm(list = ls())
options(stringsAsFactors = F)

library(export)

fpkm <- read.table('strand_specific/expr/sense_fpkm.txt',row.names = 1,header = T,
                  sep = '\t')  # use counts :  data=log(edgeR::cpm(data)+1) 
dim(fpkm)

#必要步骤
fpkm=na.omit(fpkm)                  #过滤掉包含缺失的基因
fpkm=fpkm[rowSums(fpkm)!=0,]   #预处理，过滤掉所有样品中均为0的基因

dim(fpkm)
data = t(fpkm[order(apply(fpkm,1,mad), decreasing = T),] ) #[1:5000] # all,30000,15000,10000,8000,5000,2000
dim(data)
class(data)
data=as.data.frame(data)

pheno<-data.frame(Condition = factor(c(rep('CK',12),rep('LN',12)),levels = c('CK','LN')),
                  time=factor(c(rep('0.5h',2),rep('1h',2),rep('1d',2),rep('2d',2),rep('4d',2),rep('8d',2),
                                rep('0.5h',2),rep('1h',2),rep('1d',2),rep('2d',2),rep('4d',2),rep('8d',2))),
                  stage=factor(c(rep('CK_early_stage',4),rep('CK_mid_stage',4),rep('CK_late_stage',4),
                                 rep('LN_early_stage',4),rep('LN_mid_stage',4),rep('LN_late_stage',4))),
                  row.names = colnames(fpkm))

#=======================================================================================================
prcomp(data,scale.=T) #一般需要
library(ggfortify)
library(ggplot2)
library(ggrepel)

screeplot(prcomp(data), type='lines')   #碎石图，决定选择几个主成分。
df <- cbind(data,pheno)

p <- autoplot(prcomp(data,scale. = T), data=pheno, colour='time',size=6, #ellipse.border.remove=T,ellipse.type="convex",
              shape='Condition', frame=T, frame.type = 'norm') +
  #geom_label_repel(aes(label=rownames(data), fill=factor(time),vjust=1),
   #                size=4,show.legend = F)+
  #geom_vline(xintercept = 0, colour="#990000",linetype="dashed") +
  #geom_hline(yintercept = 0, colour="#990000",linetype="dashed") +
  theme_classic(base_size = 18)+theme(#legend.position = c(0.11,0.76),
                                 legend.background = element_rect(colour = 'black') )+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        #axis.line = element_line(colour = "black",size=1)
        ) 
p 
ggsave('PCA3.png')
graph2ppt(p,file='PCA-samples-frame-plot1.pptx')
graph2ppt(p,file='PCA-samples-frame-scale-plot1.pptx')
graph2ppt(p,file='PCA-samples-frame-scale-repel-plot1.pptx')

#================================================================================================
library(FactoMineR)
library(factoextra)  
library(ggplot2)

dat.pca <- PCA(data,graph = FALSE ,scale.unit = T)   #scale.unit = T is default
dat.pca
fviz_eig(dat.pca, addlabels = TRUE, ylim = c(0, 50))
p<-fviz_pca_ind(dat.pca,axes = c(1, 2), #选择要显示的主成分
                geom.ind = c("point", "text"),# show points and text
                col.ind = pheno$time,  # color by groups
                # palette = "lancet", #颜色方案，可选"npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty"
                palette = c('red',"green",'blue','navy','purple','yellow'),
                addEllipses = F,   # Concentration ellipses
                #ellipse.type="confidence", ellipse.level=0.95, #置信区间椭圆
                legend.title = "Groups",repel = T, alpha.ind=0.5, #透明度
                invisible="none", #如需去除重心点，设置为quali
                pointsize = 2, title="PC1 vs PC2", #主标题
                legend.title = "Group",#legend标题
                )+theme_minimal()+
     theme(plot.title = element_text(hjust = 0.5)) #设置主题与标题居中
p
ggsave('PCA.png')



#==================================================================================ggplot2
library(tidyr)
library(dplyr)

pca_data <- prcomp(data, scale. = T) # 主成分计算  scale = F or T is not the same.    
screeplot(pca_data, type = "lines")  #查看合适主成分个数
summary(pca_data)
rownames(pca_data$x)

x_per <- round(summary(pca_data)$importance[2, 1]*100, 1)#提取PC1的百分比
y_per <- round(summary(pca_data)$importance[2, 2]*100, 1)#提取PC2的百分比

df_sample <- cbind(pca_data$x,pheno)

#设置适合科研的背景色
theme_set(theme_bw(base_size = 18))
#绘图
plot_1 <- ggplot(df_sample, 
                 aes(x = PC1, y=PC2,label=rownames(pheno),
                     color =time,shape= Condition)) +
  geom_point(size=6)+
  #geom_text(size=6)+ #stat_ellipse(level = 0.95, size = 1) +
  xlab(paste("PC1","(", x_per,"%)",sep=" ")) +
  ylab(paste("PC2","(", y_per,"%)",sep=" ")) +
  geom_vline(xintercept = 0,lty =2) +geom_hline(yintercept = 0,lty=2) 
  #t图例按照G1，G2，G3...排列
  #scale_color_discrete(limits = c(paste("G", seq(1:2),sep = "")))
plot_1 
ggsave('pca.png')
graph2ppt(plot_1,file='PCA-samples-no-scale-plot1.pptx')
graph2ppt(plot_1,file='PCA-samples-scale-plot1.pptx')


library(ggpubr)
pca <- ggscatter(df_sample,x='PC1',y='PC2',
                 color ='time',shape= 'Condition', 
                 size = 5 #,main = "PCA plot"
                 )+
  xlab(paste("PC1","(", x_per,"%)",sep=" ")) +
  ylab(paste("PC2","(", y_per,"%)",sep=" ")) +
  theme_classic(base_size = 18)
pca
graph2ppt(pca,file='PCA-samples-plot.pptx')

#---------------------------------------------------------------------------------3D pca
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(grDevices))
library(gmodels)
suppressPackageStartupMessages(library(pca3d)) #3D pca绘制
suppressPackageStartupMessages(library(maptools)) #避免标签重复


#调用代码：
data <- read.table('strand_specific/expr/sense_fpkm.txt',row.names = 1,header = T,
                   sep = '\t')  # use counts :  data=log(edgeR::cpm(data)+1) 
dim(data)
head(data)
groups<-data.frame(Condition = factor(c(rep('CK',12),rep('LN',12)),levels = c('CK','LN')),
                  name= names(data),row.names = names(data)) 

pca <- fast.prcomp(data,scale. = T)   #gmodels包中的fast.prcomp函数计算PCA ,比内置函数更快。
head(pca$rotation)

pca.data <-data.frame(cbind(pca$rotation,groups))

color.lib <-c('red','blue')
color <-color.lib[as.numeric(pca.data$Condition)]
shape.lib <-c(16,17)
shape <-shape.lib[as.numeric(pca.data$Condition)]

s3d <- scatterplot3d(pca.data[,c('PC1','PC2','PC3')], color = color,pch= shape,  #形状
              angle = 60, 
              #main='PCA 3D plot' ,
              cex.symbols = 3 , #cex.symbols 点的大小
              ) 
legend('bottom',legend = levels(pca.data$Condition),xpd = T,horiz = T,
       col = color.lib,pch = shape.lib,inset = -0.18)
text(s3d$xyz.convert(pca.data[,c('PC1','PC2','PC3')] +0.02),
     labels = pca.data$name,cex=1,col = 'black')

#显示部分标签
lab.name <-c('CK.48h.R2','LN.24h.R1','CK.60min.R2')
pca.data$label <-''
pca.data$label[match(lab.name,pca.data$name)] <-lab.name
s3d <- scatterplot3d(pca.data[,c('PC1','PC2','PC3')], color = color,pch= shape,  #形状
                     angle = 60, main='PCA 3D plot' ,cex.symbols = 3 , #cex.symbols 点的大小
                     ) 
legend('bottom',legend = levels(pca.data$Condition),xpd = T,horiz = T,
       col = color.lib,pch = shape.lib,inset = -0.18)
text(s3d$xyz.convert(pca.data[,c('PC1','PC2','PC3')] +0.02),
     labels = pca.data$label , cex=1,col = 'black')

dev.off()


