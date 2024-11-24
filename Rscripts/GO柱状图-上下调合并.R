rm(list=ls())
#上下调在一组展示。
library(ggplot2)
library(reshape2)
library(export)
library(XLConnect)
library(tidyr)

if(T){
Data0=readWorksheetFromFile('ssRNAseq/GO_DEGs.xlsx',2)[,c(1,3,4,10)]
Data<-spread(Data0,key='Group',value = 'GeneRatio')
Data[is.na(Data)] <- 0
colnames(Data)=c('Class','Term','Down Regulated','Up Regulated')

#将宽矩阵转换为长矩阵
dat=melt(Data,
         id.vars = c('Term','Class'),
         value.name = 'ratio',
         variable.name = 'direction')

#将ratio值转换为百分比
dat$ratio=dat$ratio*100

#设置不同分类的柱子颜色
colors=c('#00B454','#FF3900')

#设置y轴坐标点分隔（breaks）
ybreaks=pretty(dat$ratio) #自动设置
#设置主标题，x轴label，y轴label，以及legend标题
MainTitle=""
xLabel=""
yLabel="Percentage (%)"
legendTitle=""

#设置ggplot对象，映射Term为x轴，ratio为y轴
g<-ggplot(data=dat,aes(Term,ratio))

#添加条形图图层，映射direction为分组依据，设置堆叠方式为dodge，并行排列
g<-g+geom_bar(aes(fill=direction),stat="identity", position="dodge")

#设置分面，依Class按列分面，设置坐标scale与不同分面宽度space为自由调整
g<-g+facet_grid(.~Class,scales = "free",space = "free")

#设置y轴坐标分隔点
g<-g+scale_y_continuous(breaks=ybreaks)
#设置不同分组颜色
g<-g+scale_fill_manual(values=colors)
#设置标题
g<-g+labs(title=MainTitle,x=xLabel,y=yLabel,fill=legendTitle)

#设置外观主题
g<-g+ theme_bw()

#微调外观
g<-g+theme(
  #设置网格线颜色与大小
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  panel.grid.major.y = element_line(colour = "grey80",size=.25),
  
  #设置标题字体与居中
  plot.title=element_text(size=15,hjust=.5),
  
  #设置坐标轴字体，x轴label旋转
  axis.text.x = element_text(size=14,color='black',angle = 60,hjust = 1),
  axis.text.y = element_text(size=14,color="black"),
  
  #去除分标题（Class）背景
  strip.background=element_blank(),
  
  #设置分标题字体大小
  strip.text.x = element_text(size = 12, face = "bold"),
  
  #设置legend位置与标题字体
  legend.position="right", 
  legend.title=element_text(size=12,hjust=.5),
  )
}
g
graph2ppt(g,file='GO-114vsNIP.barplot.pptx')
dev.off()

