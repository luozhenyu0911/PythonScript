rm(list=ls())
options(stringsAsFactors = F)


library(ggplot2)
library(RColorBrewer)  
library(reshape2) 
library(Hmisc)
#====================方法1, R基本函数cor=======================

#pearson correlation: 线性相关
#spearman correlation: 秩相关，秩序相关, 等级相关，rank correlation

#线性关系
#y=2x
#x=c(0,1,2,3,4,5)
#y=c(0,2,4,6,8,10)
#plot(x,y)
#cor(x,y,method='pearson')
#cor(x,y,method='spearman')

#非线性关系
#x=c(0,1,2,3,4,5)
#y=c(0,0.1,100,7,18,60)
#plot(x,y)
#cor(x,y,method='pearson')
#cor(x,y,method='spearman')

#表达量相关系数
#expr=read.delim('strand_specific/expr/sense_fpkm.txt',header=T,row.names = 1,
 #               sep = "\t",na.strings = '-') #读取表达矩阵
#r=cor(expr,method='pearson',use='pairwise')  
# method可选 pearson 或者spearman
# use="pairwise" 表示剔除缺失

#write.table(r,file='pearson.r.txt',col.names = NA,quote=F,sep="\t")

#=================方法2，Hmisc函数包===========================更好
#rm(list=ls())
#if(!requireNamespace('Hmisc',quietly = TRUE))  install.packages("Hmisc",update=F)

expr=read.delim('strand_specific/expr/sense_fpkm.txt',header=T,
                row.names = 1,sep = "\t",na.strings = '-') #读取表达矩阵

expr=expr[rowSums(expr)!=0,] #
res=rcorr(as.matrix(expr),type="pearson")
#type可选 peason 或者spearman
#rcorr函数默认会剔除缺失
res$r#相关系数
res$P#P值

#rcorr返回相关系数r/rho，以及相关性检验显著性p
#write.table(res$r,file='pearson.correlation.txt',col.names = NA,quote=F,sep="\t")
#write.table(res$P,file='pearson.pvalue.txt',col.names = NA,quote=F,sep="\t")

mydata <- melt(as.matrix(res$r))
colnames(mydata)<-c("S1","S2","corr")


#------------------------------------------气泡标签图-------------------------------------------------------
#if(!requireNamespace('corrplot',quietly = TRUE))  install.packages("corrplot",update=F)
library(corrplot)
library(export)

#设置颜色渐变方案
col<- colorRampPalette(c("white","red"))

#绘制简单corrplot气泡图，type有lower，upper，full三种
#绘制类型有"circle", "square", "ellipse", "number", "shade", "color", "pie"
#png("corr.simple.png",width=3000,height=2500,units = 'px',res=300)
corrplot(as.matrix(res$r), type = "full",
         method="pie",col=col(100))

dev.off()


#绘制混合corr气泡图，
#默认上三角默认为circle（upper="circle），下三角默认为number （lower="number"），也可自定义
#order为样本排序依据，alphabet为字母顺序，还可以"AOE", "FPC", "hclust"几种聚类方式排序
#mixed方式的颜色无法修改
#png("corr.mixed.png",width=3000,height=2500,units = 'px',res=300)
corrplot.mixed(as.matrix(res$r),lower = "number",order ="alphabet",
               pch.col = "black",  cl.lim = c(0, 1))
dev.off()



#------------------------------------------相关矩阵热图 --------------------------------------------------------
#设置ggplot基本坐标系，横轴为S1列，纵轴为S2列，颜色填充为corr列
g<-ggplot(mydata, aes(x = S1, y = S2, fill = corr))

#添加矩阵块图层，并设置边缘颜色
g<-g+geom_tile(colour="black") 

#添加文本图层，设置显示的label为corr列的值，字体大小，以及字体颜色
g<-g+geom_text(aes(label=corr),size=3,colour="black")

#设置横纵坐标大小一致，保证整体为正方形
g<-g+coord_equal()

#设置颜色填充，设置最小值颜色，最大值颜色，缺失值颜色，也可用scale_fill_gradient2设置中间值及其对应颜色
g<-g+scale_fill_gradient(low='white',high='red',na.value = "grey50")  
#g<-g+scale_fill_gradient2(low='white',high='red',mid='pink',
                           #midpoint = (max(mydata$corr)+min(mydata$corr))/2,na.value = "grey50")

#设置标题，以及x轴和y轴的标题
g<-g+labs(title="",x="",y="")


#设置其他选项，背景色，字体大小，标题居中，x与y轴样品名称字体颜色与字体大小，设置legend位置，文字大小，标题文字大小
g<-g+theme(panel.background=element_rect(fill="white",colour=NA),
      text=element_text(size=15),  
      plot.title=element_text(size=15,hjust=.5),
      axis.text.x = element_text(color="black", size=11, angle=90),
      axis.text.y = element_text(color="black", size=11, angle=0),
      legend.position="right",
      legend.text=element_text(size=10),
      legend.title=element_text(size=10)
)

#png("corr.matrix.heatmap.png",width=3000,height=2500,units = 'px',res=300) 
g
#dev.off()



#------------------------------------------相关性气泡图 --------------------------------------------------------
#设置ggplot基本坐标系，横轴为S1列，纵轴为S2列
g<-ggplot(mydata, aes(x= S1 , y=S2))

#添加点图层，使用corr列的值作为点大小及其颜色填充依据，设置点的形状，边缘颜色
g<-g+geom_point(aes(size=corr,fill = corr), shape=21, colour="black")

#设置横纵坐标大小一致，保证整体为正方形
g<-g+coord_equal()

#设置颜色填充，设置最小值颜色，最大值颜色，缺失值颜色，也可用scale_fill_gradient2设置中间值及其对应颜色
g<-g+scale_fill_gradient(low='white',high='red',na.value = "grey50")  
#g<-g+scale_fill_gradient2(low='white',high='red',mid='pink',midpoint = (max(mydata$corr)+min(mydata$corr))/2,na.value = "grey50",na.value=NA)

#设置气泡的最大值，及设置不显示size的legend
g<-g+scale_size_area(max_size=12, guide=FALSE) 
#设置标题，以及x轴和y轴的标题
g<-g+labs(title="",x="",y="")


#设置其他选项，背景色，字体大小，标题居中，x与y轴样品名称字体颜色与字体大小，设置legend位置，文字大小，标题文字大小
g<-g+theme(panel.background=element_rect(fill="white",colour=NA),
      text=element_text(size=15),  
      plot.title=element_text(size=15,hjust=.5),
      axis.text.x = element_text(color="black" ,vjust=0.5,size=11, angle=45),
      axis.text.y = element_text(color="black", size=11, angle=0),
      legend.position="right",
      legend.text=element_text(size=10),
      legend.title=element_text(size=10)
)

#png("corr.bubble.heatmap.png",width=3000,height=2500,units = 'px',res=300) 
g
#dev.off()
