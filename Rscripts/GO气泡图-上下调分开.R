rm(list=ls())
options(stringsAsFactors = F)

library(ggplot2)
library(XLConnect)
library(export)
library(ggThemeAssist)
#------------------------------------------------------------------------UP regulation
#Data=read.delim('data.txt',header=T,stringsAsFactors =F)
if(T){
Data=readWorksheetFromFile('ssRNAseq/GO_DEGs.xlsx',4)
Data$FDR=-log10(Data$FDR)  #对p值取对数转换

#按照GeneRatio列进行排序
order<-order(Data$GeneRatio)
Data$GOterms<- factor(Data$GOterms, levels = Data$GOterms[order])

#设置label
MainTitle=''
xLabel='Gene Ratio'
yLabel=''
sizeLabel='Count'
colorLabel='-log10(P)'

#设置点的大小
DotRange=c(3,10)

#设置p值颜色，从小到大
cc = colorRampPalette(c('red','white'))
col=cc(100)[5:65]   #去除极颜色
col=rev(col)

colorbreaks=pretty(Data$FDR,5)
sizebreaks=pretty(Data$Counts,5)

#构建ggplot对象，添加坐标轴映射
g<-ggplot(data=Data,mapping=aes(GeneRatio,GOterms))

#添加点图层
g<-g+geom_point(aes(size=Counts, color=FDR),show.legend = T,stroke = 1) #stroke,点边界

#按照Ontology进行分面，并设置坐标与不同面的大小为free
g<-g+facet_grid(Ontology~.,scales= "free",space = "free")

#设置点大小控制属性
g<-g+scale_size(name=sizeLabel,breaks=sizebreaks,
                labels=sizebreaks,range = DotRange,guide = "legend")

#设置颜色属性
g<-g+ scale_color_gradientn(colours=col,breaks=colorbreaks)

#设置y轴位置
g<-g+scale_y_discrete(position = "left")

#设置标题
g<-g+labs(title=MainTitle,x=xLabel,y=yLabel,size=sizeLabel,color=colorLabel)

#设置主题，并微调部分外观
g<-g+theme_bw()
g<-g+theme(plot.title = element_text(hjust=0.5),
           strip.text = element_text(size = 16,
                                     colour = 'black'),#分面字体
           strip.background = element_rect(fill = 'grey', 
                                           linetype = 1),
           axis.title = element_text(size = 18),
           axis.text = element_text(colour="black", size=14),
           axis.text.y = element_text(lineheight = 0.5, size=18),
           legend.text=element_text(size=14),
           legend.title=element_text(size=16),
           legend.position="right")
}
g

graph2ppt(g,file='GO-114vsNIP.up.pptx')
dev.off()


#--------------------------------------------------------------- DOWN regulation

if(T){
  Data=readWorksheetFromFile('ssRNAseq/GO_DEGs.xlsx',5)
  Data$FDR=-log10(Data$FDR)  #对p值取对数转换
  
  #按照GeneRatio列进行排序
  order<-order(Data$GeneRatio)
  Data$GOterms<- factor(Data$GOterms, levels = Data$GOterms[order])
  
  #设置label
  MainTitle=''
  xLabel='Gene Ratio'
  yLabel=''
  sizeLabel='Count'
  colorLabel='-log10(P)'
  
  #设置点的大小
  DotRange=c(3,10)
  
  #设置p值颜色，从小到大
  cc = colorRampPalette(c('red','white'))
  col=cc(100)[5:65]   #去除极颜色
  col=rev(col)
  
  colorbreaks=pretty(Data$FDR,5)
  sizebreaks=pretty(Data$Counts,5)
  
  #构建ggplot对象，添加坐标轴映射
  gg<-ggplot(data=Data,mapping=aes(GeneRatio,GOterms))
  
  #添加点图层
  gg<-gg+geom_point(aes(size=Counts, color=FDR),show.legend = T,stroke = 1) #stroke,点边界
  
  #按照Ontology进行分面，并设置坐标与不同面的大小为free
  gg<-gg+facet_grid(Ontology~.,scales= "free",space = "free")
  
  #设置点大小控制属性
  gg<-gg+scale_size(name=sizeLabel,breaks=sizebreaks,
                  labels=sizebreaks,range = DotRange,guide = "legend")
  
  #设置颜色属性
  gg<-gg+ scale_color_gradientn(colours=col,breaks=colorbreaks)
  
  #设置y轴位置
  gg<-gg+scale_y_discrete(position = "left")
  
  #设置标题
  gg<-gg+labs(title=MainTitle,x=xLabel,y=yLabel,size=sizeLabel,color=colorLabel)
  
  #设置主题，并微调部分外观
  gg<-gg+theme_bw()
  gg<-gg+theme(plot.title = element_text(hjust=0.5),
               strip.text = element_text(size = 16,
                                         colour = 'black'),#分面字体
               strip.background = element_rect(fill = 'grey', 
                                               linetype = 1),
             axis.title = element_text(size = 18),
             axis.text = element_text(colour="black", size=14),
             axis.text.y = element_text(lineheight = 0.5, size=18),
             legend.text=element_text(size=14),
             legend.title=element_text(size=16),
             legend.position="right")
}
gg
#ggThemeAssistGadget(gg)
graph2ppt(gg,file='GO-114vsNIP.down.pptx')

dev.off()


####---------------------------------------------------------------------拼图

library(patchwork)
p <-g + gg
p  & scale_fill_continuous(limits = c(0, 13))
p+ plot_layout(guides = 'collect')
p