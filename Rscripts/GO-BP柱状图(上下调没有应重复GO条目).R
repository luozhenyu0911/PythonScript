rm(list=ls())
options(stringsAsFactors = F)

library(ggplot2)
library(XLConnect)
library(export)
#分别读取上调与下调差异基因富集结果，包括三列，GO.ID	Term	Pvalue：
up<- readWorksheetFromFile('ssRNAseq/GO_DEGs.xlsx',4)[c(1:16),c(2,4,5,8)]
down<- readWorksheetFromFile('ssRNAseq/GO_DEGs.xlsx',5)[c(1:10),c(2,4,5,8)]


#对富集的显著性p值做-log10转换,其中下调为log10 转换
#上调为正数，下调为负数，以此区别上下调
up$FDR=-log10(up$FDR)
down$FDR=log10(down$FDR)

#将上下调结果进行合并
mydata=rbind(up,down)

#根据p值正负号添加两列：
#label1：上调的填写功能名称，下调的为空（NA）
#label2：上调的为空（NA）,下调的填写功能名称
mydata<-transform(mydata,
                  label1=ifelse(FDR>0,GOterms,NA),
                  label2=ifelse(FDR>0,NA,GOterms))

#按照显著性来进行排序，区分开上下调，上调从大到小，下调从小到大（负数）
mydata$GOterms <- factor(mydata$GOterms, 
                      levels = unique(mydata$GOterms[order(mydata$FDR)]))

#分别设置最大值，最小值颜色
highColor="red"
lowColor="white"

#计算显著性p值（log转换后）最大值和最小值
maxPvalue=ceiling(max(abs(mydata$FDR)))
minPvalue=-maxPvalue

#设置ggplot对象，以功能名称为x轴，以Pvalue(log转换后)为y轴，并以Count值作为颜色填充依据
g<-ggplot(data = mydata, aes(x = GOterms, y = FDR,fill = Counts)) 

#添加条形图图层，并设置y轴的范围
g<-g+geom_bar(stat = "identity", width = 0.8,colour="black",size=0.25) + 
        ylim(minPvalue,maxPvalue)


#设置颜色填充,最大值与最小值对应的颜色
g<-g+scale_fill_gradient(low=lowColor,high=highColor) 
#g<-g+scale_fill_gradient2(low=lowColor,midpoint=50,mid="pink",high=highColor) #还可以添加设置中点颜色

#转换x和y坐标
g<-g+coord_flip()


#额外添加功能名称到y=0位置，其中label1（上调标签）填写到，label2（下调标签）
#hjust微调显示位置（0为左对齐,<0可以增加与坐标轴的空格；1为右对齐，>1可以增加与坐标轴的空格; 0.5为居中）
g<-g+geom_text(aes(y = 0,label=label1),size=5,hjust= 1.1) #添加上调标签
g<-g+geom_text(aes(y = 0,label=label2),size=5,hjust=-0.1) #添加下调标签


#设置极简主题
g<-g+ theme_minimal()

#设置标题
g <- g + labs(title="",x="",y="-log10(Pvalue)")

#其他设置，设置网格线颜色与大小，设置标题字体与居中，设置坐标轴字体，去除多余的y轴功能名称，设置legend位置与字体大小
g<-g+theme(
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x = element_line(colour = "grey80",size=.25),
        panel.grid.minor.x = element_line(colour = "grey80",size=.25),
        plot.title=element_text(size=15,hjust=.5),
        axis.text.x = element_text(color="black", size=11, angle=0),
        axis.text.y = element_blank(), #左坐标轴字清除
        legend.position="right",
        legend.text=element_text(size=10),
        legend.title=element_text(size=10))
g
graph2ppt(g,file='GO-barplot.pptx')
dev.off()
