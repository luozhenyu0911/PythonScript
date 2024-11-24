rm(list = ls())
options(stringsAsFactors = F)

library(reshape2)
library(XLConnect)
library(ggplot2)
library(ggpubr)
library(export)
windowsFonts()

num <-readWorksheetFromFile('e:/projects/rice/strand_specific/deg/deg-numbers.xlsx',2)
colnames(num) <-c('regulated','0.5h','1h','1d','2d','4d','8d')

num.m <-melt(num)
sizebreaks<-pretty(num.m$value,5)

p <- ggplot(data=num.m ,aes(x=variable,y=value,fill=regulated))+ 
  geom_bar(stat = "identity",position = "stack") + 
  #绘制条形图，position = "dodge"设置条形图不堆叠显示
  geom_text(aes(label=value),position=position_stack(),
            vjust=-0.5,color="black",size=5) + #在条形图上方0.5处(vjust=-0.5)以黑色(color="black")字体大小为5显示(size=5)数值大小
  #scale_fill_manual(values = brewer.pal(12, "Paired")[c(1,5)])+ #设置填充的颜色
  theme_classic(base_size = 18) +   #levels(c('30min','60min','24h','48h','96h','8d')) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  theme(legend.position = c(0.2,0.9),
        axis.text.x=element_text(#label=c('30min','60min','24h','48h','96h','8d'),
                                 angle=0,hjust = 0.5,colour="black",size=18), 
        axis.text.y=element_text(size=18,colour="black"), 
        axis.title.y=element_text(size = 20), #设置y轴标题的字体属性
        axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(colour="black", size=16),
        legend.title=element_blank() ,legend.background = element_rect(colour = 'black')) + 
  #scale_fill_discrete(name="", breaks = c("down", "up"),
   #                   labels = c("Down-regulated", "Up-regulated")) + 
  scale_fill_manual( values = c('#56B4E9','#ff4040'), 
                    name=" ", breaks = c("down", "up"),
                    labels = c("Down-regulated", "Up-regulated"))+
  ylab("Number of DEGs")+xlab("") 
p

graph2ppt(p,file='LN-CK-6times-barplot.pptx')
ggsave(p,filename = 'LN-CK-6times-barplot.png')



num1 <-readWorksheetFromFile('e:/projects/rice/strand_specific/deg/deg-numbers.xlsx',3)[,-c(7:10)]
name <-c('Groups','1h/0.5h','1d/0.5h','2d/0.5h','4d/0.5h','8d/0.5h')
colnames(num1) <-name
num.m1 <-melt(num1)

p1 <- ggplot(data=num.m1 ,aes(x=variable,y=value,fill=Groups))+ 
  geom_bar(stat = "identity",position = "stack") + #绘制条形图，position = "dodge"设置条形图不堆叠显示
  geom_text(aes(label=value),position=position_stack(),
            vjust=-0.5,color="black",size=5) + #在条形图上方0.5处(vjust=-0.5)以黑色(color="black")字体大小为5显示(size=5)数值大小
  #scale_fill_manual(values = brewer.pal(12, "Paired")[c(1,5)])+ #设置填充的颜色
  theme_classic(base_size = 18) +  
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  theme(legend.position = c(0.2,0.9),
        axis.text.x=element_text(
          angle=0,hjust = 0.5,colour="black",size=18), 
        axis.text.y=element_text(size=18,colour="black"), 
        axis.title.y=element_text(size = 20,face="plain"), #设置y轴标题的字体属性
        axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(colour="black", size=16),
        legend.title=element_blank() #,legend.background = element_rect(colour = 'black')
        ) + 
  #scale_fill_discrete(name="", breaks = c("down", "up"),
  #                   labels = c("Down-regulated", "Up-regulated")) + 
  scale_fill_manual( values = c('#1c86ee','#ff1493'), 
                     name=" ", breaks = c("down", "up"),
                     labels = c("Down-regulated", "Up-regulated"))+
  ylab("Number of DEGs")+xlab("") 
p1
graph2ppt(p1,file='LN-DEGs-barplot.pptx')
ggsave(p1,filename = 'LN-DEGs-barplot.png')
