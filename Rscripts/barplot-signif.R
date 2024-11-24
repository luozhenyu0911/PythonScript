#柱状图是不能用ggsignif求显著性的，但还可以用人为标注以示说明.

rm(list = ls())
options(stringsAsFactors = F)

library(ggsignif)
library(reshape2)
library(XLConnect)
library(ggplot2)
library(ggpubr)
library(export)
windowsFonts()

num <-readWorksheetFromFile('e:/projects/rice-LN-RNAseq/strand_specific/pcr-ggplot.xlsx',1,
                            rownames=1)
colnames(num) <-c('1h','1d','2d','4d','8d')

num.m <-melt(num)

num.m$condition <- c("CK1h", "LN1h","CK1d", "LN1d","CK2d", "LN2d",
                     "CK4d", "LN4d","CK8d", "LN8d")
num.m$trt <- c(rep(c('CK','LN'),5))

class(num.m[,2])

sizebreaks<-pretty(num.m$value,5)

my_comparisons <- list(c("CK1h", "LN1h"),c("CK1d", "LN1d"),c("CK2d", "LN2d"),
                       c("CK4d", "LN4d"),c("CK8d", "LN8d"))

p <- ggplot(data=num.m ,aes(x=variable,y=value,fill=trt))+ 
  geom_bar(stat = "identity",position = "dodge") + 
  theme_classic(base_size = 18) +  
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  theme(#legend.position = c(0.2,0.9),
        axis.text.x=element_text(angle=0,hjust = 0.5,colour="black",size=18), 
        axis.text.y=element_text(size=18,colour="black"), 
        axis.title.y=element_text(size = 20), 
        axis.line = element_line(colour = "black",size=1),
        legend.text=element_text(colour="black", size=16),
        legend.title=element_blank() ,legend.background = element_rect(colour = 'black')) + 
  scale_fill_manual(values = c('#56B4E9','#ff4040'))+
  ylab("Relative Expression Level") + xlab("") +
  #geom_signif(comparisons = my_comparisons, step_increase = 0.01,
   #           map_signif_level = F,test =wilcox.test )
  geom_signif(annotations = c(0.001, "***"), y_position = c(1.2,1.2), 
              xmin = c(0.75, 1.75), xmax = c(1.25, 2.25), 
            tip_length = c(c(0.02, 0.02),c(0.02, 0.02)), vjust = -1)
p

graph2ppt(p,file='LN-CK-6times-barplot.pptx')
ggsave(p,filename = 'LN-CK-6times-barplot.png')
