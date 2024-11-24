rm(list=ls())

library(XLConnect)
library(ggplot2)
library(export)

Edata<-readWorksheetFromFile('ssRNAseq/GO_DEGs.xlsx',2)
head(Edata)

#按照FDR列对GOTREMS进行排序
order<-order(Edata$FDR,decreasing = T)
Edata$GOterms<- factor(Edata$GOterms, levels =unique(Edata$GOterms[order]) )

#GO富集柱状图
p <- ggplot(data=Edata)+  theme_classic(base_size = 18)+
  geom_bar(aes(x=GOterms,y=-log10(FDR), fill=Ontology), 
           stat='identity') +
  theme(legend.position = c(-1.5,.4),
        axis.text.x=element_text(angle=0,size=14, vjust=0.5,colour = 'black'),
        axis.text.y=element_text(angle=0,size=18, vjust=0.5,colour = 'black'),
        axis.title.x = element_text(face="bold",color="black"),
        axis.title.y = element_text(face="bold",color="black"),
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =18),
        panel.background = element_rect(fill="transparent"),
        #panel.border =element_rect(fill='transparent', color='black')
        )+ coord_flip() + 
  labs(x = "", y = "-log10(pvalue)")+
  facet_grid(Ontology~. ,scales= "free",space = "free")
p
graph2ppt(p,file='GO-barplot.pptx')