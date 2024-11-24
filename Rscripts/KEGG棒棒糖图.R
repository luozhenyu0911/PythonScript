rm(list = ls())
options(stringsAsFactors = F)

library(ggplot2)
library(ggthemes)
library(XLConnect)
library(export)

tab<-readWorksheetFromFile('ssRNAseq/GO_DEGs.xlsx',3)

colnames(tab)
head(tab)

tab<-tab[order(tab$FDR,decreasing = T),]
tab<-tab[order(tab$Group,decreasing = F),]
tab$Description<-factor(tab$Description,
                        levels=unique(as.character(tab$Description)))
tab

p<-ggplot(tab,aes(x=Description,y=-log10(FDR)))+
  geom_bar(stat="identity",width=0.1,fill='grey60')+
  geom_point(aes(color=Group),size=11)+ #可选
  geom_text(aes(label=Counts),alpha=I(0.8))+ #Counts
  theme_bw(base_size = 18) +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    # element_line(size = 0.8,color="darkgray"), # element_blank(),
    axis.line.x = element_line(colour = "black", size = 0.8),
    axis.line.y = element_line(colour = "black", size = 0.8),
    axis.ticks.x = element_line(size = 0.8),
    axis.ticks.y = element_line(size = 0.8),
    axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0.5),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
    #  legend.position="NA",
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold"),
    legend.background = element_rect(fill = "transparent"),
    strip.background = element_rect(colour = "white", 
                                    fill = "white",size = 0.2),
    strip.text.x = element_text(),
    strip.text.y = element_text(),
    #text = element_text(size = 18, face = "bold",colour = 'black'),
    plot.title = element_text(face = "bold",colour = 'black'))+
  scale_color_pander()+
  xlab("KEGG Pathway")+ylab("-log10(FDR)")+
  coord_flip() 

p

graph2ppt(file='KEGG-group.pptx')
