rm(list = ls())
options(stringsAsFactors = F)

library(tidyr)
library(XLConnect)
library(ggplot2)
library(pheatmap)
library(export)
library(RColorBrewer)

ego1<- readWorksheetFromFile('e:/projects/rice/strand_specific/go/KEGG.xlsx',1)[,c(1,3,7)]
a1<-spread(ego1,key='Group',value = 'FDR')
#dat[dat==0] = NA 
a<-a1[,2:ncol(a1)]
rownames(a)<-a1[,1]

write.csv(a,'kegg-heatmap.csv')
a<-read.csv('kegg-heatmap.csv',row.names = 1)
a = -log10(a)
collab <-c('Up','Down','Up','Down','Up','Down','Up','Down')
ann<-data.frame(Time=c(rep('24h',2),rep('48h',2),rep('4d',2),rep('8d',2)),
                row.names=colnames(a))
p1<-pheatmap(a, border_color = 'black',#angle_col=315,
             cellheight = 25,cellwidth = 25,na_col = 'white',
             cluster_rows = F,cluster_cols = F,
             fontsize=18,show_rownames = T, 
             color=colorRampPalette(colors = c('white','red'))(100)[40:100] ,
             gaps_col = c(2,4,6),labels_col = collab,
             annotation_col = ann)
p1
graph2ppt(p1,file='kegg-ln_ck.heatmap.pptx',width=11,height=10.2)
#pdf(p1,file = 'a3.go.pdf')
dev.off()


ego<- readWorksheetFromFile('e:/projects/rice/strand_specific/go/KEGG.xlsx',2)[,c(1,3,7)]
a2<-spread(ego,key='Group',value = 'FDR')
#dat[dat==0] = NA 
a3<-a2[,2:ncol(a2)]
rownames(a3)<-a2[,1]

write.csv(a3,'kegg-ln.heatmap.csv')
a3<-read.csv('kegg-ln.heatmap.csv',row.names = 1)
a3 = -log10(a3)
collab <-c('Up','Down','Up','Down','Up','Down','Up','Down')
ann<-data.frame(Time=c(rep('24h/30min',2),rep('48h/30min',2),
                       rep('4d/30min',2),rep('8d/30min',2)),
                row.names=colnames(a3))
p<-pheatmap(a3, border_color = 'black',#angle_col=315,
             cellheight = 25,cellwidth = 25,
             na_col = 'white',
             cluster_rows = F,cluster_cols = F,
             fontsize=18,show_rownames = T, 
             color=colorRampPalette(colors = c('white','red'))(100)[40:100] ,
             gaps_col = c(2,4,6),
             labels_col = collab,
             annotation_col = ann)
p
graph2ppt(p,file='kegg-ln.heatmap.pptx',width=11.6,height=7.7)

