rm(list = ls())
options(stringsAsFactors = F)

library(tidyr)
library(XLConnect)
library(ggplot2)
library(pheatmap)
library(ggpubr)
library(export)

ego<- readWorksheetFromFile('ssRNAseq/GO_DEGs.xlsx',2)[,c(1,3,4,8)]
#Group  Ontology  GOterms  FDR
a1<-spread(ego,key='Group',value = 'FDR') #KEY为列名，Value为矩阵里面值

#dat[dat==0] = NA 
a<-a1[,3:ncol(a1)]
rownames(a)<-a1[,2]
#a$Ontology <- a1[,1]
#a[is.na(a)] <- 0

anno =data.frame(Ontology=a1[,1] , row.names = rownames(a))
table(anno)

write.csv(a,'GO_matrix_for_heatmap.csv')#去excel中修改列的顺序,并添加至GO_DEGs.xlsx。

aa<-read.csv('GO_matrix_for_heatmap.csv',row.names = 1)

aa = -log10(aa)
anno_col <- data.frame(Condition=c('Up','Down'),row.names=colnames(aa))
collab <-c('Up','Down')

p1<-pheatmap(aa, angle_col=270,border_color = 'black',
             #legend_breaks =c(0,2,4,6,8,10,12,14),kmeans_k = 3,
         #display_numbers = TRUE,number_color = "blue",
         #cellheight = 25,cellwidth = 25,
         show_rownames = T,na_col = 'white',
         cluster_rows = F,cluster_cols = F,fontsize=18,
         gaps_row = c(23, 26),#gaps_col = c(2,4,6),
         color=colorRampPalette(colors = c('white','orange'))(100)[30:100],
         annotation_row = anno,labels_col = collab,annotation_col=anno_col,
         annotation_names_row = F,annotation_names_col = F)
p1
#pdf(p1,file = 'a3.go.pdf')
graph2ppt(p1,file='GO.heatmap.pptx')

