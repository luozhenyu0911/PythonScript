rm(list = ls())
options(stringsAsFactors = F)

library(tidyr)
library(XLConnect)
library(ggplot2)
library(pheatmap)

ego1<- readWorksheetFromFile('DEGs/go_kegg/GO.xlsx',10)[1:46,][,c(1,3,4,11)]
a3<-spread(ego1,key='Group',value = 'log')
#dat[dat==0] = NA 
a<-a3[,3:ncol(a3)]
rownames(a)<-a3[,2]
a[is.na(a)] <- 0
anno =data.frame(Ontology=c(rep('BP',21),rep('CC',3),rep('MF',nrow(a)-24)),
                 row.names = rownames(a))
#ann_colors = list(#group = c("white", "firebrick"),
        # Ontology= c(BP = "blue", CC = "yellow",MF='firebrick'))
p1<-pheatmap(a, angle_col=315,border_color = 'black',#legend_breaks =c(0,2,4,6,8,10,12,14),kmeans_k = 3,
         #display_numbers = TRUE,number_color = "blue",cellwidth = 25, cellheight = 12,
         cluster_rows = F,cluster_cols = F,fontsize=14,show_rownames = T, gaps_row = c(21, 24),
         color=colorRampPalette(colors = c('white','firebrick'))(100),gaps_col = c(1,3),
         annotation_row = anno,#annotation_colors = ann_colors
         )
p1
#pdf(p1,file = 'a3.go.pdf')
dev.off()


ego2<- readWorksheetFromFile('DEGs/go_kegg/GO.xlsx',10)[47:130,][,c(1,3,4,11)]
b3<-spread(ego2,key='Group',value = 'log')
#dat[dat==0] = NA 
b<-b3[,3:ncol(b3)]
rownames(b)<-b3[,2]
b[is.na(b)] <- 0
anno =data.frame(Ontology=c(rep('BP',27),rep('CC',18),rep('MF',nrow(b)-45)),
                 row.names = rownames(b))
#annocol =data.frame(Group=colnames(b3)[3:length(colnames(b3))],
 #                   row.names = colnames(b3)[3:length(colnames(b3))])
ann_colors = list(#group = c("white", "firebrick"),
 Ontology= c(BP = "green", CC = "lightblue",MF='purple'))
p2<-pheatmap(b, angle_col=315,#legend_breaks =c(0,2,4,6,8,10,12,14),kmeans_k = 3,
             #display_numbers = TRUE,number_color = "blue",cellwidth = 25, cellheight = 12,
             cluster_rows = F,cluster_cols = F,fontsize=14,show_rownames = T, gaps_row = c(27, 45),
             color=colorRampPalette(colors = c('white','red'))(100),gaps_col = c(1,3),
             annotation_row = anno,annotation_colors = ann_colors,
             #annotation_col = annocol,
             border_color = 'black')
p2

ego3<- readWorksheetFromFile('DEGs/go_kegg/GO.xlsx',10)[131:214,][,c(1,3,4,11)]
c3<-spread(ego3,key='Group',value = 'log')
#dat[dat==0] = NA 
c<-c3[,3:ncol(c3)]
rownames(c)<-c3[,2]
c[is.na(c)] <- 0
anno =data.frame(Ontology=c(rep('BP',26),rep('CC',14),rep('MF',nrow(c)-40)),
                 row.names = rownames(c))
#ann_colors = list(#group = c("white", "firebrick"),
# Ontology= c(BP = "blue", CC = "yellow",MF='firebrick'))
p3<-pheatmap(c, angle_col=315,border_color = 'black',#legend_breaks =c(0,2,4,6,8,10,12,14),kmeans_k = 3,
             #display_numbers = TRUE,number_color = "blue",cellwidth = 25, cellheight = 12,
             cluster_rows = F,cluster_cols = F,fontsize=14,show_rownames = T, gaps_row = c(26, 40),
             color=colorRampPalette(colors = c('white','navy'))(100),gaps_col = c(2,4),
             annotation_row = anno,#annotation_colors = ann_colors
)
p3
#d <-merge(b,c, by='row.names',sort=F)
