a<-dd[which(dd$cluster == 9), ]
targt<-a[which(rownames(a)%in%c("Seita.1G236100.v2.2","Seita.6G055700.v2.2")),]
write.table(a,file = "D:/ARyuyan/R_data/cca1/9",sep = "\t",row.names = T,col.names = T,quote = F)

c2<-read.table("D:/ARyuyan/R_data/cca1/data/2")
c2<-c2[c2$V3<c2$V2,]
write.table(c2,file = "D:/ARyuyan/R_data/cca1/data/new2",sep = "\t",row.names = T,col.names = T,quote = F)
c3<-read.table("D:/ARyuyan/R_data/cca1/data/3")
c3<-c3[c3$V4>c3$V5,]
write.table(c3,file = "D:/ARyuyan/R_data/cca1/data/new3",sep = "\t",row.names = T,col.names = T,quote = F)
names( countdata)<- c("T1Dark4h","T2Dark12h","T3Light4h","T4Light12h")
library(pheatmap)
countdata<-correlatedA[,c(6:8,9:11)]

countdata<-scale(t(countdata))
countdata[is.na (countdata)] <- 0
pheatmap(countdata,legend_labels=F,show_rownames = F,
         cluster_row = FALSE,cluster_col = FALSE,
         main="Cluster2 Standardized FPKM",angle_col=0,
         border_color = NA,
         legend = F,
         color = colorRampPalette(colors = c("green","black","blue"))(100))
?pheatmap
