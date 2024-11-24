# data<-log2(countdata+1)
# library(pheatmap) 
# choose_gene=head(rownames(countdata),1000) ## 50 maybe better 
# choose_matrix=countdata[choose_gene,] 
# #choose_matrix=t(scale(t(choose_matrix))) 
# pheatmap(data,show_rownames = F) 
# install.packages("factoextra")
deg.count<-read.table("D:/ARyuyan/R_data/cca1/deg.count")
deg.fpkm<-read.table("D:/ARyuyan/R_data/cca1/data/newdeg.fpkm")
countsdata<-deg.fpkm
countdata<-countsdata[,2:ncol(countsdata)]
rownames(countdata)<-deg.fpkm$V1
colnames(countdata)
rownames(countdata)

names(countdata)<-c("D4",
                    "D12",
                    "L4",
                    "L12")
library(factoextra)
# # bb<-b[,1:12]
# # df <- scale(log2(countdata+0.1))
datacount<-t(countdata)
ct<-scale(datacount)
 df<-t(ct)
fviz_nbclust(df, kmeans, method = "wss") + geom_vline(xintercept = 9, linetype = 2)
set.seed(123)
#利用k-mean是进行聚类
km_result <- kmeans(df, 16, nstart = 24)
#查看聚类的一些结果
print(km_result)
#提取类标签并且与原始数据进行合并
dd <- cbind(df, cluster = km_result$cluster)
dd<-data.frame(dd)
head(dd)
#查看每一类的数目
table(dd$cluster)
a<-data.frame(km_result$cluster)
#进行可视化展示
fviz_cluster(km_result, data = df,
             palette = c("Pink","Violet","BlueViolet","RoyalBlue","DeepSkyBlue","LightCoral","SpringGreen","Tomato","Gold","Chocolate","Magenta","SlateBlue","ForestGreen","Turquoise","Wheat","Red",""),
             # ellipse.type = "t",
             # star.plot = TRUE, 
             # repel = TRUE,
             # ellipse.type = "norm",
             stand = T,
             show.clust.cent = F,
             shape = "●",
             main = "",
             xlab = "", ylab = "",
             ellipse=F,
             pointsize=2,
             geom="point",
             ggtheme = theme_classic(),
)

?fviz_cluster
a<-dd[which(dd$cluster == 1), ]

write.table(a,file = "D:/ARyuyan/R_data/cca1/data/15",sep = "\t",row.names = T,col.names = T,quote = F)

