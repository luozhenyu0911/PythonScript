#--------------------------------------------------------------------TCseq
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("TCseq")

#表观基因组学和转录组时间过程测序数据的定量和差异分析，
#时间过程数据的时间模式的聚类分析和可视化。

library(TCseq)
#svg(filename="heat.svg",width=3.5,height=15)
data <- read.csv('strand_specific/cluster10/deg.uniq.fpkm_zscore.csv',row.names = 1)[,1:12]
data <- as.matrix(data)
tca <- timeclust(data, algo = "cm",#选择聚类方法,别的聚类方法全是灰色
                 k = 3, #聚类数目
                 standardize = F # if TRUE, z-score transformation will performed 
)

p <- timeclustplot(tca, categories="time", #横坐标标签
                   value = "expression", #纵坐标标签
                   cols = 1,  #出图的列数
                   #membership.color = rainbow(30),
                   axis.text.size=16, legend.text.size=14,
                   legend.title.size=14, title.size=18,
                   axis.title.size=18 , axis.line.size = 0.6 )
p

# to plot a individual cluster
print (p[[3]]) # plot cluster 3
dev.off()

cluster1 <- clustCluster(tca) #输出基因的聚类信息

write.csv(cluster1,file = "heat_cluster.csv")

x <- matrix(sample(500, 1600, replace = TRUE), nrow = 200,
            dimnames = list(paste0('peak', 1:200), 1:8))
clust_res <- timeclust(x, algo = 'cm', k = 4, standardize = TRUE)
p <- timeclustplot(clust_res, cols =2)
# to plot a individual cluster
print (p[[2]]) # plot cluster 2
print (p[[3]]) # plot cluster 3
