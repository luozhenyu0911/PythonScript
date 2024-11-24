# kmeans一般在数据分析前期使用，选取适当的k，将数据聚类后，然后研究不同聚类下数据的特点。
#step1：导入一组具有n个对象的数据集，给出聚类个数k；
#step2：从n个对象中随机取出k个作为初始聚类中心；
#step3：根据欧几里得距离来判断相似度量，确定每个对象数据哪个簇；
#step4：计算并更新每个簇中对象的平均值，并将其定为每个簇的新的聚类中心；
#step5：计算出准则函数E；
#step6：循环step3，step4，step5直到准则函数E在允许的误差范围内；
# 轮廓系数结合了聚类的凝聚度和分离度，用于评估聚类效果。该值处于-1~1间，值越大，表聚类效果越好。

#求最佳K值的方法，如：
#轮廓系数;
#寻找SSE的拐点;
#通过分割算法来估算类的中心;
#Calinsky标准：诊断多少集群适应数据;
#贝叶斯信息准则： 基于层次聚类的期望最大化;
#基于Affinity propagation (AP)聚类;
#基于Gap Statistic.

#----------------------------------------------------------------##Z-score标准化
data_scale <- as.data.frame(t(apply(data_df,1,scale))) 
names(data_scale) <- names(data_df)
#----------------------------------------------------------------------------


rm(list=ls())
options(stringsAsFactors = F)

# data_df=t(data_df) 转置，可以对样品进行聚类
data_df = read.csv("strand_specific/cluster10/deg.uniq.fpkm.csv", row.names = 1)
#data_scale <- as.data.frame(t(apply(data_df,1,scale))) 
#names(data_scale) <- names(data_df)

#zscore <-data.frame(t(scale(t(data_df))))
#zscore[is.na(zscore)] <- 0    # NA -> 0
write.csv(zscore,file = 'deg.uniq.fpkm.csv')

data <- read.csv("strand_specific/cluster10/deg.uniq.fpkm_zscore.csv", row.names = 1)[,1:12]


#K值确定有以下方法，最常用：计算轮廓系数。
#================================================================寻找SSE的拐点
#SSE（sum of the squared error），即所有点到其所属簇中心的距离的平方和即误差的平方和
wss <- (nrow(data)-1)*sum(apply(data,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(data,centers=i)$withinss)
plot(1:20, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
#从图看出，拐点值最佳（在这拐点前，函数值下降得很快，在该点后，函数值下降很慢）

#===============================================================基于Gap Statistic
library(cluster)
gap_stat <- clusGap(data, FUN = kmeans, nstart = 25, K.max = 20, B = 300, 
                    verbose = interactive())
plot(1:20, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
#从图看出，拐点值最佳.

#=================================================================计算轮廓系数,最常用。
# fpc包用来统计不同p值的轮廓系数
#install.packages('fpc')
library(fpc)
K <- 2:20   #设置遍历得K值范围自定义，计算每个K对应得轮廓系数，一般得K值都不会太大。
round <- 30 # 每次迭代20次，可自定义，取均值，避免误差；结果显示没有收敛的时候，可以增加迭代次数
bestK <- sapply(K, function(i){    # sapply循环遍历某个集合内的所有或部分元素
  print(paste("K=",i))
  mean(sapply(1:round, function(r){
    print(paste("Round",r))
    result <- kmeans(data,i)
    stats <- cluster.stats(dist(data),result$cluster)  #cluster.stats是fpc包统计轮廓系数得函数
    stats$avg.silwidth
  }))
})

# 展示不同K值对应的轮廓系数
plot(K,bestK,type='l',main = '轮廓系数与K的关系',ylab='轮廓系数')
#选择最优的轮廓系数（最大值对应的K）,也可以根据结果自己选一个
Select_K <- K[match(max(bestK),bestK)]
Select_K
#进行kmeans聚类
result <- kmeans(data,Select_K)
Group = result$cluster    #最终分类结果
Group

# 1.====================================================================根据上述结果做热图
library(ComplexHeatmap)
library(circlize)
#按照分类结果将分到一组里面的基因进行排序，方便画图

ind <- order(Group)  # index确定了每个基因排序后的位置
# expr <- log2(data[ind,]+1)   # 输入fpkm，以log2标准化
#建立Group的dataframe，用于绘图
Group=as.data.frame(as.factor(Group[ind]))
colnames(Group)=c('Group')
expr <-  data[ind,]   # ComplexHeatmap要求输入格式为矩阵

#设置指示条, ComplexHeatmap可以用+将指示条不同热图直接合并
ha = rowAnnotation(df = Group,width = unit(5,"mm"),
                   col=list(Group=setNames(rainbow(Select_K),levels(Group$Group))))

#设置热图
hexp=Heatmap(as.matrix(data), cluster_rows = F,cluster_columns = F,
             color=cc(100), #column_title = "Gene Expression",
             column_title_gp = gpar(fontsize=20),right_annotation = ha, 
             heatmap_legend_param = list(title = "FPKM Z-score"),
             #width = unit(100,"mm"),na_col = "grey",
             show_column_names=T,show_row_names=F)
hexp
#png('draw.png',h=1000,w=500,res=72,units='px')
#draw(ha+hexp,gap = unit(c(1,5),"mm"))

# 2.=====================================================================可以直接km参数分组
p <- Heatmap(as.matrix(data),name = "FPKM z-score", 
             km = 9,# km_title = "%i",column_title = "Gene Expression",
             column_names_side = "bottom",
             col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
             row_dend_side = "left",
             show_row_names = FALSE)
p

# Combined data
cls <- plyr::ldply(row_order(p), data.frame)  #获取聚类分组信息
names(cls) <- c("id", "rowid")
cls <- dplyr::arrange(cls, rowid) #按照rowid排序
data$cluster <- cls$id 

#用均值拟合
res <- plyr::ddply(data, "cluster", function(x){ 
  colMeans(x[,1:12]) }) #样品列数

data$group <- row.names(data)
res$group <- "mean"

data2 <- dplyr::bind_rows(data, res)

df <- data.table::melt(data2, id = 13:14, #最后两列
                       variable.name = "sample")
df$col <- ifelse(df$group == "mean", "red", "grey")
df$cluster <- as.numeric(df$cluster)

# Generate plot 
library(ggplot2)
library(grid)

grid.newpage()
pushViewport(viewport(layout = grid.layout(3,3))) #出图，根据K确定几行几列

#如何实现多行一图？？？
vplayout = function(x,y){
  viewport(layout.pos.row = x,layout.pos.col = y)}

#9为K值 ,批量出图。
for (i in 1:9){ 
  p<- ggplot(dplyr::filter(df, cluster == i), 
              aes(x = sample, y = value, group = group, col = col)) +
    geom_line() + theme_light(base_size = 18) +
    scale_colour_manual(values = c("grey", "red")) +
    guides(col = FALSE) #legend取消
  ggsave(p,filename = paste('cluster-',i,'.plot.png',sep = ''))
  #print(p,vp = vplayout(1,i))  
}
p
