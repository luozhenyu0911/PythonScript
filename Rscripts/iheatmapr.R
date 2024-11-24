rm(list = ls())
options(stringsAsFactors = F)


#载入绘图包
library(iheatmapr)
library(datasets)
library(reshape2)
#使用acast调用Indometh数据集中的内容
Indometh_matrix <- acast(Indometh, Subject ~ time, value.var = "conc")
Indometh_matrix <- Indometh_matrix[as.character(1:6),]
rownames(Indometh_matrix) <- paste("Patient",rownames(Indometh_matrix))
##计算相关性矩阵
Indometh_patient_cor <- cor(t(Indometh_matrix))
##取每个样本数据中的最大值和最小值
patient_max_conc <- apply(Indometh_matrix,1,max)
patient_min_conc <- apply(Indometh_matrix,1,min)
##给每个样本随机分配一个分组
patient_groups <- c("A","A","B","A","B","A")


####绘制相关性矩阵热图
example_heatmap <- main_heatmap(Indometh_patient_cor, name = "Correlation")
example_heatmap
#可以使用colors参数修改颜色，colors = “RColorBrewer palettle/颜色名字”。
main_heatmap(Indometh_patient_cor,
             colors = "Pinks",
             name = "Correlation")



##添加第二幅热图
main_heatmap(Indometh_patient_cor, name = "Correlation") %>%
  add_main_heatmap(Indometh_matrix, name = "Indometacin<br>Concentration"
#如果name=相同名字的话，新添加的热图会与之前的热图共享相同的颜色。
#如果你想改变添加热图的位置，就使用side=c(“left”, “right”, “top”,”bottom”)。
                   


###添加标签
main_heatmap(Indometh_matrix, name="Correlation") %>%
                     add_row_labels() %>%
                     add_col_labels() %>%
                     add_row_title("Patients", buffer=0.2) %>%
                     add_col_title("Patients", buffer=0.1)
#buffer规定title文字与图之间的距离。
                   


###添加聚类关系
main_heatmap(Indometh_patient_cor) %>%
                     add_row_clustering() %>%
                     add_col_clustering()
#如果是希望以颜色表示样本的聚类关系，则需要先进行一个聚类，之后在手动进行负值。
clust_assign <- kmeans(Indometh_matrix, 3)$cluster
                   
 main_heatmap(Indometh_patient_cor) %>%
                     add_row_clusters(clust_assign) %>%
                     add_col_clusters(clust_assign)
   
                 
###添加样本注释
main_heatmap(Indometh_patient_cor) %>%
                     add_row_annotation(data.frame("Max" = patient_max_conc,
                                                   "Min" = patient_min_conc,
                                                   "Groups" = c("A","A","B","B","A","A")),
                                        colors = list("Max" = "Reds",
                                                      "Min" = "Blues", "Groups" = c("purple","pink")))
                   
#除了add_row_annotation，还可以使用add_row_signal和add_row_groups添加注释。

main_heatmap(Indometh_patient_cor) %>%
                     add_row_signal(patient_max_conc, "Max<br>Concentration", title = "Max", colors = "Reds") %>%
                     add_row_signal(patient_min_conc, "Min<br>Concentration", title = "Min", colors = "Reds") %>%
                     add_row_groups(c("A","A","B","B","A","A"), "Groups")
   
                   
                   
                                   
###给出一个完整的代码。
main_heatmap(Indometh_patient_cor,name = "Correlation") %>%
                     add_col_clustering() %>%
                     add_row_clustering(k = 3) %>%
                     add_row_title("Patients") %>%
                     add_col_title("Patients") %>%
                     add_row_annotation(data.frame("Max" = patient_max_conc,
                                                   "Min" = patient_min_conc,
                                                   "Groups" = patient_groups)) %>%
add_main_heatmap(Indometh_matrix,
                name = "Indometacin<br>Concentration") %>%
                     add_col_labels() %>%
                     add_col_title("Time") %>%
                     add_col_summary()
