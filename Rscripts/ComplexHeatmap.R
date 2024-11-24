rm(list = ls())
options(stringsAsFactors = F)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(XLConnect)
#unit(8,"cm") # 8cm
#=========================================================================ComplexHeatmap
mat <-read.csv('strand_specific/deg.uniq.fpkm_zscore.csv',row.names = 1)
mat$gene <-row.names(mat)
pheno<-data.frame(treat = factor(c(rep('CK',6),rep('LN',6))),
                  #time = factor(c(rep('30min',1),rep('60min',1),rep('24h',1),rep('48h',1),rep('96h',1),rep('8d',1),
                   #               rep('30min',1),rep('60min',1),rep('24h',1),rep('48h',1),rep('96h',1),rep('8d',1))),
                  row.names = colnames(mat)[1:12])
sheet <- readWorksheetFromFile('data/N-related_genes.xlsx',1)

gene_pos<-which(rownames(mat) %in% sheet[,1])
rownames(mat)[gene_pos]
anno <-c('OsNRT2.4','OsNRT2.3','OsAMT3.3','OsAMT1.2','OsAMT3.2',
         'OsAMT1.1','OsAMT2.2','OsAMT1.3','OsNRT1.1b','OsNRT2.1',
         'OsNRT2.2','OsNAR2.1','OsNPF2.4','OsNAR2.2')
row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, 
                                                 labels = anno))
ha = HeatmapAnnotation(df = pheno,#gp = gpar(col = "black")
                       show_legend = c(FALSE),show_annotation_name = F
                       #col = list(treat = c("CK" ="red","LN" ="blue"),#time = c(...))
                       #boxplot = anno_boxplot(mat,axis = TRUE, #注释轴
                        #                      gp = gpar(# fill = ifelse(mat > 0, "red", "green"))),
                         #                               col = rep(rainbow(6),2))),
                       #annotation_height = c(1, 1, 10)
                       #annotation_height = unit.c((unit(4, "cm"))*0.5, unit(4, "cm"), unit(2, "cm"))
                       #gap = unit(c(1, 2), "mm") 
                       #heatmap = anno_density(mat, type = "heatmap"),
                       #violin = anno_density(mat, type = "violin")
                      #  which = "row", width = unit(2, "cm")  #行注释及宽度。
                       )
ha_row <-rowAnnotation(df=pheno)  #直接行注释
collab = c(rep('0.5h',1),rep('1h',1),rep('24h',1),rep('48h',1),rep('96h',1),rep('8d',1),
           rep('0.5h',1),rep('1h',1),rep('24h',1),rep('48h',1),rep('96h',1),rep('8d',1))

#复杂注释:anno_points(),anno_barplot(),anno_boxplot(),anno_histogram(),anno_density(),anno_text()
#draw(ha)
#text<-do.call("paste0", rep(list(colnames(mat)), 1)) 
#ha_rot_tx= HeatmapAnnotation(text = anno_text(text, rot = 135, just = "left", 
 #                                             offset = unit(2, "mm")))
#draw(ha_rot_tx)
#全局设置
#names(ht_global_opt())
ht_global_opt(heatmap_row_names_gp = gpar(fontface = "bold"), fast_hclust = TRUE,
              heatmap_column_names_gp = gpar(fontsize = 16))
ht <-Heatmap(as.matrix(mat[,1:12]), cluster_rows = T,cluster_columns = F,
             show_heatmap_legend = T,
             show_row_names = FALSE, show_column_names = T,
        #cell_fun = function(j, i, x, y, width, height, fill) { #显示数字，两位小数,用于相关性
         # grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))},
        #name = ' ', #图例标题
        row_order = c('30min','60min','24h','48h','96h','8d'),
        width = unit(16, "cm"), #height = unit(22, "cm"), ##
        #column_km = 3,
        #split = 10, #按行层次聚类，为10类，与Kmeans冲突
        top_annotation = ha, #bottom_annotation =ha,
        right_annotation = row_anno, 
        col = colorRampPalette(colors = c('cyan1','black','yellow'))(10),# rev(rainbow(10)), #颜色
        row_dend_reorder = TRUE ,
        column_split =c( rep("CK",6),rep("LN", 6)),#row_split = rep(c("A", "B"), 9),
        column_title = " ",  
        #gap = unit(2, "mm"),#设置切分的间隔为5mm
        column_title_gp = gpar(fontsize = 16, fontface = "bold",col='black'), #设置标题形式
        column_names_gp = gpar(fontsize = 16), #gpar()设置
        row_dend_width  = unit(2, "cm"), row_names_side = "left", 
        column_names_side = "bottom", #column_dend_side = "bottom"，
        row_dend_side = "left", #聚类树位置
        #heatmap_legend_side = "bottom"
        column_names_rot = 45,
        row_title_gp = gpar(fill=rainbow(10),fontsize = 10, 
                            col = "black", border = rainbow(10)),
        #km = 10,# kmeans聚类数目,按行切分 row_km=2
        border = F,
        #annotation_legend_param = list(treat=list(title = "",   #注释图例设置
        #                             title_gp = gpar(fontsize = 14), ncol = 2, 
        # title_position = "topcenter")
        #                            labels_gp = gpar(fontsize = 10))),
        heatmap_legend_param = list(title = "z-score", title_gp = gpar(fontsize = 14),
                                    title_position ="leftcenter-rot",labels_gp = gpar(fontsize = 14),#热图bar更改
                           at = c(-4,-2, 0,2, 4), labels = c("-4",'-2', "0",'2', "4")) ,
        # width = unit(5, "cm"), #热图宽度固定
        # rect_gp = gpar(col = "black", lty = 1, lwd = 2), #设置热图边框
        # row_title = "I am a row title", row_title_rot = 0, #标题旋转角度
         split = 10,#data.frame(rep(c("A", "B"), 6), rep(c("C", "D"),6),#更通用的切分方式
        # col = colorRamp2(c(-3, 0, 3), c("green", "white", "red")), # col自定义颜色，colorRamp2来自circlize包
        # na_col = "orange",  #NA值颜色 
        # clustering_distance_rows = "pearson", clustering_method_rows = "single",
        show_row_dend = FALSE #是否显示聚类树，仍聚类
        # row_names_gp = gpar(col = c(rep("red", 4), rep("blue", 8)),
        # heatmap_legend_param = list(title = "legend" ,title_position = "topcenter",
                                    #legend_height=unit(8,"cm"), 
                                     #legend_direction="vertical") #修改图例
        
        # row_order = 12:1, column_order = order(as.numeric(gsub("LOC", "", colnames(mat)), #不聚类时，调整行列名称顺序
        # combined_name_fun = function(x) paste(x, collapse = "\n") # combined_name_fun = NULL 每个行切片的标题可以通过ta设定
        
        )  #+# +ha_row #添加行注释
           rowAnnotation(link = row_anno_link(at = gene_pos, labels = anno),#at为索引
                width = unit(1, "cm") + max_text_width(anno))

ht
dim(ht)
ht[1:2237,1:6] #sub heatmap
ht[1:2237,7:12]

#num = row_dend(ht) #分类的数目
idx = row_order(ht) #分类的索引

for (i in seq(1,10)){
  print(length(idx[[i]]))
  write.csv(mat[idx[[i]],] ,file = paste('cluster',i,'.csv',sep = '') ) }

column_order(ht)

#--------------------------------------------------------------------------------多个热图
library(ComplexHeatmap)

mat1 = matrix(rnorm(80, 2), 8, 10)
mat1 = rbind(mat1, matrix(rnorm(40, -2), 4, 10))
rownames(mat1) = paste0("R", 1:12)
colnames(mat1) = paste0("C", 1:10)

mat2 = matrix(rnorm(60, 2), 6, 10)
mat2 = rbind(mat2, matrix(rnorm(60, -2), 6, 10))
rownames(mat2) = paste0("R", 1:12)
colnames(mat2) = paste0("C", 1:10)

ht1 = Heatmap(mat1, name = "ht1",km=2)
ht2 = Heatmap(mat2, name = "ht2",width = unit(5, "cm"))

ht_list = ht1 + ht2 #两个热图合并在一张图里
class(ht_list)
draw(ht_list, row_title = "Two heatmaps, row title", row_title_gp = gpar(col = "red"),
     column_title = "Two heatmaps, column title", column_title_side = "bottom",
     gap = unit(1, "cm") ,#两图间隔
     main_heatmap = "ht1",#设置主热图
     show_row_dend =TRUE, cluster_rows = TRUE,
     heatmap_legend_side = "left", annotation_legend_side = "bottom",
     show_heatmap_legend = FALSE, show_annotation_legend = FALSE)
