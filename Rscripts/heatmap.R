rm(list = ls())
options(stringsAsFactors = F)

library(pheatmap)
library(RColorBrewer)
display.brewer.all()
library(circlize)

#---------------------------------------------------------------------------------大量基因表达热图
display.brewer.all() #所有配色
cc<-colorRampPalette(rev(brewer.pal(n=7,name = 'RdYlBu')))#配色

#若不聚类可以允许有缺失值，若进行聚类应该去除缺失或替换NA。
fpkm<-read.csv(file = 'strand_specific/cluster10/all.DEGs.uniq.FPKM_z-score.csv', # na.strings = '-',
              row.names = 1)

#tpm<-read.csv('ln.ck.24h-up.anti.csv',row.names = 1)
#tpm$fc <-tpm$LN.24h/tpm$CK.24h
#t <-t(tpm)
#s <-scale(t,center = T,scale = T)
#tpm_z <- t(scale(t(tpm)))

head(fpkm)
fpkm <-na.omit(fpkm)
table(fpkm$cluster)
#fpkm <- fpkm[rowSums(fpkm)>0,]

#固定顺序
fpkm$cluster<- factor(fpkm$cluster,levels= paste0("cluster",1:10),
                      ordered = T)
table(fpkm$cluster)

#若为FPKM，一般log标准化。
#matrix<-log2(fpkm+1)
#matrix<-scale(fpkm,center = T,scale = T) #Z-score
cluster <-data.frame(Cluster=c(rep('C1',632),rep('C2',685),rep('C3',82),rep('C4',210),rep('C5',89),
                               rep('C6',291),rep('C7',58),rep('C8',99),rep('C9',73),rep('C10',18)),
                     row.names = rownames(fpkm))

pheno<-data.frame(Treatment = factor(c(rep('CK',6),rep('LN',6))),
                  row.names = colnames(fpkm)[1:12])
collab = c(rep('0.5h',1),rep('1h',1),rep('1d',1),rep('2d',1),rep('4d',1),rep('8d',1),
           rep('0.5h',1),rep('1h',1),rep('1d',1),rep('2d',1),rep('4d',1),rep('8d',1))
label <-fpkm$anno
#ann_colors = list(cluster = c(rainbow(15)))
#png('heatmap.png',height = 800,width = 800,units = 'px',res = 75)
ph<-pheatmap(fpkm[,1:12], cluster_rows=F, cluster_cols=F,
             angle_col = 0, gaps_col = c(6),border_color = 'black',
             gaps_row = c(632,1317,1399,1609,1698,1989,2047,2146,2219),
             color = colorRampPalette(colors = c('cyan1','black','yellow'))(10),#c('#436eee','white','#ee0000')
             #color =  rev(rainbow(10)),#kmeans_k = 12, 
             legend = T, legend_labels = c('Low','High'),#cutree_cols = 2, 
             #na_col = 'grey',  #缺失值画为grey颜色。
             #cellwidth = 35,
             display_numbers = F,#scale = 'row',
             fontsize=20, #main = 'DEG Expression Heatmap',
             annotation_row = cluster,
             annotation_col=pheno,
             drop_levels=T,
             show_rownames=F, number_format = "%.2f", 
             number_color = "grey30", fontsize_number = 0.8* fontsize,
             legend_breaks = c(-1.5,2.5),
             #annotation_legend = T,
             #annotation_colors = ann_colors,
             annotation_names_row = F,annotation_names_col = F,
             labels_col = collab,#labels_row = label,
             #filename = NA, width = NA, height = NA,
             )
ph
library(export)
graph2ppt(file='pheatmap-cluster10.pptx')
dev.off()
#png('heatmap.png',height = 800,width = 800,units = 'px',res = 75)
ph<-pheatmap(fpkm[,1:12], cluster_cols=T, cluster_rows=T,
             angle_col = 45, gaps_col = c(6),border_color = 'black',
             #gaps_row = c(632,1317,1399,1609,1698,1989,2047,2146,2219),
             color = colorRampPalette(colors = c('cyan1','black','yellow'))(10),#c('#436eee','white','#ee0000')
             #color =  rev(rainbow(10)),
             #kmeans_k = 14, 
             #scale = 'row',
             legend = T, legend_labels = c('Low','High'), 
             legend_breaks = c(-1.5,2.5),
             cutree_rows = 10, cutree_cols = 4,
             #na_col = 'grey',  
             cellwidth = 35,display_numbers = F,
             fontsize=20, #main = 'DEG Expression Heatmap',
             annotation_col=pheno,#annotation_row = cluster,
             show_rownames=F,
             #annotation_legend = T,#annotation_colors = ann_colors,
             annotation_names_row = F,annotation_names_col = F,
             labels_col = collab,#labels_row = label
)
ph

# 行的聚类排列顺序
ph$tree_row$order 
# 得到行名的顺序
rownames(fpkm)[ph$tree_row$order]
# 查看按行聚类后的热图顺序结果
head(fpkm[ph$tree_row$order,])
# 查看按照行和列聚类之后，得到的热图顺序结果
head(fpkm[ph$tree_row$order,ph$tree_col$order])

#按聚类分割行或列,并保存。
member <-cutree(ph$tree_row,k=10) #自定义K
write.table(member,file = 'gene.10group.xls',sep='\t',quote = F)

#还可对差异gene的Log2FC做热图。列为AvsB_log2FC,行为gene。
##还可对GO富集的FDR做热图。列为FDR,行为GOterms。一般为 -log10(FDR)。






#---------------------------------------------------------------------------------------其他，循环
library(XLConnect)
#循环
for (id in seq(1,3)){
  sheet <- readWorksheetFromFile('brown-blue-1238-fpkm.xlsx',id)
  row.names(sheet) <-sheet[,1]
  #cell =data.frame(high = factor(sheet$cell),row.names = rownames(sheet))
  sample =data.frame(sample = factor(c(rep('BS',ncol(sheet[,-c(1)])/2),rep('M',ncol(sheet[,-c(1)])/2)),levels = c('BS','M')),
                   row.names = names(sheet[,-c(1)]))
  heatmp<-paste('heatmap',id,sep='')
  png(paste(heatmp,'.png',sep=''),width = 1200,height = 800)
  pheatmap(log2(sheet[,-c(1)]+1),angle_col=315,annotation_col = sample,show_rownames = F,
           cluster_rows = F,cluster_cols = F,fontsize=14,
           color=colorRampPalette(colors = c('#436eee','white','#ee0000'))(100))
  dev.off()
}





#============================================================================相关性cor heatmap plot 
#RNAseq样本间重复相关性分析可用FPKM,也可用read counts计算。
library(corrplot)
cor_matrix<-read.table(file = 'PearsonCorr-114_Nip-DRIP.txt', # na.strings = '-',
                 header = T,sep = '\t',quote = '',check.names = F,row.names = 1)  

#cor_matrix <- cor(fpkm)

#png('pheatmap-sample.png',width = 1800,height = 1200)
#par(mar=c(1, 1, 1, 10))

pheatmap(cor_matrix,fontsize=16,angle_col=270, 
         # main='sample heatmap',
         color=colorRampPalette(colors = c('#0408f2','#ffff00'))(100),
         #cellwidth = 20,cellheight = 20,
         display_numbers = T ,#annotation_col =condition,
         annotation_legend = T,
         #gaps_col = c(2,6),
         border_color = NA)
library(export)
graph2ppt(file="DRIPseq-sample-correlation.pptx")

require(ggcor)
corr <-cor(exprSet)
df <- as_cor_tbl(corr)   #as_cor_tbl(),需要先cor()
df01 <- fortify_cor(exprSet, cor.test = TRUE, cluster.type = 'upper')#不需要cor()

g <-ggcor(df01) +
  geom_segment(aes(x = x - 0.5, y = y + 0.5, xend = x + 0.5, yend = y - 0.5),
               data = get_data(type = "diag"), size = 0.5, colour = "grey60") +
  geom_colour(data = get_data(type = "upper", show.diag = FALSE)) +
  geom_mark(data = get_data(type = "upper", show.diag = FALSE), size = 3) +
  geom_circle2(data = get_data(r >= 0.9, type = "lower", show.diag = FALSE),r = 0.8, fill = "#66C2A5") +
  geom_num(aes(num = r), data = get_data(type = "lower",show.diag = FALSE), size = 2)
g


dev.off()
#----------------------------------------ComplexHeatmap corrlation
library(ComplexHeatmap)
library(circlize)
fpkm <-read.table('strand_specific/expr/sense_fpkm.txt',sep = '\t',header = T,row.names = 1)
cor_mat <- cor(fpkm)
#od = hclust(dist(cor_mat))$order
#cor_mat = cor_mat[od, od]
nm = rownames(cor_mat)
col_fun = circlize::colorRamp2(c(-1, 0, 1), c("green", "white", "red"))
Heatmap(cor_mat, name = "correlation", col = col_fun, rect_gp = gpar(type = "none"), 
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
          if(i == j) {
            grid.text(nm[i], x = x, y = y)
          } else if(i > j) {
            grid.circle(x = x, y = y, r = abs(cor_mat[i, j])/2 * min(unit.c(width, height)), 
                        gp = gpar(fill = col_fun(cor_mat[i, j]), col = NA))
          } else {
            grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp = gpar(fontsize = 8))
          }
        }, cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE, show_column_names = FALSE)
