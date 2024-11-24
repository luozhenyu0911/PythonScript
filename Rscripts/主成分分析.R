rm(list=ls())

expr=read.delim('相关系数计算与PCA分析/data/exp.txt',
                header=T,row.names = 1,sep = "\t",check.names = F,na.strings = "-") #读取表达矩阵
groupInfo=read.delim('相关系数计算与PCA分析/data/group.txt',
                     header=T,sep="\t",check.names = F)# 读取分组信息

expr=na.omit(expr)                  #过滤掉包含缺失的基因
expr=expr[rowSums(expr)!=0,]        #预处理，过滤掉所有样品中均为0的基因
group=as.factor(groupInfo[match(colnames(expr),groupInfo$Sample),'Group']) #exp矩阵与分组信息匹配

pc=prcomp(t(expr),center = TRUE, scale = TRUE) #一般对表达矩阵标准化，用prcomp计算PCA结果
head(pc$x)
#=========================绘制2D PCA图=================================
if(!requireNamespace('factoextra',quietly = TRUE))  install.packages('factoextra',update=F)

library("factoextra")
plot12<-fviz_pca_ind(pc,
       axes = c(1, 2), #选择要显示的主成分
       geom.ind =  c("point", "text"), # point，text或 c("point", "text")
       pointsize = 2,  #点的大小
       col.ind=group,
#       col.ind = "black",     #点的轮廓颜色，也可使用group使其颜色与填充相同,同时分组修改形状
       palette = "lancet", #颜色方案，可选"npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty"
       #palette = c("#999999", "#E69F00"), #手动设置颜色
       addEllipses = TRUE, #添加椭圆
       repel = TRUE, #文字自动微调，防止重叠
       alpha.ind=0.5, #透明度
       invisible="quali", #如需去除重心点，设置为quali
       title="PC1 vs PC2", #主标题
       legend.title = "Group",#legend标题
       )+theme_minimal()+theme(plot.title = element_text(hjust = 0.5)) #设置主题与标题居中
plot12

plot13<-fviz_pca_ind(pc,
       axes = c(1, 3), #选择要显示的主成分
       geom.ind =  c("point", "text"), # point，text或 c("point", "text")
       pointsize = 2,  #点的大小
#      col.ind = "black",     #点的轮廓颜色，也可使用group使其颜色与填充相同
       col.ind =group,     #点的轮廓颜色，也可使用group使其颜色与填充相同
       palette = "lancet", #颜色方案，可选"npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty"
       #palette = c("#999999", "#E69F00", "#56B4E9"), #手动设置颜色
       addEllipses = TRUE, #添加椭圆
       repel = TRUE, #文字自动微调，防止重叠
       alpha.ind=0.5, #透明度
       invisible="none", #如需去除重心点，设置为quali
       title="PC1 vs PC3", #主标题
       legend.title = "Group",#legend标题
       )+theme_minimal()+theme(plot.title = element_text(hjust = 0.5)) #设置主题与标题居中
plot13

plot23<-fviz_pca_ind(pc,
       axes = c(2, 3), #选择要显示的主成分
       geom.ind =  c("point", "text"), # point，text或 c("point", "text")
       pointsize = 2,  #点的大小
       col.ind = group,     #点的轮廓颜色，也可使用group使其颜色与填充相同
#      col.ind = "black",     #点的轮廓颜色，也可使用group使其颜色与填充相同
       palette = "lancet", #颜色方案，可选"npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty"
       #palette = c("#999999", "#E69F00", "#56B4E9"), #手动设置颜色
       addEllipses = TRUE, #添加椭圆
       repel = TRUE, #文字自动微调，防止重叠
       alpha.ind=0.5, #透明度
       invisible="none", #如需去除重心点，设置为quali
       title="PC2 vs PC3", #主标题
       legend.title = "Group",#legend标题
       )+theme_minimal()+theme(plot.title = element_text(hjust = 0.5)) #设置主题与标题居中
plot23
#查看
print(plot12)
print(plot13)
print(plot23)

#保存为PDF文件
pdf("PCA.2D.pdf",w=10,h=8)
print(plot12)
print(plot13)
print(plot23)
dev.off()


#=================================两种方法绘制3D PCA图=============================
#安装两个包用于绘制3D图：
if(!requireNamespace('pca3d',quietly = TRUE))
  install.packages("pca3d",update=F)
if(!requireNamespace('rgl',quietly = TRUE))
  install.packages("rgl")

#grpcolor=c(rep('red',8),rep('blue',8))
grpcolor=rep(NA,16)
grpcolor[groupInfo$Group=="A"]="red"
grpcolor[groupInfo$Group=="B"]="yellow"


#绘制3D PCA图：
library(pca3d)
library(rgl)
pca3d(pc,              #以pc变量为输入
      group= as.factor(groupInfo$Group),  #以grpcolor定义的颜色为分组
      col=grpcolor,    #点的颜色
      radius = 4,      #修改点的大小
#      show.labels=T,   #显示样品名,如需隐藏样品名，将本句注释掉
      show.ellipses=T, #显示包裹每组样品的椭圆
      show.centroids = F,  #显示每组样品的重心点位置
      show.group.labels=T, #在重心位置显示每组样品的组名
      show.plane=T,        #显示水平面
      legend="topright"     #legend的位置，可以是bottom,top,left,right,以及bottomleft等组合
)
#弹出窗口可拖拽
#注意，运行snapshotPCA3d的时候先不要关闭弹出的窗口，否则无法保存png文件
snapshotPCA3d("PCA.3DPlot.1.png")   #保存为png文件



#另外一种方法绘制3d图
if(!requireNamespace('scatterplot3d',quietly = TRUE))
  install.packages('scatterplot3d')
library(scatterplot3d) 

grpcolor=rep(NA,length(group))
grpcolor[group=="A"]="red"
grpcolor[group=="B"]="blue"

grpshape=rep(NA,length(group))
grpshape[group=="A"]=17  #参考http://www.cookbook-r.com/Graphs/Shapes_and_line_types/
grpshape[group=="B"]=19


pdf(file='PCA.3DPlot.2.pdf',w=10,h=10)
s3d<-scatterplot3d(pc$x[,1:3], #三维坐标
             color=grpcolor, # 分组着色
             pch=grpshape,cex.symbols = 3 
#             angle=45,   #显示角度
              )
text(s3d$xyz.convert(pc$x[,1:3]),labels=colnames(expr),pos=1)
dev.off()
