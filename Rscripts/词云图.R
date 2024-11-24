if(!requireNamespace('devtools',quietly = TRUE))
  install.packages('devtools',update=F)
if(!requireNamespace('pheatmap',quietly = TRUE))
  devtools::install_github("lchiffon/wordcloud2")

library(wordcloud2)


rm(list=ls())
setwd('C:/Users/ZQQ/Desktop/植物生信SHOP网课资料/cluster')


#读取数据包含两列，一列geneID，一列数值，数值可以为FC，也可以为p值等用来展示
data <- read.delim("data/wordCloud.data.txt",header=T,)

#绘制词云图并赋值给word变量
word <- wordcloud2(data,shape='circle',   
            #形状,circle,cardioid,diamond,triangle-forward', 'triangle', 'pentagon', and 'star'.
                color = "random-light", #'random-dark' and 'random-light' 
                backgroundColor = "white", #背景颜色
                size=1, #大小
                )

word   #把图片打印到屏幕，可以直接export保存。


#若要保存为图片，需要先用htmlwidgets报错为HTML，再借助weshot包将html转换为png
if(!requireNamespace('webshot',quietly = TRUE))
  install.packages('webshot',update=F)
if(!requireNamespace('htmlwidgets',quietly = TRUE))
  install.packages('htmlwidgets',update=F)
webshot::install_phantomjs()#安装


library(webshot)
library("htmlwidgets")

saveWidget(word,"wordcloud.html",selfcontained = F) #保存为html
webshot("wordcloud.html","wordcloud.png", delay =10, #delay保存所需时间，长点好
        vwidth = 1500, vheight=1000) #转换为png
