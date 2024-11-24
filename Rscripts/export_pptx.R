rm(list = ls())
options(stringsAsFactors = F)


if(T){
  x = EGG
  df = data.frame(x)
  ## 计算富集分数
  x@result$richFactor =x@result$Count / as.numeric(sub("/\\d+", "", x@result$BgRatio))
  y =x@result
  library(dplyr)
  library(ggplot2)
  showCategory = 20
  y %>% 
    arrange(p.adjust) %>% 
    slice(1:showCategory) %>% 
    ggplot(aes(richFactor,forcats::fct_reorder(Description, richFactor))) + 
    geom_segment(aes(xend=0, yend = Description)) +
    geom_point(aes(color=p.adjust, size = Count)) +
    ## 调整颜色的区间,begin越大，整体颜色越明艳
    scale_color_viridis_c(begin = 0.3, end = 1) +
    ## 调整泡泡的大小
    scale_size_continuous(range=c(2, 10)) +
    theme_bw() + 
    xlab("rich factor") +
    ylab(NULL) + 
    ggtitle("")
}

library(export)

#演示图形
library(effects)
fit=lm(prestige ~ type + income*education, data=Prestige)
plot(Effect(c("income", "education"), fit),multiline=T, ci.style="bands")
#图片先生成出来。
graph2ppt(file="test.pptx", width=7, height=5)  #导出为PPT ，打开后:排列 -> 取消组合，那么所有元素可以自由编辑.
graph2eps(file="dotplot2.eps") #导出为AI格式




#-------------------------------------eoffice
library(ggplotify)
library(eoffice)

f = "eoffice.pptx"
p = as.ggplot(~plot(cars, cex.lab=2, cex.main=2,
                    xlab="biobabble", ylab="biobabble",main = "演示专用"))
p

eoffice::topptx(p, f)

library(rvcheck)
o(f)

#----------------------open file in R by calling software in windows
library(rvcheck)
o(file = 'SML-毕业论文-预定稿-批注.doc')

