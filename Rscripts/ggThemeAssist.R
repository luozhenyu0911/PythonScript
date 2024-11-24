p = ggplot(data, aes(x=variable, y = value, fill = phylumpro )) +
  geom_bar(stat = "identity",position="fill", width=1)+
  scale_y_continuous(labels = scales::percent) +
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~ group, scales = "free_x", switch = "x") +  main_theme +
  # 关闭x轴刻度和标签
  theme(axis.ticks.x = element_blank(), legend.position="top", 
        axis.text.x = element_blank(), strip.background = element_blank())+
  xlab("Groups")+ylab("Percentage (%)")




#-------------Main Theme
main_theme = theme(panel.background=element_blank(),
                     panel.grid=element_blank(),
                     axis.line.x=element_line(size=1, colour="black"),
                     axis.line.y=element_line(size=1, colour="black"),
                     axis.ticks=element_line(color="black"),
                     axis.text=element_text(color="black", size=16),
                     legend.position=c(0.9,0.8), #"right"
                     legend.background=element_blank(),
                     legend.key=element_blank(),
                     legend.text= element_text(size=14),
                     text=element_text(size=16),
                     axis.title = element_text(size =18,face ='bold'))

g <- ggplot(mtcars, aes(x = hp, y = mpg, colour = as.factor(cyl))) + 
  geom_point(size=5) + main_theme
g


#--------------------------------------------------------------------
library(ggplot2)
library(ggThemeAssist)
# 使用mtcars生成一个点图示例
gg <- ggplot(mtcars, aes(x = hp, y = mpg, colour = as.factor(cyl))) + geom_point(size=5)
gg
# 开始调整主题
ggThemeAssistGadget(gg)
#ggThemeAssist生成的主题代码：
gg <- gg + theme(plot.caption = element_text(vjust  =  1)) + 
  theme(axis.line = element_line(size  =  0.8, linetype  =  'solid')) + 
  theme(axis.ticks = element_line(size  =  0.8)) + 
  theme(axis.title = element_text(size  =  17, face  =  'bold')) + 
  theme(axis.text = element_text(size  =  16, colour  =  'black')) + 
  theme(axis.text.x = element_text(size  =  16)) + 
  theme(axis.text.y = element_text(size  =  16)) + 
  theme(legend.text = element_text(size  =  12)) + 
  theme(legend.title = element_text(size  =  14)) + 
  theme(panel.background = element_rect(fill  =  NA)) + 
  theme(legend.key = element_rect(fill  =  NA)) + 
  theme(legend.background = element_rect(fill  =  NA)) + 
  theme(legend.position = c(0.94, 0.8)) + 
  labs(x  =  NULL, y  =  'Num', colour  =  'legend')
gg





library(ggplot2)
p <- ggplot(mtcars,aes(x = wt, y = mpg)) +geom_point()
ggannotate::ggannotate(p)




devtools::install_github("dreamRs/esquisse")
