#!/usr/bin/env Rscript
library(ggplot2)

args=commandArgs(T)
a <- read.table(args[1],head=T,sep="\t")
plot <- ggplot(data=a,aes(x=region,y=value))+
        geom_area(aes(y = value),fill = "blue")+
        theme(panel.background=element_rect(fill='transparent', color='white'),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.ticks.x=element_blank())
ggsave(args[2],plot,width = 8, height = 0.5)
