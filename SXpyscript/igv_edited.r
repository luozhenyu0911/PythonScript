#!/usr/bin/env Rscript
library(ggplot2)

args=commandArgs(T)
a <- read.table(args[1],head=T,sep="\t")
plot <- ggplot(data=a,aes(x=position,y=value,fill=tag))+
        geom_area()+
		guides(fill=FALSE)+
        theme(panel.background=element_rect(fill='transparent', color='white'),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.ticks.x=element_blank(),
			  strip.text.y=element_text(angle=180,size = 10),
			  strip.background=element_blank()) + facet_grid(tag~.,switch = "y")
ggsave(args[2],plot,width = 8, height = 1)
