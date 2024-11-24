library(ggplot2)
library(ggpubr)

### revigo:
#p0 <- p1 + ggtitle("biological_process") + theme(plot.title = element_text(hjust = 0.5,lineheight=.8, face="bold"))
#p2 <- p1 + ggtitle("molecular_function") + theme(plot.title = element_text(hjust = 0.5,lineheight=.8, face="bold"))
#p3 <- p1 + ggtitle("cellular_component") + theme(plot.title = element_text(hjust = 0.5,lineheight=.8, face="bold"))


#3 plot
p <- ggarrange(p0+rremove("xlab"),p2+rremove("ylab")+rremove("xlab"),
               p3+rremove("xlab")+rremove("ylab"),align = "v", 
               common.legend = TRUE,nrow = 1, legend = "top")
figure<- annotate_figure(p,
                #top = text_grob("Data source: \n mtcars data set", color = "red", face = "bold", size = 14),
                bottom = text_grob("semantic space y",hjust = 0.5),
                #left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
                #right = "I'm done, thanks :-)!",
                fig.lab = "blue module",fig.lab.pos  = 'bottom.right',fig.lab.face = "bold",fig.lab.size =13)
figure
ggsave(figure,filename = 'blue-go.png',width = 16,height = 8) #16*6

dev.off()

save(figure,figure1, file = 'agrigo-plot-0.01.Rdata')

#2 plot
p <- ggarrange(p3+rremove("xlab"),p2+rremove("ylab")+rremove("xlab"),  
               align = "v", common.legend = TRUE,nrow = 1, legend = "top")
figure <- annotate_figure(p,
                          #top = text_grob("Data source: \n mtcars data set", color = "red", face = "bold", size = 14),
                          bottom = text_grob("semantic space y",hjust = 0.5),
                          #left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
                          #right = "I'm done, thanks :-)!",
                          fig.lab = "black module",fig.lab.pos  = 'bottom.right',fig.lab.face = "bold",fig.lab.size =14)
figure
ggsave(figure,filename = 'black-go.png',width = 12,height = 8)

dev.off()
