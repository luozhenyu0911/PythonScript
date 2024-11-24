countsdata<-read.table(file = paste0(getwd(), '/all.id.txt'),header = T,sep = '\t',quote = '',check.names = F)  
Edata<-clusterGO
colnames(down123)
colnames(Edata)[4]='Count'
colnames(Edata)[5]='Generatio'
colnames(Edata)[6]='pvalue'
colnames(Edata)
head(Edata)
library(ggplot2)
p2 <- ggplot(Edata, aes(x=v2, y=v1)) + 
  geom_point(aes( size= Count , colour = -log10(pvalue))  ) + 
   # scale_y_discrete(limits=Edata$Description)+
  ggtitle("GO enrichment")  +  
  scale_color_gradient(low = 'blue', high = 'red') + 
   # xlim(range(Edata$type)) +
  theme(axis.text.x=element_text(angle=45,size=8, vjust=1,hjust=1 ,color='black'),
        axis.text.y=element_text(angle=0,size=10, vjust=0.3,color='black'),
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =11),
        panel.background = element_rect(fill="transparent"),
        panel.border =element_rect(fill='transparent', color='black'))
p=p2+labs(x = "",
          y = "")

print(p)




p1 <- ggplot(data=Edata)+  
  geom_bar(aes(x=Description,y=-log10(pvalue), fill=Ontology), stat='identity') +
  coord_flip() + scale_x_discrete(limits=Edata$Description) +
  theme(axis.text.x=element_text(angle=0,size=8, vjust=0.7),
        axis.text.y=element_text(angle=0,size=10, vjust=0.7),
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =14),
        panel.background = element_rect(fill="transparent"),
        panel.border =element_rect(fill='transparent', color='black'))

p1<-p1+labs(x = "",
             y = "-log10(pvalue)")
p1

library(pheatmap)
library(reshape2)
data<-d4d12ngo
data1<-acast(data,v1~v2)
data2<--log10(data1)
data2[is.infinite(data2)]<-0
pheatmap(data2,legend_labels="",show_rownames = T,
         cluster_row = FALSE,cluster_col = FALSE,
         main="",angle_col=0,
         border_color = "black",
         legend = T,
         color = colorRampPalette(colors = c("white","blue"))(100))
#ggplot2画热图富集


datas<-datas[2:12,]
library(ggplot2)
datas<-clusterGO
colnames(clusterGO)
datas$v3<--log10(datas$v3)
datas$v3[is.infinite(datas$v3)]<-0
ggheatmap <- ggplot(datas, aes(v2, v1, fill = v3))+
  geom_tile(color = "black")+
  scale_fill_gradient2( high = "red", mid = "yellow", 
                       midpoint = 0,  space = "Lab", 
                       name="-log10(pvalue)") +
  # facet_grid(cols = vars(v4)) +theme(axis.text.x =element_blank())+
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, face = "bold",
                                   size = 10, hjust = 1),
        axis.text.y = element_text(face = "bold"))+
  coord_fixed()+panel_cols()
# Print the heatmap
print(ggheatmap)
ggheatmap + 
  # geom_text(aes(v2, v1, label = v3), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())
    # legend.justification = c(1, 0),
    # legend.position = c(1, 0.9),
    # legend.direction = "horizontal")
# +
#   guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
#                                title.position = "top", title.hjust = 0.5))
