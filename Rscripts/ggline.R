rm(list = ls())
options(stringsAsFactors = F)

library(ggplot2)
library(reshape2)
library(dplyr)
library(plyr)
library(export)


mat <-read.csv('strand_specific/cluster10/all.DEGs.uniq.FPKM_z-score.csv')[,1:14]
table(mat$cluster)
#mat <-t(mat)
#matmelt <- melt(mat,id.vars = c('cluster','gene')) #保留id.vars，其余整合。#Cluster  gene  variable  value
#write.csv(matmelt,file = 'matrix_for_ggline.csv')#保存更改，分裂variable，为treat1.treat2。再添加一列groupgene：treat1+gene

data <-read.csv('strand_specific/cluster10/matrix_for_ggline.csv')
#data<-matmelt
#data <- data[data$cluster==1,]
table(data$time)
table(data$cluster)

#使用因子类型来固定顺序
data$time<- factor(data$time,
                   levels= c('30min','60min','24h','48h','96h','8d'),ordered = TRUE)
data$cluster<- factor(data$cluster,
                      levels= paste0("cluster",1:10),ordered = TRUE)
table(data$time)
table(data$cluster)


#--------------------------------------------------------------多个分组拟合折线图
p <- ggplot(data,aes(x=time, y=value, color=treat, group=groupgene))+ 
  theme_classic(base_size = 18) + 
  geom_line(size=0.8,alpha=0.01) +  xlab("")+ ylab("FPKM Z-score Normalized")+
  geom_hline(yintercept =0,linetype=2) + #ggtitle("Cluster") + 
  stat_summary(aes(group=treat,color=treat), #拟合的线
               fun.y= mean, geom="line", size=1.2) + 
  facet_wrap(.~cluster) +
  theme(legend.position = c(0.6,0.11),
        legend.background = element_rect(colour = 'black') )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=16,color = 'black'),
        strip.text = element_text(size = 18, 
                                  face = "bold"))
#guides(col = FALSE) #legend取消
p
#ggsave(p,filename = 'all-expr-trend.png')
graph2ppt(file='deg.cluster.pptx')

#---------------------------------------------单个分组拟合折线

ck<-data[data$treat=='CK',]
ln<-data[data$treat=='LN',]
p1<-ggplot(ln,aes(x=time, y=value, color=treat, group=groupgene)) + 
  theme_classic(base_size = 18) + 
  geom_line(aes(color=cluster),size=0.8,alpha=0.2) +   
  xlab("")+ ylab("FPKM Z-score Normalized")+
  geom_hline(yintercept =0, linetype=2 ) + #ggtitle("Cluster") + 
  stat_summary(aes(group=1),   #拟合的线
               fun.y=mean, geom="line", size=1.2, color="red") + 
  facet_wrap(cluster~.) +
  theme(legend.position = 'none',
        legend.background = element_rect(colour = 'black') )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=14,colour = 'black'),
        strip.text = element_text(size = 14, face = "bold"))

p1
#ggsave(p1,filename = 'LN-expr-trend.png')
graph2ppt(p1,file = 'LN-expr-trend.pptx')
graph2ppt(p1,file = 'CK-expr-trend.pptx')
