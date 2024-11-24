rm(list = ls())
options(stringsAsFactors = F)

library(ggpubr)
library(reshape2)
library(ggsignif)

module<-read.table('strand/antisense_fpkm.txt',row.names = 1,header = T,
                   )[,6:29][,c(15,16,17,18)]
#module <- module1[rowSums(module1)>0,]
boxplot(log2(module+1),las=2)
module <-log2(module +1)
modulemelt<-melt(module,measure.vars = 1:ncol(module))
modulemelt$condition <-c(rep('30min',nrow(modulemelt)/2),rep('60min',nrow(modulemelt)/2))

#compare_means(value~condition, data = modulemelt, group.by = "variable")  显著性比较
colnames(module)
my_comparisons <- list(c("LN.30min.R1", "LN.30min.R2"),
                       c("LN.60min.R1", "LN.60min.R2"),c("LN.30min.R1", "LN.60min.R1"),
                       c("LN.30min.R2", "LN.60min.R2"))

g<-ggplot(modulemelt,aes(variable, value,fill=condition)) + geom_boxplot(outlier.shape = NA)+
  theme_bw()+ ylim(0,2)+
  #geom_boxplot(aes(fill=condition))+#outlier.colour="red", outlier.shape=8, outlier.size=4
  labs(x="",y='log2(FPKM+1)')+ #title="boxplot",
  #annotate("text",label="my boxplot plot",x=15,y=15,color="red",size=8)+
  #scale_fill_manual(values=c( "#E69F00", "#56B4E9"))+
  #theme(legend.position='right') + guides(fill=guide_legend(title=NULL))+
  theme(axis.title.x=element_text(color="black", size=18, face="bold"))+
  theme(axis.title.y=element_text(color="black", size=18, face="bold"))+ 
  theme(axis.text.y=element_text(size=13,color="black"))+
  theme(axis.text.x=element_text(angle=45, size=13,hjust=0.4,vjust=0.4,color="black")) + #+
  geom_signif(comparisons = my_comparisons,step_increase = 0.1,
              map_signif_level = T,test =wilcox.test )
  #stat_compare_means(comparisons = my_comparisons) + 
  #stat_compare_means(method = "anova",label.x = 3, label.y = 15) +
  #stat_compare_means(aes(label = ..p.signif..) ) +#显著性
  #xlab('')  
g
ggsave('boxplot-fpkm.png')
dev.off()
