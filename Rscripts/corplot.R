#-----------------------------------------RNAseq样本间重复相关性分析可用FPKM,TPM(最好),也可用read counts计算（尽量不用）
rm(list = ls())
options(stringsAsFactors = F)

library(export)
#-------------相关性热图
library(pheatmap)

fpkm<- read.table('strand_specific/expr/sense_fpkm.txt',row.names = 1,header = T,
                sep = '\t')[,c(1:12)]

pheno<-data.frame(Time = factor(c(rep('0.5h',2),rep('1h',2),rep('1d',2),rep('2d',2),rep('4d',2),rep('8d',2),
                                  rep('0.5h',2),rep('1h',2),rep('1d',2),rep('2d',2),rep('4d',2),rep('8d',2))),
                  Treatment =factor(c(rep('CK',12),rep('LN',12))),
                  #rep=factor(c(rep('R1',1),rep('R2',1),rep('R1',1),rep('R2',1),rep('R1',1),rep('R2',1))),
                row.names = colnames(fpkm))
pheno<-data.frame(Time = factor(c(rep('0.5h',2),rep('1h',2),rep('1d',2),rep('2d',2),rep('4d',2),rep('8d',2))),
                  #Condition =factor(c(rep('CK',12),rep('LN',12))),
                  #rep=factor(c(rep('R1',1),rep('R2',1),rep('R1',1),rep('R2',1),rep('R1',1),rep('R2',1))),
                  row.names = colnames(fpkm))  
pheno<-data.frame(#time = factor(c(rep('30min',2),rep('60min',2),rep('24h',2),rep('48h',2),rep('96h',2),rep('8d',2))),
                  Treatment =factor(c(rep('CK',2),rep('LN',2))),
                  #rep=factor(c(rep('R1',1),rep('R2',1),rep('R1',1),rep('R2',1),rep('R1',1),rep('R2',1))),
                  row.names = colnames(fpkm)) 
#name <-colnames(fpkm)
#df <-cbind(pheno,name)
#condition =data.frame(Condition= factor(pheno$Condition),row.names = names(fpkm))
#cor_matrix <-cor(fpkm,method = 'spearman')
cor_matrix <-cor(fpkm)
colSums(fpkm)
#corrplot(cor_matrix, method = "color") 

p <- pheatmap(cor_matrix, #main='LN Sample Correlation', 
              fontsize=18,angle_col=315,
              cellwidth = 35 ,
              cluster_rows = T,cluster_cols = T,
              color=colorRampPalette(colors = c('#436eee','white','#ee0000'))(70),
              #color = colorRampPalette(colors = c('cyan1','black','yellow'))(15),
              #color =cm.colors(15),
              #color =colorRampPalette(colors = c('white','red'))(100),#[30:100]
              display_numbers = T,annotation_col = pheno,
              #border_color = 'black'
              )
p
graph2ppt(p,file='correlation-all-samples-cluster.heatmap.pptx')
graph2ppt(p,file='correlation-all-samples-without-cluster.heatmap.pptx')
graph2ppt(p,file='correlation-LN-samples-cluster.heatmap.pptx')
graph2ppt(p,file='correlation-CK-samples-cluster.heatmap.pptx')

dev.off()





library(ggcorrplot)
# 'Correlation' matrix
corr <- cor(fpkm)
# Plot
ggcorrplot(corr, hc.order = TRUE, type = "lower", 
          lab = TRUE, lab_size = 3, 
         method="circle", legend.title = 'Correlation',
        colors = c("tomato2", "white", "springgreen3"), 
       title="Correlogram of samples", ggtheme=theme_bw) 



#---------------------------------------------------------------相关性散点图
options(scipen=999) 
library(ggplot2)
library(ggpubr)

fpkm1 <-read.table('strand/antisense_fpkm.txt',header = T,sep='\t',
                  row.names = 1)[,6:29][,c(13,14,15,16,17,18)]
fpkm <- fpkm1[rowSums(fpkm1)>0,]
boxplot(log2(fpkm+1),las=2)

theme_set(theme_bw())  

head(fpkm[,3:4])
x<-fpkm$`LN.30min.R1`
y<-fpkm$`LN.30min.R2`
dat <-fpkm[,3:4]
# Scatterplot
gg <- ggplot(dat, aes(x=x, y=y)) + 
  geom_point(color="red") +stat_cor(data=dat, method = "pearson") +
  geom_smooth(method="lm", se=F) +  # stat_smooth(method="loess",se=FALSE)
  xlim(c(0, 6000)) + ylim(c(0, 6000)) + 
  labs(subtitle="LN.30min.R1 VS LN.30min.R2", y="LN.30min.R1", x="LN.30min.R2", 
       title='Correlation' )  #,  caption = "Source: midwest"
gg
