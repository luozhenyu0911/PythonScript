rm(list = ls())
options(stringsAsFactors = F)

library(contamDE)
library(ggfortify)
library(ggplot2)

if(T){
countsdata<-read.table(file = paste0(getwd(), '/gene_count-stringtie.txt'),
                       header = T,sep = '\t',quote = '',check.names = F,row.names = 1)  
countsdata <-na.omit(countsdata)
exprSet<-countsdata

}
table(is.na(exprSet))

#时间不变
#exprSet <- exprSet[,c(1,2,7,8)]
#exprSet <- exprSet[,c(3,4,9,10)]
exprSet <- exprSet[,c(5,6,11,12)]

if(T){
exprSet <-na.omit(exprSet)

sample=names(exprSet)
Condition=factor(c(rep('CK',ncol(exprSet)/2),rep('LN',ncol(exprSet)/2)))
pheno=data.frame(row.names = sample, Condition)  

p <- autoplot(prcomp(t(exprSet)), data=pheno, colour='Condition',size=4, #ellipse.border.remove=T,ellipse.type="convex",
              shape='Condition', frame=F) +
  #geom_text(aes(label=rownames(data), vjust=1, colour = Source),size=4)+
  #geom_vline(xintercept = 0, colour="#990000",linetype="dashed") +
  #geom_hline(yintercept = 0, colour="#990000",linetype="dashed") +
  theme_bw()
p 

exprSet <- exprSet[rowSums(exprSet)>1,]
d <- contamDE(exprSet,R=2,match=TRUE)   #默认countdata中在前列为control，在后列为treat
names(d)

res=as.data.frame(d$LR)
res$change = as.factor(ifelse(res$p.value < 0.05 & abs(res$logFC) >1, 
                              ifelse(res$logFC > 1 ,'Up','Down'),'NoDiff')) 

DEG=as.data.frame(res)
DEG=na.omit(DEG)  
nrDEG=DEG     

nrDEG$change = as.factor(ifelse(nrDEG$p.value < 0.05 & abs(nrDEG$logFC) >1, 
                                ifelse(nrDEG$logFC > 1 ,'Up','Down'),'NotSig')) 
this_tile <- paste0("24H-LNvsCK", 
                    '\nup gene : ',nrow(nrDEG[nrDEG$change =='Up',]) , 
                    '\ndown gene : ',nrow(nrDEG[nrDEG$change =='Down',])) 
library(ggplot2)
valcano <- ggplot(data=nrDEG, aes(x=logFC, y=-log10(p.value), color=change)) + 
  geom_point(alpha=0.4, size=4) + 
  xlim(c(-12,12))+ylim(c(0,20))+
  theme_set(theme_set(theme_bw(base_size=20)))+ 
  xlab("log2FoldChange") + ylab("-log10(padj)") + 
  ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+  
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
  geom_hline(yintercept=-log10(0.05), lty=2, col="gray", lwd=0.5)+
  scale_colour_manual(values = c('red','blue','green'), limits=c("Up", "Down", "NotSig")) 
valcano
}

d$W  # treat组的预估的数据的纯度的占比。
table(res$change)

ggsave(p,filename = 'pca-contamDE-all.png')
ggsave(valcano,filename = 'volcano-contamDE-all.png') 
write.csv(DEG,'resdata-contamDE-all.csv')
