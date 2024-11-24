library(DESeq2)
countsdata<-read.table(file = paste0(getwd(), '/all.id.txt'),header = T,sep = '\t',quote = '',check.names = F)  
head(countsdata)
countdata<-countsdata[,7:ncol(countsdata)]
rownames(countdata)<-countsdata[,1]
write.csv(countdata, file = "countdata2.csv", row.names = F, quote = F)
names(countdata)<-c("1D4","2D4","3D4","1L4","2L4","3L4","1D12","2D12","3D12","1L12","2L12","3L12")
group_list=c("D4","D4","D4","L4","L4","L4")
exprSet<-countdata[,1:6]
colData<-data.frame(row.names = colnames(exprSet),
                    group_list=group_list)
dim(colData)
dds<-DESeqDataSetFromMatrix(countData = exprSet, 
                            colData = colData,  
                            design = ~group_list)     
dds  
dds <- dds[rowSums(counts(dds))> 1,]    
dds<-DESeq(dds) 
resultsNames(dds)
res #一开始是一些描述性信息，跟矩阵是不一样的
res <- results(dds, contrast=c("group_list","D4","L4")) #返回差异分析结果
resOrdered <- res[order(res$padj),] 
head(resOrdered) 
DEG=as.data.frame(resOrdered) 
#..............................差异分析到这里就可以了，差异基因也有了

#画个热图
DEG=na.omit(DEG)  #去除一些NA的值
nrDEG=DEG       #要改一下=后后边的数据
## heatmap 热图
library(pheatmap) 
choose_gene=head(rownames(nrDEG),50) ## 50 maybe better 
choose_matrix=exprSet[choose_gene,] 
choose_matrix=t(scale(t(choose_matrix))) 
pheatmap(choose_matrix,filename = 'DEG_top50_heatmap.png') 

## volcano plot 
colnames(nrDEG) 
#plot(nrDEG$logFC,-log10(nrDEG$pvalue)) 


DEG=nrDEG 



#logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) ) 
logFC_cutoff=1 


DEG$change = as.factor(ifelse(DEG$padj < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff, 
                              ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT') 
) 
this_tile <- paste0('Cutoff for log2FoldChange is ',round(logFC_cutoff,3), 
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) , 
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',]) 
) 

library(ggplot2)
g = ggplot(data=DEG,  
           aes(x=log2FoldChange, y=-log10(padj),  
               color=change)) + 
  geom_point(alpha=0.4, size=1.75) + 
  theme_set(theme_set(theme_bw(base_size=20)))+ 
  xlab("log2 fold change") + ylab("-log10 padj") + 
  ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+ 
  scale_colour_manual(values = c('blue','black','red')) ## corresponding to the levels(res$change) 
print(g) 
ggsave(g,filename = 'volcano.png') 