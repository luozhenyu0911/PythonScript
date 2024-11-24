#GSVA：Gene Set Variation Analysis，被称为基因集变异分析。用来评估芯片和转录组的基因集富集结果。
#主要是通过将基因在不同样品间的表达量矩阵转化成基因集在样品间表达量ES矩阵，来评估不同的sets或通路在不同样品间是否富集。
#其实就是研究这些感兴趣的基因集在不同样品间的差异，或者寻找比较重要的基因集。
#它的主要输入文件为基因表达量矩阵和基因集文件，通过gsva的方法就可以得出结果。
#GSVA将表达矩阵转换成通路富集分数(ES)矩阵，再借用limma包的 lmFit 分析得到差异通路。
#计算每个样本在基因集中的富集程度，然后计算Fold change 和P value。

#对于RNA-seq数据，如果是read count可以选择Possion泊松分布;
#如果是均一化后的表达值log2(cpm/fpkm/tpm+1)，选择默认参数高斯正态分布就可.

#GSVA和GSEA的不同是:GSEA是先根据分组对基因进行差异分析,再将基因富集到gene set中.
#GSVA是先将基因对应到gene set上,得到不同的gene set的得分,然后再根据分组,对gene set的得分进行差异分析.

#gsva()重要的参数是： mx.diff
#default argument mx.diff=TRUE 来获得近似正态分布的ES（二项分布：read counts ??），
#如果设置为false，那么通常是每个基因的GSVA富集得分ES的双峰分布（高斯分布:fpkm ??）.

#差异gene sets 筛选条件一般仅用padj=0.001 or 0.01 or others,即-log2(padj) ,一般不用FC筛选。

rm(list = ls())
options(stringsAsFactors = F)


#browseVignettes("GSVA")
if(T){
  #library(GSVAdata)
  library(genefilter)
  library(RColorBrewer)
  library(Biobase)
  library(GSVA)
  library(GSEABase)
  library(limma)
  library(pheatmap)
  library(ggplot2)
  library(ggrepel)
  
#读取基因集文件
#gmt格式，每个基因集是一行，第一列是基因集的ID，第二列是名字，后面的列就是基因集所含有的基因。
geneSets <- getGmt("genesets-log.txt.gmt")
class(geneSets) # GeneSetCollection
names(geneSets)
#geneSets[[1]]

#读取表达量文件并去除重复
mydata <- read.csv('combat_data_batch5.csv',row.names = 1)
table(duplicated(rownames(mydata)))
exp=mydata[!duplicated(rownames(mydata)),]

#exp <- exp[,c(1,2,9,10)]
#exp <- exp[,c(3,4,11,12)]
#exp <- exp[,c(5,6,13,14)]
#exp <- exp[,c(7,8,15,16)]


#将数据框转换成矩阵,并且CPM/FPKM需要log().
mydata= as.matrix(log2(exp+1))

#使用gsva分析, count:kcdf="Poisson"，默认mx.diff=TRUE, parallel.sz 线程数
#min.sz=1,max.zs=Inf，设置geneset包含基因的最小和最大值。max.sz小于最大数目的基因集时，gsva()会舍弃该set.

#mx.diff=TURE, es值是一个近似的正态分布
es.dif <- gsva(mydata, geneSets, mx.diff=TRUE, verbose=FALSE, parallel.sz=4)

p <-pheatmap(es.dif,angle_col = 315,fontsize=16,cluster_cols = F,cluster_rows = F,
         color=colorRampPalette(colors = c('#436eee','white','#ee0000'))(100))#annotation_col=df
dev.off()

#mx.diff=FALSE, es值是一个双峰的分布
es.max <- gsva(mydata, geneSets, mx.diff=FALSE, verbose=FALSE, parallel.sz=4)

pp <- pheatmap(es.max,angle_col = 315,fontsize=16,cluster_cols = F,cluster_rows = F,
         color =colorRampPalette(colors = c('#436eee','white','#ee0000'))(100))#annotation_col=df
dev.off()

#看一下两种不同分布的效果,es.max 是高斯正态分布，es.dif 是二项分布.
png('density-group4.png')
par(mfrow=c(1,2), mar=c(4, 4, 4, 1))
plot(density(as.vector(es.dif)), main="Difference between largest\npositive and negative deviations",
     xlab="GSVA score", lwd=2, las=1, xaxt="n", xlim=c(-0.75, 0.75), cex.axis=0.8)
axis(1, at=seq(-0.75, 0.75, by=0.25), labels=seq(-0.75, 0.75, by=0.25), cex.axis=0.8)
plot(density(as.vector(es.max)), main="Maximum deviation from zero",xlab="GSVA score", 
     lwd=2, las=1, xaxt="n", xlim=c(-0.75, 0.75), cex.axis=0.8)
axis(1, at=seq(-0.75, 0.75, by=0.25), labels=seq(-0.75, 0.75, by=0.25), cex.axis=0.8)
dev.off()

#limma设置分组的方法
#grouplist <- factor(c(rep('BS',ncol(exp)/2),rep('M',ncol(exp)/2)))
#design <- model.matrix(~0+grouplist)   #design <- model.matrix(~group_list) 
#colnames(design) <- levels(grouplist)
#rownames(design) <- colnames(mydata)

#设置分组
sample=names(exp)
group=factor(c(rep('BS',ncol(exp)/2),rep('M',ncol(exp)/2)),levels=c('BS','M'))
phno=data.frame(sample,group)

design <- model.matrix(~0 + factor(phno$group))
colnames(design)=levels(factor(phno$group))
rownames(design)=colnames(exp)
design

cont.matrix=makeContrasts(contrasts=c('M-BS'),levels = design)  # contrasts:trt-untrt
cont.matrix

#获取需要进行差异分析的分组
res=es.max  #mix.diff= F
#定义阈值
#logFCcutoff <- 1   #一般仅用于gene水平，而不用于sets或pathway
adjPvalueCutoff <- 0.01  #可调节
#进行差异分析
fit <- lmFit(res, design)
fit <- contrasts.fit(fit, cont.matrix) ##这一步很重要，自行看看效果
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef="M-BS", number=Inf)  # coef=treat     
#allGeneSets1 <- topTable(fit,  number=Inf)
#df=allGeneSets[allGeneSets$adj.P.Val<adjPvalueCutoff & abs(allGeneSets$logFC) > logFCcutoff,]
DEgeneSets <- topTable(fit, coef="M-BS", number=Inf, p.value=adjPvalueCutoff, adjust="BH")
#DEgeneSets1 <- topTable(fit, number=Inf, p.value=adjPvalueCutoff, adjust="BH", lfc=logFCcutoff)

res <- decideTests(fit, p.value=adjPvalueCutoff)
summary(res)
df=allGeneSets[allGeneSets$adj.P.Val< adjPvalueCutoff ,]
#df1=allGeneSets1[allGeneSets$adj.P.Val<0.001 & abs(allGeneSets$logFC) > 0.5,]
df
#allGeneSets0 <- topTable(fit, coef="M", number=Inf) 
#DEgeneSets0 <- topTable(fit, coef="M", number=Inf, p.value=adjPvalueCutoff, adjust="BH")

#画火山图
# Hide some of the labels, but repel from all data points
allGeneSets$label <- rownames(allGeneSets)
allGeneSets$label[(nrow(df)+1): nrow(allGeneSets)] <- ""

allGeneSets$change = as.factor(ifelse(allGeneSets$adj.P.Val <adjPvalueCutoff & abs(allGeneSets$logFC) >0, 
                                ifelse(allGeneSets$logFC > 0 ,'Up','Down'),'NoDiff')) 
this_tile <- paste0("M vs BS", 
                    '\nup genesets : ',nrow(allGeneSets[allGeneSets$change =='Up',]) , 
                    '\ndown genesets : ',nrow(allGeneSets[allGeneSets$change =='Down',])) 
g <- ggplot(allGeneSets,aes(logFC,-1*log10(adj.P.Val),label = rownames(allGeneSets)))+
  geom_point(size=4,aes(color = change)) + 
  #xlim(-1,1) + ylim(0,9)+
  theme(plot.title = element_text(size=15,hjust = 0.5))+
  labs(title=this_tile,x="log2(FC)", y="-log10(FDR)")+
  scale_color_manual(values =c("#00ba38","#619cff","#f8766d"))+
  geom_hline(yintercept=-log10(adjPvalueCutoff ),linetype=4)+
  #geom_vline(xintercept=c(-0.5,0.5),linetype=4)+
  geom_text_repel(aes(logFC,-1*log10(adj.P.Val),label = label), #fill=factor(allGeneSets$label)),
                  size = 5,box.padding = unit(0.5, "lines"),
                  show.legend = F,segment.color = "black",point.padding = unit(0.8, "lines"))+ 
  theme_classic()   

#png('pheatmap-group1-diffsets-mx.diff=F-M_VS_BS.png')
#pheatmap(degsetsp,angle_col = 315,cluster_cols = F,cluster_rows = F,
 #        color =colorRampPalette(colors = c('#436eee','white','#ee0000'))(100))#差异基因集
#dev.off()

#--------------------------------------------------------------------------------mix.diff=T
res1=es.dif

fit1 <- lmFit(res1, design)
fit1 <- contrasts.fit(fit1, cont.matrix) 
fit1 <- eBayes(fit1)
allGeneSets1 <- topTable(fit1, coef="M-BS", number=Inf) 
DEgeneSets1 <- topTable(fit1, coef="M-BS", number=Inf, p.value=adjPvalueCutoff, adjust="BH")
res1 <- decideTests(fit1, p.value=adjPvalueCutoff)
summary(res1)
df1=allGeneSets1[allGeneSets1$adj.P.Val<adjPvalueCutoff ,]
df1

allGeneSets1$label <- rownames(allGeneSets1)
allGeneSets1$label[(nrow(df1)+1):nrow(allGeneSets1)] <- ""  #调节
   
allGeneSets1$change = as.factor(ifelse(allGeneSets1$adj.P.Val <adjPvalueCutoff & abs(allGeneSets1$logFC) >0, 
                                      ifelse(allGeneSets1$logFC > 0 ,'Up','Down'),'NoDiff')) 
this_tile <- paste0("M vs BS", 
                    '\nup genesets : ',nrow(allGeneSets1[allGeneSets1$change =='Up',]) , 
                    '\ndown genesets : ',nrow(allGeneSets1[allGeneSets1$change =='Down',])) 

gg  <- ggplot(allGeneSets1,aes(logFC,-1*log10(adj.P.Val),label = rownames(allGeneSets1)))+
  geom_point(size=4,aes(color = change )) + 
  #xlim(-1,1) + ylim(0,10)+
  theme(plot.title = element_text(size=15,hjust = 0.5))+
  labs(title=this_tile,x="log2(FC)", y="-log10(FDR)")+
  scale_color_manual(values =c("#00ba38","#619cff","#f8766d"))+
  geom_hline(yintercept=-log10(adjPvalueCutoff ),linetype=4)+
  #geom_vline(xintercept=c(-0.5,0.5),linetype=4)+
  geom_text_repel(aes(logFC,-1*log10(adj.P.Val),label = label), 
                  size = 5,box.padding = unit(0.5, "lines"),#fill=
                  show.legend = F,segment.color = "black",point.padding = unit(0.8, "lines"))+ 
  theme_classic()   
dev.off()
#DEgeneSetspkm1 = merge(DEgeneSets1,es.dif,by=0,all.x=TRUE)[,-c(2:7)]
#degsetsp1=DEgeneSetspkm1[,-1]
#row.names(degsetsp1)=DEgeneSetspkm1[,1]
}

ggsave(g ,filename = 'volcano-module-group-mx.diff=F-0.01fdr.png') 
ggsave(gg,filename = 'volcano-module-group-mx.diff=T-0.01fdr.png') 
#summary(res)
#summary(res1)
#save(g,gg,allGeneSets1,allGeneSets,file = 'gsva-plot-1-fdr0.01.Rdata')
#save(g,gg,allGeneSets1,allGeneSets,file = 'gsva-plot-2-fdr0.01.Rdata')
#save(g,gg,allGeneSets1,allGeneSets,file = 'gsva-plot-3-fdr0.01.Rdata')
#save(g,gg,allGeneSets1,allGeneSets,file = 'gsva-plot-4-fdr0.01.Rdata')
save(g,gg,allGeneSets1,allGeneSets,file = 'gsva-plot-fdr0.01.Rdata')
