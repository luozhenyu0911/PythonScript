#limma的核心函数是lmFit和eBayes， 前者是用于线性拟合，后者根据前者的拟合结果进行统计推断。
#lmFit RNAseq需要:表达矩阵，分组对象, 差异比较矩阵
#与芯片数据不同，limma使用voom方法来对RNA-seq数据raw count进行normalization.
#limma包接受多种数据类型：counts（需要用voom进行normalization），
##芯片数据,log2(rpkm/FPKM+1),需要log.
#edgeR 和 DESeq2 主要接受raw counts矩阵数据。

#limma为高斯分布(正态分布)，可接受标准化后数据，如fpkm.  
#edgeR,DESeq2为负二项分布，接受read counts。

##理解boxplot的重要性，来看数据集是否需要log，以便后面才能用limma包进行差异分析


#差异比较矩阵con.matrix是否需要？两种都可以，只要区别在design矩阵写法：
#-----不需要比较矩阵con.matrix的design矩阵写法：
# design=model.matrix(~factor(sCLLex$Disease))  ##  ~ A
#fit=lmFit(sCLLex,design)
#fit=eBayes(fit)
#topTable(fit,coef=2,adjust='BH')  #coef=2

#-----需要比较矩阵con.matrix的design矩阵写法：
#design=model.matrix(~0+factor(sCLLex$Disease)) ## ~ 0 + A
#cont.matrix=makeContrasts('M-BS',levels = design)
#fit=lmFit(sCLLex,design)
#fit2=contrasts.fit(fit,cont.matrix)
#fit2=eBayes(fit2)
#topTable(fit2,adjust='BH',coef="M-BS")  #coef="M-BS"

boxplot()

rm(list = ls()) 
options(stringsAsFactors = F)
library(ggplot2)
library(limma)

countsdata<-read.table(file =  'expr-matrix.txt',
                       header = T,sep = '\t',quote = '',check.names = F)  
head(countsdata)
countsdata <-na.omit(countsdata)
exprSet<-countsdata[,7:ncol(countsdata)]
rownames(exprSet)<-countsdata[,1]
exprSet <-na.omit(exprSet)
countdata<-exprSet
countdata <- countdata[rowSums(countdata)>0,]

#时间不变
#countdata <- countdata[,c(1,2,7,8)]
countdata <- countdata[,c(3,4,9,10)]
#countdata <- countdata[,c(5,6,11,12)]

#CK恒定
#countdata <- countdata[,c(1,2,3,4)]
#countdata <- countdata[,c(1,2,5,6)]
#countdata <- countdata[,c(3,4,5,6)]

#LN恒定
#countdata <- countdata[,c(7,8,9,10)]
#countdata <- countdata[,c(7,8,11,12)]
#countdata <- countdata[,c(9,10,11,12)]
head(countdata)

countdata <-na.omit(countdata)

#第一步：构建分组矩阵
sample=names(countdata)
group=factor(c(rep('CK',ncol(countdata)/2),rep('LN',ncol(countdata)/2)),levels=c('CK','LN'))#control在前
phno=data.frame(sample,group)

design <- model.matrix(~0 + factor(phno$group))
colnames(design)=levels(factor(phno$group))
rownames(design)=colnames(countdata)
design

cont.matrix=makeContrasts(contrasts=c('LN-CK'),levels = design)  # contrasts:trt-untrt
cont.matrix

#===========================================================================================针对count

#样本间测序深度相近，最大与最小的不超过3倍，推荐该方法
#library(edgeR)
#dge <- DGEList(counts=countdata)
#dge <- calcNormFactors(dge)
#logCPM <- cpm(dge, log=TRUE, prior.count=3)
#fit <- lmFit(logCPM, design)
#fit2=contrasts.fit(fit,cont.matrix)
#fit2=eBayes(fit2)


#第二步：根据分组信息和raw counts表达矩阵进行normalization，样本间测序深度差别较大时，推荐该方法
#voom函数将counts（所有计数加0.5以避免对数取零）转换为logCPM值，之后将logCPM值矩阵进行标准化。
v <- voom(countdata,design,plot=TRUE, normalize="quantile")

#---------------------------------------------------------针对voom，无需运行，比较NORM前后的数据分布。
png("limma_voom_RAW_vs_NORM.png",height = 1000,width = 1400)
exprSet_new=v$E
par(cex = 0.7)
n.sample=ncol(countdata)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
par(mfrow=c(2,2))
boxplot(countdata, col = cols,main="expression value",las=2)
boxplot(exprSet_new, col = cols,main="expression value",las=2)
hist(as.matrix(countdata))
hist(exprSet_new)
dev.off()

#第三步：做差异分析，提取差异分析结果
fit <- lmFit(v, design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)

tempOutput = topTable(fit2, coef='LN-CK', n=Inf)
DEG_limma_voom = na.omit(tempOutput)
head(DEG_limma_voom)
write.csv(DEG_limma_voom,"DEG_limma_voom-LNvsCK-24h.csv",quote = F) #自行筛选DEG

DEgenes <- topTable(fit2, coef="LN-CK", number=Inf,p.value=0.05, adjust="BH", lfc=1)
res <- decideTests(fit2, p.value=0.05, lfc=1)
summary(res)

resdata<-merge(as.data.frame(DEG_limma_voom), as.data.frame(countdata), by='row.names',sort=F)                     
head(resdata)
resdata$significant <- 'unchanged'  
resdata$significant[resdata$adj.P.Val <= 0.05 & resdata$logFC>=1] <- 'upregulated'     
resdata$significant[resdata$adj.P.Val <= 0.05 & resdata$logFC<=-1] <- 'downregulated'  
head(resdata)
table(resdata$significant)

DEG=as.data.frame(DEG_limma_voom)
DEG=na.omit(DEG)  
nrDEG=DEG   
nrDEG$change = as.factor(ifelse(nrDEG$adj.P.Val < 0.05 & abs(nrDEG$logFC) >1, 
                                ifelse(nrDEG$logFC > 1 ,'Up in LN','Up in CK'),'NotSig')) 
this_tile <- paste0("LN vs CK", 
                    '\nup gene : ',nrow(nrDEG[nrDEG$change =='Up',]) , 
                    '\ndown gene : ',nrow(nrDEG[nrDEG$change =='Down',])) 
#this_tile <- 'blue'
valcano <- ggplot(data=nrDEG, aes(x=logFC, y=-log10(adj.P.Val), color=change)) + 
  geom_point(alpha=0.4, size=4) + 
  theme_set(theme_set(theme_bw(base_size=14)))+ 
  xlab("log2FoldChange") + ylab("-log10(padj)") + 
  ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+  
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
  geom_hline(yintercept=-log10(0.05), lty=2, col="gray", lwd=0.5)+
  scale_colour_manual(values = c('red','blue','green'), limits=c("Up", "Down", "NoDiff")) 
valcano

summary(res)
table(resdata$significant)
ggsave(valcano,filename = 'volcano-limma-LNvsCK-24h.png') 
write.csv(resdata,'resdata-limma-LNvsCK-24h.csv')

#save(DEG_limma_voom,group_list,exprSet,file='DEG_result.Rdata')  #保存工作环境
#load(file='DEG_result.Rdata')   #加载环境
