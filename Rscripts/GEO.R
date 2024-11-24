#测序分为芯片测序和测序仪测序。芯片的DGE主要使用limma包，NGS测序使用DESeq2.
#此处主要针对GEO基因芯片数据。
#=========================================================================Download GEO 
rm(list=ls())
options(stringsAsFactors = F)

#method 1 :R包GEOquery下载
  #1 下载并保存GEO数据
    #下载有error可以：设置镜像、翻墙、rm(list=ls())、重试…
library(GEOquery)
?getGEO()
GSE_name = 'GSE115796'
options( 'download.file.method.GEOquery' = 'libcurl' ) #windows系统
gset <- getGEO( GSE_name, getGPL = F, AnnotGPL = F) #仅下载GSE，不下载和注释GPL。也可选择true。
gset
save( gset, file = 'gset.Rdata' )
  #2  加载GEO数据
    #gset包含下载的所有信息
    #由于gset是列表，故将其转为可操作的数据结构Gset对象
load("gset.Rdata")
class(gset)  # method 1 读入的是list
length(gset)
Gset <- gset[[1]] #转换为Gset对象
class(Gset)      # ExpressionSet
Gset            #Annotation: GPL6244
      
  #3 用GEOquery里的pdata函数获取样本信息,exprs()获取表达矩阵。
    #看一下pdata和exprs的结构，很明显是数据框
exprs <- exprs(Gset)  #matrix
sample <-sampleNames(Gset)
group_list <-as.character(pdata[,1]) #pdata[,1],根据分组信息所在列确定。
pdata<-pData(Gset)    #data.frame

write.table(exprs,'exprsdata.txt',sep = '\t')
write.table(pdata,'pdata.txt',sep = '\t',row.names = T)

dim(pdata)    # dim查看行列
colnames(pdata)  # colnames查看列名


#method 2 ：网页下载Series Matrix File,读入R
g  <- read.table('GSE115796_series_matrix.txt',sep='\t',quote = '',
                    fill = T,comment.char = '!',header = T) #comment.char = '!',告知注释信息为 ！，不读入。   
class(g) # method 2读入的是data.frame.
rownames(g) <- g[,1]
g <- g[,-1]


#method 3:网页下载GSE42872_RAW.tary原始数据，不建议。

#==========================================================================probe ID转换symbol
  # 1 根据ExpressionSet 中Annotation: GPL6244 下载相应探针id转换包。
    #如生信菜鸟团搜索：GPL6244 ，得到hugene10sttranscriptcluster,此是前缀，需加上.db
BiocManager::install('hugene10sttranscriptcluster.db')
library(hugene10sttranscriptcluster.db)
?hugene10sttranscriptcluster.db
ls("package:hugene10sttranscriptcluster.db")

  #若上述搜索不到或平台偏僻，可以从GEO中platforms查看GPL，使用GEOquery:
gpl <- getGEO('GPL6244', destdir = '.')
colnames(Table(gpl))
head(Table(gpl)[,c(1,12)])     #c(1,12),根据colnumbers
probe <-Table(gpl)[,c(1,12)]

   # 2 转换包中所有探针id
ids <-toTable(hugene10sttranscriptclusterSYMBOL)
length(unique(ids$symbol))
tail(sort(table(ids$symbol)))
table(sort(table(ids$symbol)))
plot(table(sort(table(ids$symbol))))

   # 3 过滤、匹配前面exprs矩阵中探针
table(rownames(exprs) %in% ids$probe_id)
dim(exprs)
exprs=exprs[rownames(exprs) %in% ids$probe_id,]
dim(exprs)
exprSet <-exprs
ids=ids[match(rownames(exprSet),ids$probe_id),]
head(ids)
exprSet[1:5,1:6]
  # 4 转换exprs的probe id
tmp = by(exprSet,ids$symbol,function(x) rownames(x)[which.max(rowMeans(x))] )
probes = as.character(tmp)
dim(exprSet)
exprSet=exprSet[rownames(exprSet) %in% probes ,]
dim(exprSet)
rownames(exprSet)=ids[match(rownames(exprSet),ids$probe_id),2]
exprSet[1:5,1:6]


#===============================================================检验exprset是否正确.
  #1 根据生物学知识,如：管家基因GAPDH、ACTB的中值应显著高于总体的基因中值。
boxplot(exprSet,las=2) #las=2,样本名竖排。
exprSet['GAPDH',]
boxplot(exprSet['GAPDH',])
exprSet['ACTB',]
boxplot(exprSet['ACTB',])
  #2 根据exprSet表达矩阵的分布图
library(reshape2)
exprSet_L=melt(exprSet)
exprSet_L
colnames(exprSet_L)=c('probe','sample','value')
group_list <-as.character(pdata[,1])
#可手动构造分组信息。需要时用stringr包,如：str_split(pdata$title,'_',simplify=T)
group_list=c(rep('RB_139',3),rep('RB_535',3),   
             rep('RB_984',3),rep('Healthy',2))
exprSet_L$group=rep(group_list,each=nrow(exprSet))
head(exprSet_L)
     # ggplot2 
library(ggplot2)
p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
print(p)
p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_violin()  #小提琴
print(p)
p=ggplot(exprSet_L,aes(value,fill=group))+geom_histogram(bins = 200)+facet_wrap(~sample, nrow = 4)
print(p)
p=ggplot(exprSet_L,aes(value,col=group))+geom_density()+facet_wrap(~sample, nrow = 4)
print(p)
p=ggplot(exprSet_L,aes(value,col=group))+geom_density() 
print(p)

p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
p=p+stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
p=p+theme_set(theme_set(theme_bw(base_size=20)))
p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank())
print(p)

   #hcluster
colnames(exprSet)=paste(group_list,1:ncol(exprSet),sep=' ')
        # Define nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), cex = 0.7, col = "blue")
hc=hclust(dist(t(exprSet)))    #dist(),计算矩阵行向量之间的距离
par(mar=c(5,5,5,10)) 
plot(as.dendrogram(hc), nodePar = nodePar,  horiz = TRUE)

   # PCA 
BiocManager::install('ggfortify')
library(ggfortify)
df=as.data.frame(t(exprSet))
df$group=group_list 
png('pca.png',res=120)
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')
dev.off()


#T.TEST (只能针对两组分组，多个分组的group_list不行)
dat <-exprSet
group_list <- as.factor(group_list)
group1 <- which(group_list==levels(group_list)[1])
group2 <- which(group_list==levels(group_list)[2])
dat1 <-dat[,group1]
dat2 <-dat[,group2]

#对单基因做T.test
t.test(exprSet[3,]~group_list) #第三个基因
#对所有基因做T.test
pvals <-apply(exprSet,1,function(x){    #str(t.test(as.numeric(x)~group_list)),为list
  t.test(as.numeric(x)~group_list)$p.value})
padj <-p.adjust(pvals,method = 'BH')
avg1<-rowMeans(dat1)
avg2<-rowMeans(dat2)
log2FC <-avg2-avg1
DEG_t.test <-cbind(avg1,avg2,log2FC,pvals,padj)


#-----------------------------ggpubr--------------------------------------
library(ggpubr)
data("ToothGrowth")
df1 <- ToothGrowth
head(df1)
p <- ggboxplot(df1, x="dose", y="len", color = "dose",
               palette = c("#00AFBB", "#E7B800", "#FC4E07"),
               add = "jitter", shape="dose") #增加了jitter点，点shape由dose映射
p

#=========================================================================芯片差异分析 limma包
# DEG by limma :3个矩阵：exprSet、design、contrast.matrix
library(limma)
design <- model.matrix(~0+factor(group_list))  #设计分组的矩阵，数据少时可手动
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
design
class(design)
     # 下面的 contrast.matrix 矩阵非常重要，制定了谁比谁这个规则
contrast.matrix <- makeContrasts(paste0(unique(group_list),
                                 collapse = "-"),levels = design)
contrast.matrix  #  -1（case）是用来比较的， 1 （control） 是用来被比较的。

     # 这个矩阵声明，我们要把case组跟control组进行差异分析比较
  ##step1
fit <- lmFit(exprSet,design)
  ##step2
fit2 <- contrasts.fit(fit, contrast.matrix) #这一步很重要
fit2 <- eBayes(fit2)  ## default no trend !!!  eBayes() with trend=TRUE
  ##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
    #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(nrDEG)


#============================================================================plot
  ## heatmap 
library(pheatmap)
choose_gene=head(rownames(nrDEG),25) #choose 25 genes
choose_matrix=exprSet[choose_gene,]
choose_matrix=t(scale(t(choose_matrix))) #scale（），normlized
pheatmap(choose_matrix)

  ## volcano plot:火山图只要求有logFC、adj.P.Val
colnames(nrDEG)
plot(nrDEG$logFC,-log10(nrDEG$adj.P.Val))

DEG=nrDEG
logFC_cutoff <- with(DEG,mean(abs( logFC)) + 2*sd(abs( logFC)) )
DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))
this_title <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',]))
this_title
head(DEG)
g = ggplot(data=DEG, aes(x=logFC, y=-log10(P.Value), color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_title  ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))  ## corresponding to the levels(res$change)
print(g)


#=====================================================================富集分析
library(clusterProfiler)
gene=head(rownames(nrDEG),1000) 
gene.df <- bitr(gene,fromType ='SYMBOL',toType = c('ENSEMBL', 'ENTREZID'),
                OrgDb = org.Hs.eg.db)
head(gene.df)
  ## KEGG pathway analysis
kk <- enrichKEGG(gene=gene.df$ENTREZID,organism= 'hsa',pvalueCutoff = 0.05)
head(kk)[,1:6]
