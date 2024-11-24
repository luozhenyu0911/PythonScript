#WGCNA旨在寻找协同表达的基因模块(module)，并探索基因网络与关注的表型间的关联关系，以及网络中的核心基因。
#WGCNA适用于复杂的数据模式，推荐5组(或者15个样品)以上的数据.
#WGCNA分为表达量聚类分析和表型关联两部分.
#主要包括基因间相关系数计算、基因模块的确定、共表达网络、模块与性状关联四个步骤。
#一般来说，TOM就是WGCNA分析的最终结果，后续的只是对TOM的下游注释。
#下游分析：得到模块之后的分析有：1.模块的功能富集  2.模块与性状之间的相关性 3.模块与样本间的相关系数
           #挖掘模块的关键信息：1.找到模块的核心基因  2.利用关系预测基因功能

#基因表达矩阵: 如是芯片数据，常规的归一化矩阵即可，如是转录组数据，最好是RPKM/TPM/CPM值或者其它归一化好的表达量。
               #常规表达矩阵即可，即基因在行，样品在列，进入分析前做一个转置。
               #CPM、FPKM或其它标准化方法影响不大，推荐使用Deseq2的varianceStabilizingTransformation
               #或log2(fpkm+1)对标准化后的数据做个转换。如果数据来自不同的批次，需要先移除批次效应。
               #如果数据存在系统偏移(boxplot)，需要做下quantile normalization。
#性状与分组矩阵：用于关联分析的性状必须是连续性数值型特征 (如Height, Weight,Diameter)。
                #如果是样品分组信息或分类变量，需要转换为0-1矩阵的形式
                #1表示属于此组或有此属性，0表示不属于此组或无此属性。

#无向网络在power小于15或有向网络power小于30内，
#没有一个power值可以使无标度网络图谱结构R^2达到0.8或平均连接度降到100以下，
#可能是由于部分样品与其他样品差别太大造成的。
#这可能由批次效应、样品异质性或实验条件对表达影响太大等造成, 
#可以通过绘制样品聚类查看分组信息、关联批次信息、处理信息和有无异常样品 。
#如果这确实是由有意义的生物变化引起的，也可以使用后面程序中的经验power值.

rm(list = ls())
options(stringsAsFactors = F)

#=================================================================================step1: 输入数据的准备
library(reshape2)
library(stringr)
library(WGCNA)
library(export)
enableWGCNAThreads(nThreads = 7)  #打开多线程 

#使用 Log(FPKM/TPM/CPM+1)  ,或者count标准化后的矩阵。不用原始counts。
#如果通过基因表达量上的差异来过滤基因，就相当于人为地去划分模块了，
#而我们要利用未经差异筛选后的表达矩阵来通过表达量高低与否将基因分在不同模块。
fpkm <- read.table('strand_specific/expr/sense_fpkm.txt',
                    row.names = 1,header = T,sep = '\t')  #允许表达量为负数。
#fpkm <-log2(fpkm+1) #或不进行log()

#table(fpkm<=0)
#fpkm <- fpkm[rowMeans(fpkm)>1,]
#那么如果不进行基因表达量上的差异来筛选基因，那么现在有几万个基因，
#而且又有那么多表达量在很多样本中都为零的基因，该如何过滤呢？
fpkm<-fpkm[!apply(fpkm,1,function(x){sum(floor(x)==0)> 20 }),]  #floor(),向下取整
#过滤基因，一个基因在所有24样本中如果有超过20个样本表达量为0，那这基因就不要。90%？？

#floor(nrow(fpkm)*0.75)

## 因为WGCNA针对的是基因进行聚类，而一般我们的聚类是针对样本用hclust即可，所以这个时候需要转置。
WGCNA_matrix = t(fpkm[order(apply(fpkm,1,mad), decreasing = T),]) # 15000,10000,8000,5000

#apply(merged,1,mad) > max(quantile(apply(merged,1,mad), probs=seq(0, 1, 0.25))[2],0.01)

#merged <- t(fpkm[which(apply(fpkm,1,mad)> max(quantile(apply(fpkm,1,mad), 
                                                 #probs=seq(0, 1, 0.25))[2],0.01)),])
dim(WGCNA_matrix)
#dim(merged)

datExpr <- WGCNA_matrix 


#检查样本表达值是否OK
gsg<-goodSamplesGenes(datExpr,verbose = 3)
gsg$allOK  #应该为TRUE。
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(as.data.frame(datExpr))[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(as.data.frame(datExpr))[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}
#Flagging genes and samples with too many missing values...
#  ..step 1

collectGarbage()

#找到离群样本
sampletree <-hclust(dist(datExpr), method = "average")
sizeGrWindow(13,10)
#par(cex=0.6) 字体大小
par(mar=c(0,5,2,0))
plot(sampletree, main = "Sample clustering to detect outliers",sub = '',
     xlab="",cex.lab = 2, cex.axis=2, cex.main = 3,cex=2 )
dim(datExpr)
#clust = cutreeStatic(sampletree, cutHeight = 25000, minSize = 10)
#table(clust)
#keepSamples = (clust==1)
#datExpr = datExpr[keepSamples, ]

pheno <-read.csv('strand_specific/WGCNA/pheno1.csv',row.names = 1)

save(datExpr,pheno,file='data_input.Rdata')
load('data_input.Rdata')
#==================================================================================step2:确定最佳beta值

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, RsquaredCut = 0.85,
                        networkType = "unsigned",
                        powerVector = powers, verbose = 5)

sizeGrWindow(10,5)
par(mfrow = c(1,2))
par(mar=c(5,6,2,2))
cex1 = 1.5
#png("power-beta-value.png",width = 800,height = 600)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"),
     cex.lab = 1.5, cex.axis=1.5, cex.main = 2,cex=2 );
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
#abline(h=0.85,col="blue")
#abline(h=0.80,col="green")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"),
     cex.lab = 1.5, cex.axis=1.5, cex.main = 2,cex=2)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#dev.off()
graph2ppt(file='power-beta-value.pptx')


sft$powerEstimate     #此值即为最佳β次加权值。#rice=14
#软阈值也可以自己挑选，那么如何挑选？挑选R.sq的值尽量高（0.8），
#同时最大连通性mean.k.又不能太低（100）。

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

if (is.na(sft$powerEstimate)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))))}

#============================================================================step3：一步构建共表达矩阵
net = blockwiseModules(datExpr, power = sft$powerEstimate, #0.9   #最佳β次加权值
                       maxBlockSize = nGenes, minModuleSize = 30,  #设置模块最大最小基因数maxBlockSize = nGenes
                       TOMType = "unsigned", corType = "pearson",  #无标度网络
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs =F,saveTOMFileBase = "FPKM-TOM",
                       loadTOMs = T,verbose = 3)
table(net$colors) 
#一般来说，结果包含十几个到二十几个模块是比较正常的，此外一个模块中的基因数量不宜过多。
#像模块1的基因数量达到了大几千，就有点太多，这主要是因为powers不同导致的。（正常分析，不必在意）
#当然如果grey模块中的基因数量比较多也是不太好的，表示样本中的基因共表达趋势不明显，不同特征的样本之间差异性不大，或者组内基因表达一致性比较差。

#=====================================================================================step4: 模块可视化

# Convert labels to colors for plotting

mergedColors = labels2colors(net$colors)
table(mergedColors)

# 绘制树状图和模块颜色    # hclust for the genes.
expColor=t(numbers2colors(log2(datExpr+1),colors=blueWhiteRed(100),naColor="grey"))
colnames(expColor)=rownames(datExpr)
#png("genes-modules.png",width = 800,height = 600)
sizeGrWindow(14,9)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],# 将所有基因分配给其相应模块
                    "Module colors",
                    dendroLabels = F, hang = 0.03,#聚类树枝设置
                    frame.plot = F,guideCount = 40,#guideAll = T,
                    addGuide = T, guideHang = 0.05,#聚类树指示虚线
                    cex.colorLabels = 1.5, #cex.dendroLabels = 2.5,
                    cex.lab = 1.5,cex.main = 2,cex.axis=1.3 ,
                    marAll = c(1, 8, 3, 0.1))#下左上右

graph2ppt(file='genes-modules.pptx')
dev.off()

#with all samples
plotDendroAndColors(net$dendrograms[[1]], 
                    colors=cbind(mergedColors[net$blockGenes[[1]]],expColor),
                    c("Module",colnames(expColor)),
                    dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05,
                    cex.rowText=0.5)
dev.off()
#绘制两两模块间的邻接矩阵
png("wgcna.adjacency.heatmap.png",height = 1000,width = 900)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap",plotDendrograms = F,
                      marDendro = c(4,4,2,4))
dev.off()

modulelabels <-net$colors
moduleColors <-labels2colors(net$colors)
genetree <- net$dendrograms[[1]]
MEs <- net$MEs
save(modulelabels,moduleColors,genetree,MEs ,file = 'net_construction.Rdata')


load('strand_specific/WGCNA/Rdata/net_construction.Rdata')
load('strand_specific/WGCNA/Rdata/data_input.Rdata')
pheno <-read.csv('strand_specific/WGCNA/pheno-treatment.csv',row.names = 1)
#================================================================================step5:模块和性状的关系
## 最重要的一步。主要是针对于连续变量(性状)，如果是分类变量（分组信息），需转换成连续变量方可使用。
table(pheno)  
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
moduleColors <- labels2colors(modulelabels)  #每基因对应的颜色
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes # ME值
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs,pheno, use = "p")   #模块与性状相关性
colnames(MEs)
modules <-MEs
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)#模块与性状相关性的Pvalue
  
# 输出每个基因所在的模块，以及与该模块的KME值
file.remove('All_Gene_KME.txt')
for(module in substring(colnames(MEs),3)){
  if(module == "grey") next
  ME=as.data.frame(MEs[,paste("ME",module,sep="")])
  colnames(ME)=module
  datModExpr=datExpr[,moduleColors==module]
  datKME = signedKME(datModExpr, ME)
  datKME=cbind(datKME,rep(module,length(datKME)))
  write.table(datKME,quote = F,row.names = T,append = T,file = "All_Gene_KME.txt",col.names = F)
}

#绘制所有模块的表达值热图与特征值条形图
for(module in substring(colnames(MEs),3)){
  if(module == "grey") next
  ME=MEs[,paste("ME",module,sep="")]
  dev.off()
  #p(paste("wgcna.", module, ".express.barplot.png", sep=""),height = 700,width = 900)
  par(mfrow=c(2,1),mar=c(0.3,5.5,3,2))
  da=t(scale(datExpr[,moduleColors==module]))
  colname <-colnames(da)
  plotMat(t(scale(datExpr[,moduleColors==module])),
          clabels=colname,cex.lab = 1.5,
          rlabels=F,cex.main=2,cex=2,cex.axis=1.5)
  
  par(mar=c(5,5,0,0.7))
  barplot(ME,col=module,main="",cex.main=2,ylab="Eigengene Expression",
          xlab=paste(module, "Module"),
          cex.lab = 1.5, cex.axis=1.5, cex.main = 2,cex=2)
  graph2ppt(file=paste(module, '.express.barplot.pptx', sep=""))
  
}

#单个出图
module <- 'greenyellow'#更改这里
ME=MEs[,paste("ME",module,sep="")]
#png(paste("wgcna.", module, ".express.barplot.png", sep=""),height = 700,width = 900)
par(mfrow=c(2,1),mar=c(0.3,5.5,3,2))
da=t(scale(datExpr[,moduleColors==module]))
colname <-colnames(da)
plotMat(t(scale(datExpr[,moduleColors==module])),#nrgcols=30,rcols=module,
        rlabels=F,cex.main=2,clabels=colname,cex.lab = 1.5, 
        cex.axis=1.5,cex=2)
par(mar=c(5,5,0,0.7))
barplot(ME,col=module,main="",ylab="Eigengene Expression",
        xlab="Greenyellow Module",#更改这里
        cex.lab = 1.5, cex.axis=1.5, cex.main = 2,cex=2)
graph2ppt(file=paste(module, '.express.barplot.pptx', sep=""))
dev.off()

#----------------------------------------------------------------  
#各模块间的相关性图
sizeGrWindow(5,7.5);
par(cex = 1.5)
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
MET = orderMEs(MEs)
#MET = orderMEs(cbind(MEs, M))  #样本与模块关联相关性图
sizeGrWindow(6, 8) 
#png("Eigengene.png",width = 800,height = 600)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), 
                      plotDendrograms = T, xLabelsAngle = 90,marDendro = c(0,2.5,1,2), 
                      cex.lab = 1.5, cex.axis=1.5, cex.main = 2)
library(export)
graph2ppt(file='Eigengene.pptx')
dev.off()



sizeGrWindow(10,6)
#display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");  #signif 四舍五入
dim(textMatrix) = dim(moduleTraitCor)


#png("Module-pheno-relationships.png",width = 600,height = 800)
par(mar = c(4, 13, 1.5, 2.2))
#par(mar = c(6, 8.8, 3, 2.2))  
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(pheno),  
               yLabels = names(MEs),ySymbols = names(MEs),
               colorLabels = F,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,  
               xLabelsPosition = "bottom",
               xLabelsAngle = 45,
               xLabelsAdj = 1,
               yLabelsPosition = "left",
               zlim = c(-1,1),
               main = paste(""),
               cex.lab = 1.5, cex.axis=1.5, cex.main = 2)
graph2ppt(file='Module-pheno-treatment-relationship.pptx')
dev.off()



## 通过模块与各种表型的相关系数，可以很清楚的挑选自己感兴趣的模块进行下游分析了。


#=============================================================step6:感兴趣性状或样本的模块的具体基因分析
# 性状跟模块虽然求出了相关性，可以挑选最相关的那些模块来分析，
# 但是模块本身仍然包含非常多的基因，还需进一步的寻找最重要的基因。
# 所有的模块都可以跟基因算出相关系数，所有的连续型性状也可以跟基因的表达值算出相关系数。
# 如果跟性状显著相关的基因也跟某个模块显著相关，那么这些基因可能就非常重要。


#----------------------------------------------------首先计算模块与基因的相关性矩阵geneModuleMembership 
# names (colors) of the modules
modNames = substring(names(MEs), 3) #从第三个字符开始切分字符
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
## 算出每个模块跟基因的皮尔森相关系数矩阵，MEs是每个模块在每个样本里面的值

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))  #相关性p值
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM.", modNames, sep="")

nGenes = ncol(datExpr)
intmodules <- modNames
for (module in intmodules){
  moduleGenes = (moduleColors==module) #判断模块基因
  column = match(module, modNames) 
  geneInModule<-colnames(datExpr)[moduleGenes]
  filename=paste('modulegenes-',module,'-',nGenes,'.genes.txt',sep = '')
  write.table(as.data.frame(geneInModule),file = filename,row.names = F,col.names = F)
  print(length(geneInModule))
}

#-------------再计算性状或样本与模块的相关性矩阵geneTraitSignificance
## 只有连续型性状才能只有计算，这里把是否属于 M 这个样本用0,1进行数值化。
M=as.data.frame(pheno$LN.8d)
names(M)<-'LN.8d'
geneTraitSignificance = as.data.frame(cor(datExpr,M, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))#相关性p值
names(geneTraitSignificance) = paste("GS.", names(M), sep="")
names(GSPvalue) = paste("p.GS.", names(M), sep="")

intmodules <- modNames
for (module in intmodules){
#module = "blue"  
column = match(module, modNames)  #匹配blue所在列
moduleGenes =(moduleColors==module) #判断模块包括的基因，布尔值
names(M)<-'LN.8d'
sizeGrWindow(7, 7)
par(mfrow = c(1,1),mar=c(5,5,4,2))
pdf(paste('MM-GS-',module,'.pdf',sep=''),
    width = 800,height = 600)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module membership in", module, "module"),
                   ylab = paste("Gene significance for",names(M)),
                   main = paste("Module Membership vs. Gene Significance\n"),
                   cex.main = 2, cex.lab = 1.8, cex=2,
                   cex.axis = 1.5, col= module)
dev.off()
#这些基因不仅是跟其对应的模块高度相关，且跟其对应的性状高度相关，进一步说明这些基因值得深度探究。



#  基因与性状关系（GS）& 基因与模块关系（MM）
## 1 计算 module membership (MM): 基因（TMP）与模块（MEs）的相关性
##2 计算 Gene Significance (GS): 基因（TMP）与性状的相关性
#   
#=====================================================================================
{
  modNames = substring(names(MEs), 3)
  traitNames = names(pheno)
  
  # 计算 MM
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 
                                            nSamples))
  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="")
  
  # 计算 GS
  geneTraitSignificance = as.data.frame(cor(datExpr, pheno, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
                                            nSamples))
  names(geneTraitSignificance) = paste("GS.", traitNames, sep="");
  names(GSPvalue) = paste("p.GS.", traitNames, sep="");
  
  geneInfo<-cbind(geneModuleMembership, MMPvalue, geneTraitSignificance, GSPvalue)
  write.table(geneInfo, file = "geneInfo.txt", 
              sep = "\t", 
              quote = F)
}

save(datExpr, pheno, MEs, geneModuleMembership, geneTraitSignificance, 
     file = "trait_analysis.RData")



#模块基因的筛选
FilterGenes= abs(geneModuleMembership[moduleGenes, column])> 0.8 & abs(geneTraitSignificance[moduleGenes, 1])>0.2
table(FilterGenes)
GsAndMMGenes<-colnames(datExpr)[moduleGenes][FilterGenes]
filtername =paste('MM-GS-filter-',module,'.txt',sep='')
write.table(as.data.frame(GsAndMMGenes),filtername,row.names = F,col.names = F)
}

#====================================================================================step7:网络的可视化
#---------------------------------------------------------------------------------首先针对所有基因画热图
# 主要是可视化 TOM矩阵，WGCNA的标准配图
# 然后可视化不同模块的相关性热图、不同模块的层次聚类图
# 还有模块诊断，主要是 intramodular connectivity

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
geneTree = net$dendrograms[[1]] 
#读入构建时保存的TOM矩阵
load('FPKM-TOM-block.1.RData') 
TOM=as.matrix(TOM)
#TOM = TOMsimilarityFromExpr(datExpr, power =14 )#改变power
dissTOM = 1-TOM #非相似性矩阵
plotTOM = dissTOM^7                                
diag(plotTOM) = NA   #对角线值设为NA，使图更好看。
#png("Network-heatmap-allgenes.png",width = 800,height = 600)
#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes") 
#dev.off()

#15000基因，极为耗时，费时费资源，建议选取其中部分基因即可

nSelect = 2000
# 为了重现性，我们设置随机seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
sizeGrWindow(9,9)
plotDiss = selectTOM^7
diag(plotDiss) = NA
png("Network-heatmap-2000genes.png",width = 800,height = 600)
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
TOMplot(plotDiss, selectTree, selectColors, #terrainColors=TRUE, 
        main = "Network heatmap plot, selected genes",
        col = myheatcol)

dev.off()

#============================================================step8:提取指定模块的基因名
# 主要是关心具体某个模块内部的基因
#循环选择感兴趣的过滤模块

#library(clusterProfiler)
#library(maize.db)
library(ggplot2)
#maize.db<-loadDb(file = "maize.orgdb") 
                            
intmodules <- modNames
for (module in intmodules){
  moduleGenes = (moduleColors==module) #判断模块基因
  column = match(module, modNames) 
  geneInModule<-colnames(datExpr)[moduleGenes]
  filename=paste('modulegenes-',module,'.txt',sep = '')
  write.table(as.data.frame(geneInModule),file = filename,row.names = F,col.names = F)
  print(length(geneInModule))
}

#-----------------------------------------------------------------------


#单个模块
# Select module
module = "blue"
# Select module genenames
genename = colnames(datExpr) 
inModule = (moduleColors==module) #判断
modgenename = genename[inModule]  #提取
length(modgenename)
head(modgenename)


# 使用WGCNA包自带的热图很丑。
library(pheatmap)
dat=t(datExpr[,moduleColors==which.module ] ) 
#pheatmap(dat ,show_colnames =F,show_rownames = F) 
n=t(scale(t(log(dat+1)))) 
n[n> 2]= 2 
n[n< -2]= -2
n[1:4,1:4]
#pheatmap(n,show_colnames =F,show_rownames = F)
group_list=pheno
ac=data.frame(g=group_list)
rownames(ac)=colnames(n) 
dev.off()
png(paste('pheatmap-',which.module,'.png',sep = ''),width = 1200,height = 800)
pheatmap(n,show_colnames =T,show_rownames = F,angle_col = 315)
dev.off()
# 有了基因信息，可进行下游分析。包括GO/KEGG等功能数据库的注释。


#=============================================================Step9: 模块的导出

#读入构建时保存的TOM矩阵
load('strand_specific/WGCNA/Rdata/FPKM-TOM-block.1.RData') 
TOM=as.matrix(TOM)
#TOM = TOMsimilarityFromExpr(datExpr, power =sft$powerEstimate )#改变power

# Select module
intmodules <- modNames
for (module in intmodules){
# Select module genename
genename = colnames(datExpr) 
inModule = (moduleColors==module)
modgenename = genename[inModule]
## 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modgenename, modgenename)
#------------------------------------------------------------导出到cytoscape，常用
cyt = exportNetworkToCytoscape(modTOM,
  edgeFile = paste("CytoscapeInput-edges-0.2-", paste(module, collapse="-"), ".txt", sep=""), #边文件
  nodeFile = paste("CytoscapeInput-nodes-0.2-", paste(module, collapse="-"), ".txt", sep=""), #节点文件
  weighted = TRUE,threshold = 0.2, #权重阈值筛选,default=0.5
  nodeNames = modgenename, nodeAttr = moduleColors[inModule])
}
#如果模块包含的基因太多，网络太复杂，还可以进行筛选，比如：
nTop = 30  #前50基因
IMConn = softConnectivity(datExpr[, modgenename]);
top = (rank(-IMConn) <= nTop)
filterTOM <- modTOM[top, top]
cyt = exportNetworkToCytoscape(filterTOM,
        edgeFile = paste("CytoscapeInput-30edges-", paste(module, collapse="-"), ".txt", sep=""),
        nodeFile = paste("CytoscapeInput-30nodes-", paste(module, collapse="-"), ".txt", sep=""), 
        weighted = TRUE,threshold = 0.02,
        nodeNames = modgenename[top], nodeAttr = moduleColors[inModule][top])
#}

#所有模块的网络
module <- unique(moduleColors)
table(moduleColors)
genename = colnames(datExpr) 
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = "Cytoscape-edges-0.2-all-module.txt", #边文件
                               nodeFile = "Cytoscape-nodes-0.2-all-module.txt", #节点文件
                               weighted = TRUE,threshold = 0.2, #TOM权重阈值筛选,default=0.5
                               nodeNames = genename, #所有模块基因
                               nodeAttr = moduleColors)#所有模块基因对应的所有颜色

#install.packages("https://github.com/xuzhougeng/org.Osativa.eg.db/releases/download/v0.01/org.Osativa.eg.db.tar.gz", repos = NULL,  type="source")
