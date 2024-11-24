
#================清空当前环境,加载软件包=======
rm(list=ls())
library(topGO)
library(Rgraphviz)

#==============设置输入文件====================
input="DEGlist.txt"  #差异基因名称的列表
mapfile="Zea_mays.AGPv4.GO.list"    #所有基因GO map结果

#==============开始分析========================
geneID2GO<-readMappings(file=mapfile) #topGO自带函数
geneNames<- names(geneID2GO)
myInterestingGenes<- read.table(input)[,1]

geneList<-rep(1,length(geneID2GO))
names(geneList) <- names(geneID2GO)
geneList[match(myInterestingGenes,names(geneList))]=0
topDiffGenes<-function(allScore){return(allScore<0.01)}

#BP节点的富集分析
sampleGOdata <- new("topGOdata",nodeSize = 6,ontology="BP", 
                    allGenes = geneList, annot = annFUN.gene2GO, 
                    gene2GO = geneID2GO,geneSel=topDiffGenes)

#使用fisher检验
result <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
allRes <- GenTable(sampleGOdata,Fisher = result, orderBy = "Fisher", 
                   topNodes = attributes(result)$geneData[4])

gene=NULL
GOs=sampleGOdata@graph@nodeData@data
for( i in allRes$GO.ID){
  e=paste('g=names(GOs$`',i,'`$genes)',sep="")
  eval(parse(text=e))
  d=intersect(g,myInterestingGenes)
  e=paste('gene[\'',i,'\']=toString(d)',sep="")
  eval(parse(text=e))
}

allRes$gene=gene
write.table(allRes, file=paste(input,".BP.xls",sep=""), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
pdf(paste(input,".BP.pdf",sep=""))
showSigOfNodes(sampleGOdata, score(result), firstSigNodes = 10, useInfo = "all")
dev.off()
png(paste(input,".BP.png",sep=""))
showSigOfNodes(sampleGOdata, score(result), firstSigNodes = 10, useInfo = "all")
dev.off()

#MF节点的分析
sampleGOdata <- new("topGOdata",nodeSize = 6,ontology="MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO,geneSel=topDiffGenes)
result <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
allRes <- GenTable(sampleGOdata,Fisher = result, orderBy = "Fisher", topNodes = attributes(result)$geneData[4])
gene=NULL
GOs=sampleGOdata@graph@nodeData@data
for( i in allRes$GO.ID){
  e=paste('g=names(GOs$`',i,'`$genes)',sep="")
  eval(parse(text=e))
  d=intersect(g,myInterestingGenes)
  e=paste('gene[\'',i,'\']=toString(d)',sep="")
  eval(parse(text=e))
}
allRes$gene=gene
write.table(allRes, file=paste(input,".MF.xls",sep=""), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
pdf(paste(input,".MF.pdf",sep=""))
showSigOfNodes(sampleGOdata, score(result), firstSigNodes = 10, useInfo = "all")
dev.off()
png(paste(input,".MF.png",sep=""))
showSigOfNodes(sampleGOdata, score(result), firstSigNodes = 10, useInfo = "all")
dev.off()


#CC节点的富集分析
sampleGOdata <- new("topGOdata",nodeSize = 6,ontology="CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO,geneSel=topDiffGenes)
result <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
allRes <- GenTable(sampleGOdata,Fisher = result, orderBy = "Fisher", topNodes = attributes(result)$geneData[4])
gene=NULL
GOs=sampleGOdata@graph@nodeData@data
for( i in allRes$GO.ID){
  e=paste('g=names(GOs$`',i,'`$genes)',sep="")
  eval(parse(text=e))
  d=intersect(g,myInterestingGenes)
  e=paste('gene[\'',i,'\']=toString(d)',sep="")
  eval(parse(text=e))
}
allRes$gene=gene
write.table(allRes, file=paste(input,".CC.xls",sep=""), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
pdf(paste(input,".CC.pdf",sep=""))
showSigOfNodes(sampleGOdata, score(result), firstSigNodes = 10, useInfo = "all")
dev.off()
png(paste(input,".CC.png",sep=""))
showSigOfNodes(sampleGOdata, score(result), firstSigNodes = 10, useInfo = "all")
dev.off()

