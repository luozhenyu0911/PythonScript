#富集分析：
  #A：ORA  差异基因富集分析（只需要导入差异基因，不需要全部基因）         GO/KEGG
  #B: FCS  基因集(gene set)富集分析（不管有无差异，需要全部genes表达值）  GSEA

#一代ORA富集分析 (Over Represention Analysis)
   #ORA只是对差异基因感兴趣，而对差异基因是否上调或下调不关注。
#其目的是要判断我们找到的差异基因更多是在目标注释集中，而不是非注释基因组。
   #ORA的缺点 ：#1.已经选出了DEGs，需要主观的过滤。  
               #2.在统计检验的时候不考虑基因的表达情况  
              #3.一些微弱的却具有效力的基因集被过滤掉了。

#------------------------------------------------------------------------------------  GO功能富集
library(clusterProfiler)
library(org.Osativa.eg.db)

  # -----------------------------------------------------------------------------id type translate
library(GenomicAlignments)
library(GenomicFeatures)
library(annotation)

keytypes(org.Osativa.eg.db)
for (i in keytypes(org.Osativa.eg.db)) {print(head(keys(org.Osativa.eg.db, keytype = i)))}
head(keys(org.Osativa.eg.db, keytype = 'GID'))
keys(org.Osativa.eg.db)

resdata<-read.csv(file = 'resdata-LNvsCK-24h.csv',row.names = 1)  #使用DEseq得到的 resdata结果。
head(resdata)
     #列出、查看包中feature id类型
translate_db <- select(org.Osativa.eg.db, keys=as.character(resdata$Row.names), 
                       columns=c('ENTREZID','SYMBOL','GO'), keytype='ENSEMBL')
head(translate_db)

#ID转换时，ENSEMBL ID如:ENSMUSG00000096474.2 ，其中.2表示版本，
#在转换时应去掉小数点及其后数字,应为ENSMUSG00000096474 ，否则报错，无法转换。
#针对有版本号的ensembl_id的转换
     # library(stringr); data$ensembl_id <- str_split(data,'[.]',simplify = T)[,1]
# or 
     # data$ensembl_id <- unlist(lapply(data,function(x){strsplit(as.character(x),'[.]')[[1]][1]}))    



#GO富集：一个基因的功能可以从三个角度定义：MF分子功能 Molecular Function; 
                                          #CC细胞成分 Cellular Component; 
                                          #BP生物学进程 Biological Process。
#一般可以三类混合分析 'ALL '，也可单独分析某一类 BP/CC/MF,  三者中BP最常用。

#目前19类物种有 org.db 的feature id 转换包：http://bioconductor.org/packages/release/BiocViews.html  # AnnotationData > PackageType > OrgDb
  #若不存在于19个GO注释包之间，则通过AnnotationHub自行构建：
#BiocManager::install("AnnotationHub")
library(AnnotationHub)
hub <- AnnotationHub()      #构建所有，耗时长
q <- query(hub,'sativa')  #在hub中搜索需要的物种关键词，如水稻sativa
q
id <- q$ah_id[length(q)]    #此为选择最后一个行名来构建，也可选别的。如已知某行名，则直接：Musculus <- hub[['AH70580']]
Musculus <- hub[[id]]       #双切片索引至sativa变量，其即为水稻的注释包。
      #--------------------------------------------------------------------------------------------------------------------

database <- org.Hs.eg.db
ensembl <- as.character(resdata$Row.names[resdata$significant != 'unchanged']) #提取差异表达的ensembl id，需要factor to character
ids <- bitr(ensembl,fromType = 'ENSEMBL',toType = c('ENTREZID','SYMBOL'),OrgDb = database)   #bitr(),clusterProfiler包中转换id函数
#ids与此作用相同：translate_db <- select(org.Hs.eg.db,keys=as.character(resdata$Row.names),columns=c('ENTREZID','SYMBOL'),keytype='ENSEMBL')
head(ids)

ego <- enrichGO(gene = ids[,1],       #ENSEMBL, GO一般使用ENSEMBL id
                OrgDb = database,
                ont = 'ALL',          #可选BP,CC,MF
                keyType = 'ENSEMBL',  #此处要为ids[,1]列名类型。若gene = ids[,2],则为ENTREZID
                pAdjustMethod = 'BH',minGSSize = 1,maxGSSize = 500,    #default is ok
                pvalueCutoff = 0.1,
                qvalueCutoff = 0.1,
                pool = T,             # ont为ALL时，pool为T。若ont为单个时，不需要此参数。
                readable = T) 
head(ego[,1:8])

ego <-as.data.frame(ego)
write.table(ego,file = 'ego.txt',quote = F,sep = '\t',row.names = F,col.names = T)

#--------------------------------------------------------------------------------------------------plot
  #作图前不能把ego结果保存为数据框，会报错。
barplot(ego,showCategory = 10)    #展示10个富集结果
library(ggplot2)
ggsave('go_barplot.png')          #save picture

cnetplot(ego)   #网络图
heatplot(ego)#,showCategory = 10

library(upsetR)
upsetplot(ego)

dotplot(ego,font.size=5)        #气泡图
enrichMap(ego, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai) #网络图
plotGOgraph(ego)                #GO图


#-----------------------------------------------------------------------------------------KEGG Enrichment
#一般接触较多的是KEGG Pathway数据集

#查看物种的KEGG注释包：https://www.genome.jp/kegg/catalog/org_list.html 
#里面常用模式生物缩写：人 hsa ;鼠 mmu
spe <- 'hsa'
kk <- enrichKEGG(gene = ids[,2],    #ENTREZID,  KEGG一般为ENTREZID
                 organism = spe,    # or 'hsa'
                 keyType = 'kegg', 
                 pAdjustMethod = 'BH',
                 minGSSize = 10,    #or 1,根据需要
                 maxGSSize = 500, 
                 pvalueCutoff = 0.1,
                 qvalueCutoff = 0.5,
                 use_internal_data = F)
head(kk[,1:8])

#-----------------------------------------------------------------plot
barplot(kk,showCategory = 10)    #展示10个富集结果
ggsave('kegg_barplot.png') 
dotplot(kk,showCategory = 10)

barplot(kk,showCategory = 25, title="The KEGG enrichment analysis of all DEGs")+
        scale_size(range=c(2, 12)) + scale_x_discrete(labels=function(kk) str_wrap(kk,width = 25))
dotplot(kk,showCategory = 25, title="The KEGG enrichment analysis of all DEGs")+
        scale_size(range=c(2, 12)) + scale_x_discrete(labels=function(kk) str_wrap(kk,width = 25))

browseKEGG(kk,'mmu01100')  # 显示通路图

#-------------------------------------------------------------------------------------------------ggplot2简介
library(ggplot2)
gg <-data.frame(pathway=c('test1','test2','test3','test4'),
                pvalue=c(0.01,0.02,0.03,0.01),
                generatio=c(0.05,0.06,0.08,0.01))
gg
ggplot(data = gg,    #必要的2个参数：data, mapping=aes(此中属性可改可选)
       mapping = aes(x=pvalue,y=pathway,size=generatio)) + 
       geom_point()   #  geom_point()  指明作图类型，此为点图
ggplot(data = gg, mapping = aes(x=pvalue,y=pathway,color=generatio)) + geom_point()   
ggplot(data = gg, mapping = aes(x=pvalue,y=pathway,color=generatio,size=generatio)) + geom_point()   #以generatio值区分size和color
