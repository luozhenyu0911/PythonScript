rm(list = ls())
options(stringsAsFactors = F)

library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(clusterProfiler)
library(org.Osativa.eg.db)
#library(GenomicAlignments)
#library(GenomicFeatures)
#library(annotation)
#library(AnnotationHub)
library(ggrepel)
library(cowplot)
library(pheatmap)
library(corrplot)
library(export)
#register(MulticoreParam(4))#服务器多线程加速

countsdata<-read.table('counts/expr-matrix-s2.txt',row.names = 1,header = T,sep='\t')[,6:9]
head(countsdata)
countsdata <-na.omit(countsdata)
countdata<-countsdata
#countdata <- countdata[rowSums(countdata)>0,]
head(countdata)

#level中 Control在Treat组前
coldata<-data.frame(Condition = factor(c(rep('Nip',2),rep('B114',2)),levels = c('Nip','B114')),
                    row.names = colnames(countdata))  
coldata

#--------------------------------RNAseq样本间重复相关性分析可用FPKM,TPM.(最好)
#也可用read counts计算（尽量不用）。

if(T){
fpkm<- read.table('sense_fpkm.txt',header = T,sep = '\t',quote = '',
                   check.names = F,row.names = 1)
#cor_matrix <-cor(fpkm) #,method ="spearman"
cor_matrix <-cor(fpkm)

pheno<-data.frame(Condition = factor(c(rep('Nipponbare',2),rep('114',2)),
                                     levels = c('Nipponbare','114')),
                  row.names = colnames(fpkm))
  
p <- pheatmap(cor_matrix, fontsize=18,border_color = NA,#angle_col=45,
              #color=colorRampPalette(colors = c('#436eee','white','#ee0000'))(10),
              #color=colorRampPalette(colors = c('white','#ee0000'))(100)[30:100],
              #color=colorRampPalette(colors = c('cyan1','black','yellow'))(100),
              color=colorRampPalette(colors = c('#0408f2','#ffff00'))(100),
              display_numbers = T,annotation_col = pheno ,annotation_row = pheno,
              annotation_names_row = F,annotation_names_col = F,
              cellwidth = 60,cellheight = 60 #gaps_col = c(8)
              )
p
dev.off()
}
p
graph2ppt(file="sense-sample-correlation.pptx")
#-------------------------------------------------------------------------------------------------------

dds<-DESeqDataSetFromMatrix(countData=countdata,colData = coldata,
                            design =~Condition) # 批次效应design =~batch+Condition
dds  
#dds <- dds[rowSums(counts(dds))> 1,]    
dds<-DESeq(dds)

# 标准化后的数据
normalized_counts <- counts(dds, normalized= TRUE)
head(normalized_counts)
#normalized_counts_mad <- apply(normalized_counts, 1, mad)
#normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]

resultsNames(dds) 
res<-results(dds,name ="Condition_B114_vs_Nip")  # a vs b: 则a为处理组，b为control组。
head(res)

#res$padj[is.na(res$padj)] <- 1 # 校正后padj为NA的赋值为1
res<-res[order(res$padj),]

diff_gene <- subset(res,padj<0.05 & (log2FoldChange > 1 | log2FoldChange < -1))  
diffdata <- as.data.frame(diff_gene)
#write.csv(diffdata,'diffdata-1.csv')  

gene_up <- row.names(diff_gene[diff_gene$log2FoldChange > 1, ])
length(gene_up)   
#write.table(x=gene_up,file='gene_up.txt',quote = F,sep = '\t',row.names = T,col.names = T)
gene_down <- row.names(diff_gene[diff_gene$log2FoldChange < -1, ])
length(gene_down) 
#write.table(x=diffdata,file='diffdata.txt',quote = F,sep = '\t',row.names = T,col.names = T)

resdata<-merge(as.data.frame(res), as.data.frame(counts(dds,normalized=F)), 
               by='row.names',sort=F)                     
head(resdata)

resdata$significant <- 'unchanged'  
resdata$significant[resdata$padj < 0.05 & resdata$log2FoldChange> 1] <- 'upregulated'     
resdata$significant[resdata$padj < 0.05 & resdata$log2FoldChange< -1] <- 'downregulated'  
head(resdata)


DEG=as.data.frame(res)
DEG=na.omit(DEG)  
nrDEG=DEG 
##logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) ) # logFC_cutoff=1 
nrDEG$significant = as.factor(ifelse(nrDEG$padj < 0.05 & abs(nrDEG$log2FoldChange) >1, 
                       ifelse(nrDEG$log2FoldChange > 1 ,'Up in 114','Up in Nipponbare'),
                       'Not Significant')) 
#this_tile <- paste0("24H-LN vs CK", 
 #                   '\nup gene : ',nrow(nrDEG[nrDEG$significant =='Up',]) , 
  #                  '\ndown gene : ',nrow(nrDEG[nrDEG$significant =='Down',])) 
threshold <- nrDEG$significant
valcano <- ggplot(data=nrDEG, aes(x=log2FoldChange, y=-log10(padj), 
                                  color=significant)) + 
  geom_point(alpha=0.5,size=5) + #size=abs(nrDEG$log2FoldChange)
  #scale_size(range=c(15,18))+
  theme_set(theme_set(theme_classic(base_size=18)))+ 
  xlim(c(-14, 14))+ 
  scale_fill_discrete(labels=c("Up in 114", "Up in Nipponbare", "Not Significant"))+ #ylim(c(-15, 150))+
  xlab("log2FoldChange") + ylab("-log10(padj)")+ 
  #theme(legend.position = c(0.9,0.8))+ ggtitle("") + 
  theme(plot.title = element_text(size=15,hjust = 0.5))+  
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  #theme_bw()+theme(panel.background=element_rect(fill='transparent',color='black'))+
  geom_vline(xintercept=c(-1,1), lty=2, col="gray30", lwd=0.5) + 
  geom_hline(yintercept=-log10(0.05), lty=2, col="gray30", lwd=0.5)+
  scale_colour_manual(values = c('#eb5352','#319d42','grey60'),
                      limit=c("Up in 114", "Up in Nipponbare", "Not Significant"))+
  theme(plot.caption = element_text(vjust  =  1)) + 
  theme(axis.line = element_line(size  =  0.8, linetype  =  'solid')) + 
  theme(axis.ticks = element_line(size  =  0.8)) + 
  theme(axis.title = element_text(size  =  17, face  =  'bold')) + 
  theme(axis.text = element_text(size  =  16, colour  =  'black')) + 
  theme(axis.text.x = element_text(size  =  16)) + 
  theme(axis.text.y = element_text(size  =  16)) + 
  theme(legend.text = element_text(size  =  14)) + 
  theme(legend.title = element_text(size = 16, face  =  'bold')) + 
  theme(legend.position = c(0.18, 0.88),
        legend.background = element_rect(colour = 'black')) + 
  labs(x = "log2FoldChange", y ='-log10(padj)',color='Change')

valcano

library(ggThemeAssist)
ggThemeAssistGadget(valcano)

dev.off()
table(resdata$significant)
#write.csv(normalized_counts,file="normalized_counts-LN-96Hvs48H.csv")
write.csv(resdata,'sense-resdata-114vsNIP.csv')

ggsave(valcano,filename = 'antisense-volcano-114vsNIP.png')



#------------------------------------------------------------------------------------  GO功能富集
#hub <- AnnotationHub()    
#q <- query(hub,'orizy')  
#q
#id <- q$ah_id[length(q)]    
#rice <- hub[[id]]    
if(T){
database <- org.Osativa.eg.db
fc <- nrDEG$log2FoldChange
names(fc) <-  rownames(nrDEG)

ego_up <- enrichGO(gene_up,       #水稻GO需要MSU id: LOC_Os01g50820
                OrgDb = database,
                ont = 'ALL',          #可选BP,CC,MF
                keyType = 'GID',  
                pAdjustMethod = 'BH',  minGSSize = 1,maxGSSize = 500,  
                pvalueCutoff = 0.1,
                qvalueCutoff = 0.1,
                pool = T) 
head(ego_up[,1:8])
p1 <- dotplot(ego_up,font.size=5)
h1<-heatplot(ego_up,foldChange = fc)#,showCategory = 10
}

if(T){
ego_down <- enrichGO(gene_down,       #ENSEMBL, GO一般使用ENSEMBL id
                   OrgDb = database,
                   ont = 'ALL',          #可选BP,CC,MF
                   keyType = 'GID',  
                   pAdjustMethod = 'BH',minGSSize = 1,maxGSSize = 500,   
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.1,
                   pool = T) 
head(ego_down[,1:8])
p2 <- dotplot(ego_down,font.size=5)
h2<- heatplot(ego_down,foldChange = fc)#,showCategory = 10
}
#作图前不能把ego结果保存为数据框，会报错。

#----------------------------------------------------------------------------------KEGG Enrichment
if(T){
spe <- 'dosa'

rap_id <- mapIds(x = database, keys = gene_up, 
                 column = "RAP","GID")
rap_id <- paste0(rap_id[!is.na(rap_id)], "-01")
rap_id <- gsub("g","t",rap_id)

kkup <- enrichKEGG(gene = rap_id,    #水稻kegg一般为RAP id。
                   organism = spe,   
                   keyType = 'kegg', 
                   pAdjustMethod = 'BH',
                   minGSSize = 10,    #or 1,根据需要
                   maxGSSize = 500, 
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.5,
                   use_internal_data = F)
head(kkup[,1:8])
bb1<-barplot(kkup,showCategory = 10)    #展示10个富集结果
dd1<-dotplot(kkup)

rice_kegg <- clusterProfiler::download_KEGG("dosa")

kegg_df <- rice_kegg$KEGGPATHID2EXTID
kegg_df <- kegg_df[kegg_df$to %in% rap_id,]
kegg_df <- merge(kegg_df, rice_kegg$KEGGPATHID2NAME,
                 by.x="from",by.y="from")

kegg_class <- as.data.frame(sort(table(kegg_df$to.y), decreasing = T)[1:10])
colnames(kegg_class) <- c("pathway","times")
keggupplot<-ggplot(kegg_class,aes(x=pathway, y = times)) +
  geom_bar(fill="#ca0020",stat="identity") + coord_flip() +
  theme_bw() + geom_text(aes(y = times+1, label = times))


rap_id1 <- mapIds(x = database, keys = gene_down, 
                 column = "RAP","GID")
rap_id1 <- paste0(rap_id1[!is.na(rap_id1)], "-01")
rap_id1 <- gsub("g","t",rap_id1)

kkdown <- enrichKEGG(gene = rap_id1,    #水稻kegg一般为RAP id。
                   organism = spe,   
                   keyType = 'kegg', 
                   pAdjustMethod = 'BH',
                   minGSSize = 10,    #or 1,根据需要
                   maxGSSize = 500, 
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.5,
                   use_internal_data = F)
head(kkdown[,1:8])
bb2<-barplot(kkdown,showCategory = 10)    #展示10个富集结果
dd2<-dotplot(kkdown)


kegg_df1 <- rice_kegg$KEGGPATHID2EXTID
kegg_df1 <- kegg_df1[kegg_df1$to %in% rap_id1,]
kegg_df1 <- merge(kegg_df1, rice_kegg$KEGGPATHID2NAME,by.x="from",by.y="from")

kegg_class1 <- as.data.frame(sort(table(kegg_df1$to.y), decreasing = T)[1:10])
colnames(kegg_class1) <- c("pathway","times")
keggdownplot <- ggplot(kegg_class1,aes(x=pathway, y = times)) +
  geom_bar(fill="#2b83ba",stat="identity") + coord_flip() +
  theme_bw() + geom_text(aes(y = times+1, label = times))


rice_kegg <- clusterProfiler::download_KEGG("dosa")
kegg_df <- rice_kegg$KEGGPATHID2EXTID
# 获取上调的KEGG PATHWAY
up_rap_id <- mapIds(x = database, keys = gene_up, column = "RAP","GID")
up_rap_id <- paste0(up_rap_id[!is.na(up_rap_id)], "-01")
up_rap_id <- gsub("g","t",up_rap_id)

up_df <- rice_kegg$KEGGPATHID2EXTID
up_df <- kegg_df[kegg_df$to %in% up_rap_id,]
up_df <- merge(kegg_df, rice_kegg$KEGGPATHID2NAME,
               by.x="from",by.y="from")

# 获取下调的KEGG PATHWAY
down_rap_id <- mapIds(x = database, keys = gene_down,  column = "RAP", "GID")
down_rap_id <- paste0(down_rap_id[!is.na(down_rap_id)], "-01")
down_rap_id <- gsub("g","t",down_rap_id)
down_df <- rice_kegg$KEGGPATHID2EXTID
down_df <- down_df[down_df$to %in% down_rap_id,]
down_df <- merge(down_df, rice_kegg$KEGGPATHID2NAME,
                 by.x="from",by.y="from")

# 统计
kegg_class_up <- as.data.frame(sort(table(up_df$to.y),  decreasing = T)[1:10])
kegg_class_down <- as.data.frame(sort(table(down_df$to.y),  decreasing = T)[1:10])
#合并
kegg_class <- rbind(kegg_class_up, kegg_class_down)
colnames(kegg_class) <- c("pathway","times")
kegg_class$source <- rep(c("up","down"),times=c(nrow(kegg_class_up),nrow(kegg_class_down)))
#画图
kegg_up_down <-ggplot(kegg_class,aes(x=pathway, y = times)) +
  geom_bar(aes(fill=source),stat="identity",position = "dodge") + 
  scale_fill_manual(values = c(up="red",down="blue")) +coord_flip() + theme_bw()
}


write.csv(as.data.frame(ego_up),file = 'ego_up-LNvsCK-24h.csv')
write.csv(as.data.frame(ego_down),file = 'ego_down-LNvsCK-24h.csv')
write.csv(as.data.frame(kkup),file = 'kegg_up-LNvsCK-24h.csv')
write.csv(as.data.frame(kkdown),file = 'kegg_down-LNvsCK-24h.csv')
dev.off()
ggsave(ph,filename = 'gene-heatmap-LNvsCK-24h.png')
ggsave(v,filename = 'volcano_with_filtergene-LNvsCK-24h.png') 

b1
h1
b2
h2
bb1
dd1
bb2
dd2
keggupplot
keggdownplot
kegg_up_down

godot<-plot_grid(p1,p2)
gpbar<-plot_grid(h1,h2)
keggdot<-plot_grid(bb1,bb2)
keggbar<-plot_grid(dd1,dd2)
godot
gpbar
keggdot
keggbar
#save(ph,valcano,b1,h1,b2,h2,bb1,dd1,bb2,dd2,godot,gobar,keggdot,keggbar,keggupplot,keggdownplot,kegg_up_down,file = '24h-LNvsCK.Rdata')
save(ph,valcano,b1,h1,b2,bb1,dd1,bb2,dd2,keggdot,keggbar,keggupplot,keggdownplot,
     kegg_up_down,file = '24h-LNvsCK.Rdata')
