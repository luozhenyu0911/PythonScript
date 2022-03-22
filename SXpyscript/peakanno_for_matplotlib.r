library(ChIPseeker)
require(ChIPseeker)
require(GenomicFeatures)

#########peak annotation

args=commandArgs(T)
rice <- loadDb(args[1])
x1 = annotatePeak(args[2], tssRegion=c(-1000, 0),level="gene",overlap="all", TxDb=rice,assignGenomicAnnotation = TRUE,genomicAnnotationPriority = c("Promoter","5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"))
write.table(x1,args[3],sep="\t")
a<-as.data.frame(table(gsub("\\(.*\\)","",as.data.frame(x1)$annotation)))
write.table(a,args[4],sep="\t",quote = FALSE,row.names = FALSE, col.names = FALSE)
