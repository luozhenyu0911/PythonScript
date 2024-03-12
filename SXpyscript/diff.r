library("DESeq2")
args=commandArgs(T)
sample1<-as.numeric(args[4])
sample2<-as.numeric(args[5])
sum<-sample1+sample2
countTRData<-as.matrix(read.csv(args[1],row.names = "gene_id"))
condition<-factor(c(rep("group1",sample1),rep("group2",sample2)))
colData<-read.csv(args[2],row.names = 1)
colData_TPvsDL<-data.frame(colData[1:sum,],condition=condition)
TRdds <- DESeqDataSetFromMatrix(countData = countTRData[,1:sum], colData = colData_TPvsDL, design = ~ condition)
TRdds <- DESeq(TRdds)
TRres <- results(TRdds)
#TRres <- TRres[order(TRres$padj), ]
TRresdata <- merge(as.data.frame(TRres), as.data.frame(counts(TRdds, normalized=TRUE)), by="row.names", sort=FALSE)
names(TRresdata)[1]<-"Gene"
write.table(TRresdata, file=args[3],sep="\t",quote=FALSE)
