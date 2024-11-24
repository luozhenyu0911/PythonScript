rm(list = ls())
options(stringsAsFactors = F)

#===============================================================================About CPM and TPM
#Counts by featurecounts
countsdata<-read.table(file = paste0(getwd(), '/expr-matrix.txt'),
                       header = T,sep = '\t',quote = '',check.names = F)  
exprSet<-countsdata[,7:ncol(countsdata)]
rownames(exprSet)<-countsdata[,1]
prefix<-"counts"

#CPM
cpm <- t(t(exprSet)/colSums(exprSet) * 1000000)
avg_cpm <- data.frame(avg_cpm=rowMeans(cpm))
lcpm<-edgeR::cpm(exprSet,log=TRUE, prior.count=3)

#TPM 
kb <- countsdata$Length / 1000
rpk <- exprSet/ kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
avg_tpm <- data.frame(avg_tpm=rowMeans(tpm))
colSums(tpm)
#write.csv(avg_tpm,paste0(prefix,"_avg_tpm.csv"))
#write.csv(avg_cpm,paste0(prefix,"_avg_cpm.csv"))
write.csv(tpm,paste0(prefix,"_tpm.csv"))
write.csv(cpm,paste0(prefix,"_cpm.csv"))
#=====================================================================================================

#### https://github.com/AAlhendi1707/countToFPKM

#convert the feature counts of paired-end RNA-Seq into FPKM normalised values 
#by library size and feature effective length

library(countToFPKM)
####USEAGE EXAMPLE
#file.readcounts <- system.file("extdata", "RNA-seq.read.counts.csv", package="countToFPKM")
#file.annotations <- system.file("extdata", "Biomart.annotations.hg38.txt", package="countToFPKM")
#file.sample.metrics <- system.file("extdata", "RNA-seq.samples.metrics.txt", package="countToFPKM")

# Import the read counts matrix data into R.
countsdata<-read.table(file ='olddata/expr-matrix.txt',header = T,sep = '\t',quote = '',check.names = F)  
head(countsdata)
exprSet<-countsdata[,7:ncol(countsdata)]
rownames(exprSet)<-countsdata[,1]
counts <- exprSet

# Import feature annotations. 
# Assign feature lenght into a numeric vector.
gene.annotations <- countsdata[,1:6]
featureLength <- gene.annotations$Length  
head(featureLength)

# Import sample metrics. 
# Assign mean fragment length into a numeric vector.
    #meanFragmentLength 平均片段长度，可以用CollectInsertSizeMetrics(Picard) 计算。
    #https://broadinstitute.github.io/picard/command-line-overview.html#CollectInsertSizeMetrics
meanFragmentLength.metrics <- read.table('olddata/meanFragmentLength.txt', sep="\t", header=F)
meanFragmentLength <- meanFragmentLength.metrics$V3
meanFragmentLength

# Return FPKM into a numeric matrix.
fpkm_matrix <- fpkm(counts, featureLength, meanFragmentLength)  
head(fpkm_matrix)  #gene变少了？
write.csv(fpkm_matrix,'fpkm_matrix_by_R.csv')
#write.table(fpkm_matrix,file='fpkm_matrix.txt',quote = F,sep = '\t',row.names = T,col.names = T)

cor_fpkm_matrix <-cor(fpkm_matrix)
head(cor_fpkm_matrix)
pheatmap(cor_fpkm_matrix,display_numbers = T)

# Plot log10(FPKM+1) heatmap of top 30 highly variable features
fpkmheatmap(fpkm_matrix, topvar=30, showfeaturenames=TRUE, return_log = TRUE)
