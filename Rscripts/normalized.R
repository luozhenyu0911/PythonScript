rm(list=ls())
options(stringsAsFactors = F)


#============================================================================== Z-score
#所谓数据的标准化是指中心化之后的数据在除以数据集的标准差，
#即数据集中的各项数据减去数据集的均值再除以数据集的标准差,即Z-score.
#数据标准化是为了消除量纲对数据结构的影响。可以使用scale方法来对数据进行中心化和标准化.

#scale(x, center = TRUE, scale = TRUE) 默认对列col 标准化,对于基因表达矩阵需要转置。

#options(digits=3) #限定输出小数点后数字的位数为3位
tpm<-read.table('strand_specific/deg.uniq.fpkm.txt',row.names = 1,header = T,sep = '\t')
dim(tpm)
t_tpm <-data.frame(t(tpm))
dim(t_tpm)
#tpm_zscore <-apply(tpm, 1, scale)
zscore <-data.frame(scale(t_tpm,center = T,scale = T))
class(zscore)
tpm_zscore <-data.frame(t(zscore))
dim(tpm_zscore)
tpm_zscore[is.na(tpm_zscore)] <- 0 # NA -> 0
#tpm_zscore[tpm_zscore==0] = NA 
write.csv(tpm_zscore,file = 'deg.uniq.fpkm_zscore.csv')


#=============================================================================About CPM and TPM
#Counts
countsdata<-read.table(file = paste0(getwd(), '/counts/expr-matrix-s2.txt'),
                       header = T,sep = '\t',quote = '',check.names = F)  
exprSet<-countsdata[,7:ncol(countsdata)]
rownames(exprSet)<-countsdata[,1]
prefix<-"counts-sas"

#CPM   edger::logcpm
CPM <- t(t(exprSet)/colSums(exprSet) * 1000000)
#avg_cpm <- data.frame(avg_cpm=rowMeans(cpm))

#TPM 
kb <- countsdata$Length / 1000
rpk <- exprSet/ kb
TPM <- t(t(rpk)/colSums(rpk) * 1000000)
#avg_tpm <- data.frame(avg_tpm=rowMeans(tpm))
colSums(TPM)#应是一百万

#write.csv(avg_tpm,paste0(prefix,"_avg_tpm.csv"))
#write.csv(avg_cpm,paste0(prefix,"_avg_cpm.csv"))
write.csv(TPM,paste0(prefix,"_TPM.csv"))
write.csv(CPM,paste0(prefix,"_CPM.csv"))

#----------------------------------------FPKM to TPM ，limma分析。
#重点就是：RPKM/FPKM矩阵可以转为TPM后,再使用limma进行差异分析!
fpkmexpMatrix <-read.table(file = paste0(getwd(), '/data/fpkmmerge.txt'),row.names = 1,
                                   header = T,sep = '\t',quote = '',check.names = F)  
fpkmToTpm <- function(fpkm)
{exp(log(fpkm) - log(sum(fpkm)) + log(1e6))}
tpms <- apply(fpkmexpMatrix,2,fpkmToTpm)  #fpkmexpMatrix 替换成上述rpk，结果一样
colSums(tpms)

#------------------------------------------------------
#! /usr/bin/env Rscript

## functions for rpkm and tpm
## from https://gist.github.com/slowkow/c6ab0348747f86e2748b#file-counts_to_tpm-r-L44
## from https://www.biostars.org/p/171766/

#此处RPKM即为FPKM。
rpkm <- function(counts, lengths) {   
  rate <- counts / lengths
  rate / sum(counts) * 1e9
}

#计算FPKM
FPKM_cal <- function(counts,lengths){
  rate<-counts/lengths
  return(rate/sum(counts)*10^9)
}

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

## read table from featureCounts output
ftr.cnt <- read.table('ssRNAseq/expression/count.sas.txt', sep="\t", 
                      stringsAsFactors=FALSE, header=TRUE)

library(dplyr)
library(tidyr)

ftr.rpkm <- ftr.cnt %>%
  gather(sample, cnt, 7:ncol(ftr.cnt)) %>%
  group_by(sample) %>%
  mutate(rpkm=rpkm(cnt, Length)) %>%
  select(-cnt) %>%
  spread(sample, rpkm)
write.table(ftr.rpkm, file= "sas_fpkm.txt", sep="\t", row.names=FALSE, quote=FALSE)

ftr.tpm <- ftr.cnt %>%
  gather(sample, cnt, 7:ncol(ftr.cnt)) %>%
  group_by(sample) %>%
  mutate(tpm=tpm(cnt, Length)) %>%
  select(-cnt) %>%
  spread(sample, tpm)
write.table(ftr.tpm, file="sas_tpm.txt", sep="\t", row.names=FALSE, quote=FALSE)
colSums(ftr.tpm[7:ncol(ftr.cnt)])#应是一百万


#----------------------------------------------------------？？？？

countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))}

countToFpkm <- function(counts, effLen)
{N <- sum(counts)
exp( log(counts) + log(1e9) - log(effLen) - log(N) )}

fpkmToTpm <- function(fpkm)
{exp(log(fpkm) - log(sum(fpkm)) + log(1e6))}

countToEffCounts <- function(counts, len, effLen)
{counts * (len / effLen)}


countsdata<-read.table(file = paste0(getwd(), '/expr-matrix.txt'),
                       header = T,sep = '\t',quote = '',check.names = F)  
cnts <- countsdata[,7:ncol(countsdata)]
rownames(cnts)<-countsdata[,1]
lens <- countsdata$Length
countDf <- data.frame(count = cnts, length = lens)
count <-cnts

countDf$effLength <- countDf$length  
countDf$tpm <- with(countDf, countToTpm(count, effLength))
tpm_expr <- data.frame(countDf$tpm)
colSums(tpm_expr)
countDf$fpkm <- with(countDf, countToFpkm(count, effLength))
with(countDf, all.equal(tpm, fpkmToTpm(fpkm)))
countDf$effCounts <- with(countDf, countToEffCounts(count, length, effLength))

write.csv(countDf,file = 'countDf.csv')

