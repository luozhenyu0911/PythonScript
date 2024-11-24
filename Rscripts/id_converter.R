rm(list = ls())
options(stringsAsFactors = F)


library(AnnotationHub)
ah <- AnnotationHub()      #构建所有，耗时长
unique(ah$species)        #查看AnonotationHub里面包括那些物种
unique(ah$rdataclass)    #看看注释格式类型,OrgDB只是其一种数据类型,但是clusterprofiler必需的类型。
display(ah)              
q <- query(ah,'Setaria_italica')  #在hub中搜索需要的物种关键词，如鼠 musculus
q
display(q)
id <- q$ah_id[length(q)]
org.Sitalica.eg.db<- ah[[id]]  #双切片索引至musculus变量，其即为鼠的注释包。
keytypes(org.Sitalica.eg.db)

keys(org.Sitalica.eg.db, keytype = "REFSEQ")  
keytypes() 

#--------------------------------------------------------------------------------
library(org.Osativa.eg.db)
library(AnnotationDbi)
msu <-read.table('strand_specific/counts-sense_TPM.csv',header = T)[,1]
length(msu)

keytypes(org.Osativa.eg.db)
for (i in keytypes(org.Osativa.eg.db)){
  print(head(keys(org.Osativa.eg.db, keytype = i)))
}

translate_db <- select(org.Osativa.eg.db,                            
                       keys=as.character(msu),    
                       columns=c('GO'),          
                       keytype='GID')                        #根据keytype来搜索columns
head(translate_db)
translate_db <-na.omit(translate_db)
head(translate_db)
write.csv(translate_db,file = 'msu2GOid.csv')


#--------------------------------------------------------------------------------------------------
library(clusterProfiler)
database<-org.Osativa.eg.db
ids <- bitr(msu,fromType = 'GID',toType = c('GO','RAP'),
            OrgDb = database)   #bitr(),clusterProfiler包中转换id函数
#ids与此作用相同：translate_db <- select(org.Hs.eg.db,keys=as.character(resdata$Row.names),columns=c('ENTREZID','SYMBOL'),keytype='ENSEMBL')
head(ids)
write.csv(translate_db,file = 'msu2rap.csv')
