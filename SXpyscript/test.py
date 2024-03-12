# coding=UTF-8
from __future__ import print_function
from interval import Interval, IntervalSet #Compute the union of exons of different transcripts of the same gene
import re
per_gene_exon=[]
per_len_exon = 0
exon_len = []
id = []
gene_len = []
mean_exon_num = []
i,j,k = 0,0,0
#out = open('/gss1/home/lxx20180918/out.txt','w')
f = open('/gss1/home/lxx20180918/index/rice7.gff3','r')
file = f.readlines() #将文件读成list
for index in range(len(file)): # 遍历文件每一行
    line = file[index].strip().split()  # 去掉结尾制表符
    if len(line) != 2:  #排除第一行
        if line[2] == 'gene' or index == len(file)-1: # 为基因时，保存基因id，保存上一个基因的外显子长度和个数
            id.append(''.join(re.findall(r"(?<=ID=)\w+?(?=\;)",line[8])))
            gene_len.append((int(line[4])-int(line[3])))
            i += 1
            if len(per_gene_exon) != 0: #判断外显子坐标是不是空的，如不是，计算外显子总长度
                intervals = IntervalSet([Interval.between(min, max) for min, max in per_gene_exon])
                for _ in intervals:
                    per_len_exon += _.upper_bound - _.lower_bound
                exon_len.append(int(per_len_exon)) # 将计算好的外显子长度追加到列表中
                per_len_exon=0 # 清零
                per_gene_exon = []
            if j != 0 and k != 0: # 计算每个基因的平均外显子个数
                mean_exon_num.append(round(float(k)/j,2))
                j,k = 0,0
            if index == len(file)-1: # 读到最后一行时，代码会把最后一行的基因id和坐标追加进去，所以要pop掉
                id.pop()
                gene_len.pop()
        elif line[2] == 'mRNA': # 转录本个数
            j += 1
        elif line[2] == 'exon': # 每个转录本外显子个数
            k += 1
            exon_id = ''.join(re.findall(r"(?<=ID=)\w+?(?=\.)",line[8]))
            if exon_id == id[i-1]: # 判断，将同一个基因下的外显子保存到临时变量中
                per_gene_exon.append((int(line[3]),int(line[4]))) 

for a,b,c,d in zip(id,gene_len,exon_len,mean_exon_num):
    print(a,b,c,d,sep='\t',end='\n')#,file=out)
