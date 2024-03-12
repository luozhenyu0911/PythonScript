# -*- coding: UTF-8 -*-
from __future__ import print_function
from __future__ import division
import subprocess
import sys
import gzip

"""
This tool to convert CGmap file to bigwig file
2019.10.23
"""

for i in range(1,len(sys.argv)):
    if sys.argv[i]=='-b':  #CGmap file
        CGmap=sys.argv[i+1]
    elif sys.argv[i]=='-i': # Genome size file
        Gfile=sys.argv[i+1]
    elif sys.argv[i]=='-c': # smallest coverage
        miniCoverage=int(sys.argv[i+1])

#读取染色体文件 2020.10.16
chr = []
file = open(Gfile,'r')
for line in file:
    chr.append(line.strip().split()[0])
file.close()

wigout = CGmap.split('.')[0]+".wig"
wig = open(wigout,'w')

#遍历CGmap文件，把甲基化率写入wig文件 2020.10.16
Chr=0
with gzip.open(CGmap,'rb') as CG:
    for line in CG:
        line=line.strip().split()
        if Chr != line[0]:
            Chr = line[0]
            wig.write("variableStep chrom="+Chr+"\n")
            if int(line[7])>=miniCoverage:
                if line[1]=="C":
                    level = '%.3f' % (int(line[6])/int(line[7]))
                    wig.write(line[2]+"\t"+level+"\n")
                else:
                    level = '%.3f' % (0-int(line[6])/int(line[7]))
                    wig.write(line[2]+"\t"+level+"\n")
        else:
            if int(line[7])>=miniCoverage:
                if line[1]=="C":
                    level = '%.3f' % (int(line[6])/int(line[7]))
                    wig.write(line[2]+"\t"+level+"\n")
                else:
                    level = '%.3f' % (0-int(line[6])/int(line[7]))
                    wig.write(line[2]+"\t"+level+"\n")

#wig文件转换为bw
bw = CGmap.split('.')[0]+".bw"
subprocess.call(["wigToBigWig",wigout,Gfile,bw])
subprocess.call(["rm",wigout])
