# -*- coding: UTF-8 -*-
from __future__ import print_function
from __future__ import division
import subprocess
import sys
import gzip


for i in range(1,len(sys.argv)):
    if sys.argv[i]=='-b':  #CGmap file
        CGmap=sys.argv[i+1]
    elif sys.argv[i]=='-i': # Genome size file
        Gfile=sys.argv[i+1]

chr = []
file = open(Gfile,'r')
for line in file:
    chr.append(line.strip().split()[0])
file.close()

wigout = CGmap.split('.')[0]+".wig"
wig = open(wigout,'w')

Chr=0
with gzip.open(CGmap,'rb') as CG:
    for line in CG:
        line=line.strip().split()
        if Chr != line[0]:
            Chr = line[0]
            wig.write("variableStep chrom="+Chr+"\n")
            level = '%.3f' % (float(line[5]))
            wig.write(line[2]+"\t"+level+"\n")
        else:
            pass

#wig文件转换为bw
bw = CGmap.split('.')[0]+".bw"
subprocess.call(["wigToBigWig",wigout,Gfile,bw])
subprocess.call(["rm",wigout])
