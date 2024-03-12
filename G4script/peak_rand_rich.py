# -*- coding: UTF-8 -*-
'''
python this.py bedlist dmrlist genome outfile
主要用来就算多个peak与重要bed文件做bedtools交集取得bed中交的数目1，在比上随机peak中bed的数目2，求foldchange（数目1/数目2），并求得sd
'''
from __future__ import print_function
from __future__ import division
import sys
import linecache
import subprocess
# import scipy.stats
import numpy as np
import pandas as pd
from pybedtools import BedTool
bed=open(sys.argv[1],'r')
dmrF=open(sys.argv[2],'r')
genome=sys.argv[3]
outf=open(sys.argv[4],'w')
bedlist=[]
for i in bed:
    bedlist.append(i.strip())
    
dmrlist=[]
for i in dmrF:
    dmrlist.append(i.strip())

def get_file_row(file1):
    count = 0
    thefile = open(file1, 'rb')
    while True:
        buffer = thefile.read(8192*1024)
        if not buffer:
            break
        count += buffer.count('\n')
    thefile.close( )
    return count
print("sample",'mean','sd','pvalue',sep='\t',file=outf)
values_list={}
for bedf in bedlist:
    bedf=bedf.strip()
    prefix=bedf.strip().split('.')[0]
    InBed=BedTool(bedf)
    for dmr in dmrlist:
        dmr=dmr.strip()
        values=[]
        valuesP=[]
        dmr_prefix=dmr.strip().split('.')[0]
        jiaobed=open(prefix+'jiao'+'.bed',"w")
        subprocess.call(["bedtools","intersect","-a",dmr,"-b",bedf,"-wa",'-u'],stdout=jiaobed)
        observed=get_file_row(prefix+'jiao'+'.bed')
        for i in range(0,10):
            InBed.shuffle(g=sys.argv[3],noOverlapping=True).moveto(str(prefix)+"R"+str(i)+".bed")
            jiaobedr=open(prefix+dmr_prefix+'jiao'+str(i)+'.bed',"w")
            subprocess.call(["bedtools","intersect","-a",dmr,"-b",str(prefix)+"R"+str(i)+".bed","-wa",'-u'],stdout=jiaobedr)
            expectation=get_file_row(prefix+dmr_prefix+'jiao'+str(i)+'.bed')
            fold_change=observed/expectation
            values.append(fold_change)
            subprocess.call(["rm",prefix+dmr_prefix+'jiao'+str(i)+'.bed'])
            subprocess.call(["rm",str(prefix)+"R"+str(i)+".bed"])
            # subprocess.call(["rm",prefix+'jiao'+'.bed'])
        # tmp=[1]
        # values_R=tmp*10
        # pvalue=scipy.stats.ranksums(values,values_R).pvalue
        arr_mean = np.mean(values)
        arr_std = np.std(values)
        valuesP.extend((arr_mean,arr_std))
        values_list[prefix+dmr_prefix]=valuesP
        
for K in values_list:
    print(K,*values_list[K],sep='\t',file=outf)
    
bed.close()
dmrF.close()
# genome=sys.argv[3]
outf.close()
