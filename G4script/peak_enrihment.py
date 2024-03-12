# -*- coding: UTF-8 -*-
'''
python this.py A_bedlist1 B_dmrlist genome outfile
看bedlist2
算B在A中的数目，再看B随机在A中的数目，然后用fisher做差异分析
'''
from __future__ import print_function
from __future__ import division
import sys
import linecache
import subprocess
import scipy.stats
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

values_list={}
for bedf in bedlist:
    bedf=bedf.strip()
    prefix=bedf.strip().split('.')[0]
    InBed=BedTool(bedf)
    for dmr in dmrlist:
        dmr=dmr.strip()
        observed_all=get_file_row(dmr)
        values=[]
        valuesP=[]
        dmr_prefix=dmr.strip().split('\t')[0]
        jiaobed=open(prefix+'jiao'+'.bed',"w")
        subprocess.call(["bedtools","intersect","-a",dmr,"-b",bedf,"-wa",'-u'],stdout=jiaobed)
        observed=get_file_row(prefix+'jiao'+'.bed')
        DEG=get_file_row(dmr)
        for i in range(0,10):
            InBed.shuffle(g=sys.argv[3],noOverlapping=True).moveto(str(prefix)+"R"+str(i)+".bed")
            jiaobedr=open(prefix+dmr_prefix+'jiao'+str(i)+'.bed',"w")
            subprocess.call(["bedtools","intersect","-a",dmr,"-b",str(prefix)+"R"+str(i)+".bed","-wa",'-u'],stdout=jiaobedr)
            expectation=get_file_row(prefix+dmr_prefix+'jiao'+str(i)+'.bed')
            values.append(expectation)
            subprocess.call(["rm",prefix+dmr_prefix+'jiao'+str(i)+'.bed'])
            subprocess.call(["rm",str(prefix)+"R"+str(i)+".bed"])
        arr_mean = int(np.mean(values))
        pvalue=scipy.stats.fisher_exact([[observed,observed_all-observed],[arr_mean,observed_all-arr_mean]])[1]
        valuesP.extend((observed_all,observed,observed_all-observed,arr_mean,observed_all-arr_mean,pvalue))
        values_list[prefix+dmr_prefix]=valuesP
print("sample",'all','observed','unobserved','arr_mean','unarr_mean','pvalue',sep='\t',file=outf)
for K in values_list:
    print(K,*values_list[K],sep='\t',file=outf)
    
bed.close()
dmrF.close()
outf.close()
