'''
python tmp.py bamlist D12h_TPM_high.bed D12h_TPM_high.txt zscore.txt
'''
from __future__ import print_function 
from __future__ import division
from pybedtools import BedTool
from scipy.stats import zscore
import pybedtools as pybt
import sys
import pandas as pd
import numpy as np
import pysam
import subprocess


def readdict(Input):
    ReadInfo={}
    for read in Input.fetch():
        ID=read.query_name
        if ID in ReadInfo:
            if int(read.reference_start) < int(ReadInfo[ID][1]):
                ReadInfo[ID][1]=int(read.reference_start)
                ReadInfo[ID][3]=int(ReadInfo[ID][2])-int(ReadInfo[ID][1])+1
            else:
                ReadInfo[ID][2]=int(read.reference_end)
                ReadInfo[ID][3]=int(ReadInfo[ID][2])-int(ReadInfo[ID][1])+1
        else:
            a=[]
            a.append(str(read.reference_name))
            a.append(int(read.reference_start))
            a.append(int(read.reference_end))
            a.append(int(a[2])-int(a[1])+1)
            ReadInfo[ID]=a
    region={}
    for value in ReadInfo.values():
        length=value[3]
        chrom=value[0]
        remain=length%2
        if int(remain) == 0:
            centerPoint=int((int(value[1]) + int(value[2]))/2)
            if chrom in region:
                pass
            else:
                region[chrom]={}
            if str(centerPoint) in region[chrom]:
                region[chrom][str(centerPoint)]+=2
            else:
                region[chrom][str(centerPoint)]=2
        else:
            centerPoint1=int((int(value[1]) + int(value[2])+1)/2)
            centerPoint2=int((int(value[1]) + int(value[2])-1)/2)
            if chrom in region:
                pass
            else:
                region[chrom]={}
            if str(centerPoint1) in region[chrom]:
                region[chrom][str(centerPoint1)]+=1
            else:
                region[chrom][str(centerPoint1)]=1
            if str(centerPoint2) in region[chrom]:
                region[chrom][str(centerPoint2)]+=1
            else:
                region[chrom][str(centerPoint2)]=1

    num = (int(Input.count()))/1000000
    return region,num

def bed_get(bed1,region,num,bamf1):
    sumReads=0
    obvalue = []
    obvalue.append(str(bamf1))
    # exvalue = []
    for line2 in bed1:
        col=line2.rstrip().split()
        start=int(col[1])
        end=int(col[2])
        length=end-start+1
        if col[0] in region:
            for i in range(start,end+1):
                if str(i) in region[col[0]]:
                    sumReads+=int(region[col[0]][str(i)])
                else:
                    pass
        bizhi = float(sumReads/(length*num))
        obvalue.append(bizhi)
        sumReads=0
    return obvalue

bamf=open(sys.argv[1],'r')
bamlist=[]
for i in bamf:
    bamlist.append(i.strip())
outf1=sys.argv[3]
outf2=sys.argv[4]

value_list=[]
for bamf1 in bamlist:
    bed2=open(sys.argv[2],'r')
    bamf1=bamf1.strip()
    Input1=pysam.AlignmentFile(bamf1,"rb")
    region,num=readdict(Input1)
    obvalue=bed_get(bed2,region,num,bamf1)
    value_list.append(obvalue)
    bed2.close()

df1=pd.DataFrame(value_list).T
df1.to_csv(outf1,sep='\t',header=None,index=None)
f2=pd.read_csv(outf1,sep='\t')
df3=f2.T
df3=df3.apply(zscore).T
df3.to_csv(outf2,sep='\t',index=None)
bamf.close()



