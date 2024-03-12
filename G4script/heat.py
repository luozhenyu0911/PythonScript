from __future__ import print_function 
from __future__ import division
from pybedtools import BedTool
import pybedtools as pybt
import sys
import pysam
import numpy
import subprocess
import scipy
import scipy.stats

Input=pysam.AlignmentFile(sys.argv[1],"rb")

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
Input.close()
sumReads=0

bed1=open(sys.argv[2],"r")

obvalue = []
exvalue = []
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
bed1.close()
tmpplist=[]
name = sys.argv[1].split(".")[0]
InBed=BedTool(sys.argv[2])
for i in range(0,10):
    extmpv = []
    InBed.shuffle(g=sys.argv[3],noOverlapping=True).moveto(name+"tmp.bed")
    randomb=open(name+"tmp.bed","r")
    sumReads=0
    for line2 in randomb:
        col=line2.rstrip().split("\t")
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
        extmpv.append(bizhi)
        sumReads=0
    extm=numpy.array(extmpv)
    tmpmean=extm.sum()/len(extmpv)    
    exvalue.append(tmpmean)
    tmpplist.append(extmpv)
    randomb.close()
    subprocess.call(["rm",name+"tmp.bed"])


ob=numpy.array(obvalue)
ex=numpy.array(exvalue)
mean1=ob.sum()/len(obvalue)
mean2=ex.sum()/len(exvalue)
print(sys.argv[1],mean1/mean2,sep='\t')
aa=numpy.array(tmpplist)
average1=list(aa.mean(axis=0))
print(scipy.stats.mannwhitneyu(obvalue,average1))
print(scipy.stats.ranksums(obvalue,average1))
