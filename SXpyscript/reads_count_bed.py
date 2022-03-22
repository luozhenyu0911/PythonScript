# -*- coding: UTF-8 -*-
from __future__ import print_function
from __future__ import division
import sys
import pysam

bedfile=open(sys.argv[2],"r")

region={}

def bamInfo(Input):
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
    return region

def Readscount(bed,readdict,rn):
    sumReads=0
    if bed[0] in readdict:
        for k in range(int(bed[1]),int(bed[2])):
            if str(k) in readdict[bed[0]]:
                sumReads+=int(readdict[bed[0]][str(k)])
            else:
                pass
        length=int(bed[2])-int(bed[1])
        norReads=sumReads/(rn*length)*1000000000
        values=norReads
    else:
        values=0
    return values

bedpos=[]
for genepos in bedfile:
    row=genepos.rstrip().split()
    #输入的bed四列
    [chrom,start,end,geneID]=row
    bedpos.append([chrom,start,end,geneID])

bamlist=open(sys.argv[1],"r")
allbamreads={}
for line in bamlist:
    bamdir=line.rstrip()
    chipfile=pysam.AlignmentFile(bamdir,"rb")
    N1=chipfile.count()
    chipdict=bamInfo(chipfile)
    for region in bedpos:
        chip=Readscount(region,chipdict,N1)
        if region[3] in allbamreads:
            allbamreads[region[3]].append(chip)
        else:
            allbamreads[region[3]]=[]
            allbamreads[region[3]].append(chip)
  
outputfile=open(sys.argv[4],"w")
labels=sys.argv[3].split(',')
print("GeneID",*labels,sep="\t",end="\n",file=outputfile)
for region in bedpos:
    if region[3] in allbamreads:
        print(region[3],*allbamreads[region[3]],sep="\t",end="\n",file=outputfile)  
