from __future__ import print_function
from __future__ import division
import sys
import numpy as np


def bamInfo(Input):
    region={}
    for read in Input:
        col=read.rstrip().split()
        if int(col[7]) >= 5:
            if str(col[0]) in region:
                region[str(col[0])][str(col[1])]=(int(col[6]),int(col[7]))
            else:
                region[str(col[0])]={}
                region[str(col[0])][str(col[1])]=(int(col[6]),int(col[7]))
    return region


def Readscount(bedlist,readdict):
    metC=0
    sumC=0
    allgene=[]
    for row2 in bedlist:
        values=[]
        if row2[0] in readdict:
            for k in range(int(row2[1]),int(row2[2])):
                if str(k) in readdict[row2[0]]:
                    sumC+=readdict[row2[0]][str(k)][1]
                    metC+=readdict[row2[0]][str(k)][0]
                else:
                    pass
            if sumC == 0:
                methyratio=0
            else:
                methyratio=metC/sumC
            values.append(methyratio)
        else:
            pass
        metC=0
        sumC=0

        if row2[0] in readdict:
            for x in range(int(row2[2]),int(row2[3])):
                if str(x) in readdict[row2[0]]:
                    sumC+=readdict[row2[0]][str(x)][1]
                    metC+=readdict[row2[0]][str(x)][0]
                else:
                    pass
            if sumC == 0:
                methyratio=0
            else:
                methyratio=metC/sumC
            values.append(methyratio)
        else:
            pass
        metC=0
        sumC=0

        if row2[0] in readdict:
            for k in range(int(row2[3]),int(row2[4])):
                if str(k) in readdict[row2[0]]:
                    sumC+=readdict[row2[0]][str(k)][1]
                    metC+=readdict[row2[0]][str(k)][0]
                else:
                    pass
            if sumC == 0:
                methyratio=0
            else:
                methyratio=metC/sumC
            values.append(methyratio)
        else:
            pass
        metC=0
        sumC=0
        allgene.append(values)
    return allgene
    bedfile.close()
 

def bamtocount(bed,Input):
    bamfile=open(Input,"r")
    region=bamInfo(bamfile)
    geneall=Readscount(bed,region)
    return geneall

bedlist=sys.argv[2].split(",")
upstream = int(sys.argv[3])
downstream=int(sys.argv[4])

for bedIn in bedlist:
    legend=bedIn.split(".")[0]
    bedfile=open(bedIn,"r")
    bedpos=[]
    for genepos in bedfile:
        row=genepos.rstrip().split()
        [chrom,start,end,geneID,score,direct]=row
        startnew=int(start)-int(upstream)
        endnew=int(end)+int(downstream)
        bedpos.append([chrom,startnew,start,end,endnew,geneID,direct])
    genevalues=bamtocount(bedpos,sys.argv[1])
    for i in genevalues:
        print(legend,i[0],i[1],i[2],sep="\t",end="\n")
    bedfile.close()
    
    
    
    
