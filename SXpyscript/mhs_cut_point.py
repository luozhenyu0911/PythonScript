from __future__ import print_function
from __future__ import division
import sys
import argparse
import pysam
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


parser = argparse.ArgumentParser(description="MHS cut point", formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('--bamfile', '-b',type= str,help="bam file",required= True)
parser.add_argument('--bedfile','-i',type= str,help="bed file",required= True)
parser.add_argument('--output','-o',type=str,help="output")
args = parser.parse_args()

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
        chrom=value[0]
        length=int(value[3])
        cutPoint1=int(value[1])
        cutPoint2=int(value[2])
        if chrom in region:
            pass
        else:
            region[chrom]={}
        if str(cutPoint1) in region[chrom]:
            region[chrom][str(cutPoint1)]+=1
        else:
            region[chrom][str(cutPoint1)]=1
        if str(cutPoint2) in region[chrom]:
            region[chrom][str(cutPoint2)]+=1
        else:
            region[chrom][str(cutPoint2)]=1
    return region


def Readscount(bedlist,readdict,rn,windows):
    sumReads=0
    allgene={}
    revwins=0-windows
    for row2 in bedlist:
        values=[]
        if str(row2[4]) is "+":
            if row2[0] in readdict:
                for k in range(int(row2[1]),int(row2[2]),windows):
                    for j in range(k,(k+windows)):
                        if str(j) in readdict[row2[0]]:
                            sumReads+=int(readdict[row2[0]][str(j)])
                        else:
                            pass
                    norReads=sumReads/(rn*windows)*1000000
                    values.append(norReads)
                    sumReads=0
            else:
                for j in range(int(row2[1]),int(row2[2]),windows):
                    values.append("0")
        else:
            if row2[0] in readdict:
                for k in range(int(row2[2]),int(row2[1]),revwins):
                    for j in range(k,(k+revwins),-1):
                        if str(j) in readdict[row2[0]]:
                            sumReads+=int(readdict[row2[0]][str(j)])
                        else:
                            pass
                    norReads=sumReads/(rn*windows)*1000000
                    values.append(norReads)
                    sumReads=0
            else:
                for j in range(int(row2[1]),int(row2[2]),windows):
                    values.append("0")
        values=map(float,values)
        allgene[str(row2[3])]=values
    return allgene
    bedfile.close()
 
def Listgene(geneID,genedict):
    Input=open(geneID,"r")
    Listmat=[]
    for line in Input:
        Id=line.rstrip()
        Listmat.append(genedict[str(Id)])
    return Listmat

def listg(g,gdict):
    Listmat=[]
    for ID in g:
        Listmat.append(gdict[str(ID)])
    return Listmat


def bamtocount(bed,Input,windowsize):
    bamfile=pysam.AlignmentFile(Input,"rb")
    region=bamInfo(bamfile)
    N=bamfile.count()
    geneall=Readscount(bed,region,N,windowsize)
    return geneall

def splist(l,n):
    lenth = len(l)
    sz = lenth // n
    c= lenth % n
    lst = []
    i = 0
    while i < n:
        if i < c:
            bs=sz+1
            lst.append(l[i*bs : i*bs+bs])
        else:
            lst.append(l[i*sz+c : i*sz + c + sz])
        i+=1
    return lst

def plotsmooth(values):
    n=len(values)
    num=range(0,n)
    num=np.array(num)
    xnew=np.linspace(num.min(),num.max(),5*n)
    func=interp1d(num,values,kind="cubic")
    ynew=func(xnew)
    plt.plot(xnew,ynew)

bamfile=pysam.AlignmentFile(args.bamfile,"rb")
bedfile=open(args.bedfile,"r")

readdict=bamInfo(bamfile)

bedpos=[]
for genepos in bedfile:
    row=genepos.rstrip().split()
    for i in range(int(row[1]),int(row[2])):
        if row[0] in readdict:
            if str(i) in readdict[row[0]]:
                readscount=readdict[row[0]][str(i)]
            else:
                readscount=0
            print(row[0],i,readscount,sep="\t",end="\n")
