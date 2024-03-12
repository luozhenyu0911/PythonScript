# -*- coding: UTF-8 -*-
from __future__ import print_function
from __future__ import division
import argparse
import pandas as pd
import numpy as np
import pysam
from math import sqrt,ceil
#import matplotlib.pyplot as plt
#import seaborn as sns
from scipy.interpolate import make_interp_spline, BSpline
#plt.switch_backend('PDF')
#plt.rcParams['pdf.fonttype'] = 42
##@jit

example_text = '''example:
 2020.7.9
 python **.py -b **.bam,**.bam... -G **.bed -L bamlabel1,bamlabel2 -o normalized_read_counts
'''
parser = argparse.ArgumentParser(description="This tool creates mutiple bam normalized reads counts within a bed.", 
                                 formatter_class= argparse.RawTextHelpFormatter,
                                 prog='base_maker',
                                 epilog=example_text)

parser.add_argument('--bam','-b',type= str,help="<bam or sam files>"+"\t"+"data files(Support for multiple files, Separated by commas), should be sorted format (Note: Need to filter out low-quality reads first).",required= True)
parser.add_argument('--bamlabel','-L',type=str,help="<bam legend names>"+"\t"+"Labels of Bam file. (Support for multiple legends).",required= True)
parser.add_argument('--bed','-G',type=str,help="<bedfiles>",required= True)
parser.add_argument('--outcounts','-o',type=str,help="<normalized_read_counts>",required= True)
args = parser.parse_args()

bamlist=args.bam.split(",")
Bamlabels=args.bamlabel.split(",")
Bed=args.bed
outf=args.outcounts


def bamTodict(Input):
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
            readInfo=[]
            readInfo.append(str(read.reference_name))
            readInfo.append(int(read.reference_start))
            readInfo.append(int(read.reference_end))
            readInfo.append(int(readInfo[2])-int(readInfo[1])+1)
            ReadInfo[ID]=readInfo
    region={}
    for value in ReadInfo.values():
        start=value[1]
        end=value[2]
        length=value[3]
        chrom=value[0]
        remain=length%2
        #偶数长度的中点有两个，长度n-1的中点和长度n+1的中点
        if int(remain) == 0:
            centerPoint1=int((start+end-1)/2)
            centerPoint2=int((start+end+1)/2)
            if chrom in region:
                pass
            else:
                region[chrom]={}
            if str(centerPoint1) in region[chrom]:
                region[chrom][str(centerPoint1)]+=0.5
            else:
                region[chrom][str(centerPoint1)]=0.5
            if str(centerPoint2) in region[chrom]:
                region[chrom][str(centerPoint2)]+=0.5
            else:
                region[chrom][str(centerPoint2)]=0.5
        #奇数长度的中点等于首尾之和的一半
        else:
            centerPoint=int((start+end)/2)
            if chrom in region:
                pass
            else:
                region[chrom]={}
            if str(centerPoint) in region[chrom]:
                region[chrom][str(centerPoint)]+=1
            else:
                region[chrom][str(centerPoint)]=1
    return region

def NormalizedReadCounts(readdict,bed,totalreads):
    Bedf=open(bed,'r')
    norReadslist=[]
    for line in Bedf:
        line=line.strip().split()
        chr=line[0]
        start=int(line[1])
        end=int(line[2])
        length=end-start+1
        counts=0
        norReads=0
        if chr in readdict:
            for i in range(start,end):
                if str(i) in readdict[chr]:
                    counts+=readdict[chr][str(i)]
                else:
                    pass
            norReads=(counts/(length*totalreads))*1000000000
            norReadslist.append(norReads)
        else:
            norReadslist.append(0)
    Bedf.close()
    return norReadslist

countsDataframe=pd.DataFrame()
for Bam in range(len(bamlist)):
    bamf=pysam.AlignmentFile(bamlist[Bam],'rb')
    bamTotalReads=bamf.count()
    bamrCoverageDict=bamTodict(bamf)
    values=NormalizedReadCounts(bamrCoverageDict,Bed,bamTotalReads)
    countsDataframe[Bamlabels[Bam]]=values
countsDataframe.to_csv(outf,sep='\t',index=True)
