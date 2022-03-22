# -*- coding: UTF-8 -*-
from __future__ import print_function
from __future__ import division
from scipy.signal import savgol_filter
#from pyfasta import Fasta
import pyBigWig
import pysam
import argparse
#import pandas as pd
import numpy as np

example_text = '''example:
 2020.08.25
 python *py -b *.bam -bs 1 -o *.bw
'''
parser = argparse.ArgumentParser(description="This tool creat genome PAIRED END coverage bigwig file using bam file by counting 60 bp nucleotides in the center of each paired fragment.",
                                 formatter_class= argparse.RawTextHelpFormatter,
                                 prog='base_maker',
                                 epilog=example_text)
parser.add_argument('--bam','-b',type= str,help="<bam or sam files>"+"\t"+"data files, should be sorted format (Note: Need to filter out low-quality reads first).",required= True)
parser.add_argument('--bs',type= int,help="<bin size>",required= True)
parser.add_argument('--output','-o',type=str,help="<output bigwig file>",required= True)
args = parser.parse_args()

bam=args.bam
bs=args.bs
bw=args.output

def bamTodict(bam):
    ReadInfo={}
    BAM=pysam.AlignmentFile(bam,'rb')
    for read in BAM.fetch():
        ID=read.query_name
        if ID in ReadInfo:
            if int(read.reference_start) < int(ReadInfo[ID][1]):
                ReadInfo[ID][1]=int(read.reference_start)
                ReadInfo[ID][3]=int(ReadInfo[ID][2])-int(ReadInfo[ID][1])
            else:
                ReadInfo[ID][2]=int(read.reference_end)
                ReadInfo[ID][3]=int(ReadInfo[ID][2])-int(ReadInfo[ID][1])
        else:
            readInfo=[]
            readInfo.append(str(read.reference_name))
            readInfo.append(int(read.reference_start))
            readInfo.append(int(read.reference_end))
            readInfo.append(int(readInfo[2])-int(readInfo[1]))
            ReadInfo[ID]=readInfo
    region={}
    for value in ReadInfo.values():
        start=value[1]
        end=value[2]
        length=value[3]
        chrom=value[0]
        remain=length%2
        centerlen=60 #要保留的中心长度
        #偶数长度的中点有两个，长度n-1的中点和长度n+1的中点
        if int(remain) == 0:
            #sam是1-based,pysam是0-based
            #使用中心60bp，每bp加1,共60次
            newstart=start + (length//2) - (centerlen//2)
            newend=newstart + centerlen
            if chrom in region:
                pass
            else:
                region[chrom]={}
            for base in range(newstart,newend):
                if str(base) in region[chrom]:
                    region[chrom][str(base)]+=1
                else:
                    region[chrom][str(base)]=1
        #奇数长度的中点等于首尾之和的一半
        else:
            if chrom in region:
                pass
            else:
                region[chrom]={}
            #使用中心61bp,中间59加1，两端各加0.5，共61次
            newstart=start + (length//2) - (centerlen//2) + 1
            newend=newstart + centerlen - 1
            #中心59bp
            for base in range(newstart,newend):                
                if str(base) in region[chrom]:
                    region[chrom][str(base)]+=1
                else:
                    region[chrom][str(base)]=1
            #左右两端
            if str(newstart-1) in region[chrom]:
                region[chrom][str(newstart-1)]+=0.5
            else:
                region[chrom][str(newstart-1)]=0.5
            if str(newend) in region[chrom]:
                region[chrom][str(newend)]+=0.5
            else:
                region[chrom][str(newend)]=0.5
    return region

#BAM/SAM format is 1-based,pysam是0-based
def computeBamCenter(bam,bamdict,binsize,outbigwig):
    bigwig=pyBigWig.open(outbigwig,'w') #创建一个bw
    BAM=pysam.AlignmentFile(bam,'rb')
    total=BAM.count()
    Chrlist=list(BAM.references) #str
    Chrlenlist=list(BAM.lengths) #int
    #向bw中添加header，必须要一次性添加才不报错
    header=[(key,chrlen) for key,chrlen in zip(Chrlist,Chrlenlist)]
    bigwig.addHeader(header)
    for key,chrlen in zip(Chrlist,Chrlenlist):
        chroms=[]
        starts=[]
        ends=[]
        values=[]
        #将每条染色体均分等长窗口，并遍历每个窗口
        for start in range(0,chrlen,binsize):
            sumCuts=0
            if start+binsize<chrlen:
                chroms.append(key)
                starts.append(start)
                end=start+binsize
                ends.append(end)
                if key in bamdict:
                    #计算每个窗口切点数
                    for k in range(start,end):
                        if str(k) in bamdict[key]:
                            sumCuts+=float(bamdict[key][str(k)])
                        else:
                            pass
                    normCuts=(sumCuts*1000000000)/(total*binsize)
                    values.append(normCuts)
                else:
                    values.append(0.0) #保持数据类型一致
            else:                
                chroms.append(key)
                starts.append(start)
                end=chrlen
                ends.append(end)
                if key in bamdict:
                    for k in range(start,end):
                        if str(k) in bamdict[key]:
                            sumCuts+=float(bamdict[key][str(k)])
                        else:
                            pass
                    normCuts=(sumCuts*1000000000)/(total*(end-start))
                    values.append(normCuts)
                else:
                    values.append(0.0)
        #Savitzky-Golay卷积平滑算法，smooth reads
        values_fit = savgol_filter(values, 7, 5)
        bigwig.addEntries(chroms, starts, ends=ends, values=values_fit)
    bigwig.close()

if __name__ == '__main__':
    Bamcoverage=bamTodict(bam)
    computeBamCenter(bam,Bamcoverage,bs,bw)
