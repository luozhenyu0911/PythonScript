# -*- coding: UTF-8 -*-
from __future__ import print_function
from __future__ import division
#from pyfasta import Fasta
import pyBigWig
import pysam
import argparse
#import pandas as pd
import numpy as np

example_text = '''example:
 2020.08.20
 python *py -b *.bam -bs 1 -o *.bw
'''
parser = argparse.ArgumentParser(description="This tool creat genome PAIRED END cut point (5' ends and 3 'ends) bigwig from a bam.", 
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
            else:
                ReadInfo[ID][2]=int(read.reference_end)
        else:
            a=[]
            a.append(str(read.reference_name)) #染色体
            a.append(int(read.reference_start)) #左切点
            a.append(int(read.reference_end)) # 右切点
            #a.append(int(a[2])-int(a[1])+1)
            ReadInfo[ID]=a
    region={}
    for value in ReadInfo.values():
        chrom=value[0]
        start=int(value[1])
        end=int(value[2])
        if chrom in region:
            pass
        else:
            region[chrom]={}
        if str(start) in region[chrom]:
            region[chrom][str(start)]+=1
        else:
            region[chrom][str(start)]=1
        if str(end) in region[chrom]:
            region[chrom][str(end)]+=1
        else:
            region[chrom][str(end)]=1
    return region

#BAM/SAM format is 1-based
def computeBamCuts(bam,bamdict,binsize,outbigwig):
    chroms=[]
    starts=[]
    ends=[]
    values=[]
    bigwig=pyBigWig.open(outbigwig,'w') #创建一个bw
    BAM=pysam.AlignmentFile(bam,'rb')
    total=BAM.count()
    Chrlist=list(BAM.references) #str
    Chrlenlist=list(BAM.lengths) #int
    #向bw中添加header，必须要一次性添加才不报错
    header=[(key,chrlen) for key,chrlen in zip(Chrlist,Chrlenlist)]
    bigwig.addHeader(header)
    for key,chrlen in zip(Chrlist,Chrlenlist):
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
                        onebased=k+1
                        if str(onebased) in bamdict[key]:
                            sumCuts+=int(bamdict[key][str(onebased)])
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
                        onebased=k+1
                        if str(onebased) in bamdict[key]:
                            sumCuts+=int(bamdict[key][str(onebased)])
                        else:
                            pass
                    normCuts=(sumCuts*1000000000)/(total*(end-start))
                    values.append(normCuts)
                else:
                    values.append(0.0)
    bigwig.addEntries(chroms, starts, ends=ends, values=values)
    bigwig.close()

if __name__ == '__main__':
    Bamcoverage=bamTodict(bam)
    computeBamCuts(bam,Bamcoverage,bs,bw)