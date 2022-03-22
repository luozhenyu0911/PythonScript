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
 2020.09.09
 python *py -b *.bam -bs 1 -o *.bw
'''
parser = argparse.ArgumentParser(description="This tool creat genome SINGLE END cut point (5' ends) bigwig from a bam.", 
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
    BAM=pysam.AlignmentFile(bam,'rb')
    region={}
    for read in BAM.fetch():
        a=[]
        a.append(read.reference_name)
        a.append(read.reference_start)
        a.append(read.reference_end)
        chrom=str(a[0])
        if chrom in region:
            pass
        else:
            region[chrom]={}
        if read.is_reverse == False:
            if str(a[1]) in region[chrom]:
                region[chrom][str(a[1])]+=1
            else:
                region[chrom][str(a[1])]=1
        elif read.is_reverse == True:
            if str(a[2]) in region[chrom]:
                region[chrom][str(a[2])]+=1
            else:
                region[chrom][str(a[2])]=1
        else:
            pass
    return region

#BAM/SAM format is 1-based,pysam是0-based
def computeBamCuts(bam,bamdict,binsize,outbigwig):
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
        bigwig.addEntries(chroms, starts, ends=ends, values=values)
    bigwig.close()

if __name__ == '__main__':
    Bamcoverage=bamTodict(bam)
    computeBamCuts(bam,Bamcoverage,bs,bw)