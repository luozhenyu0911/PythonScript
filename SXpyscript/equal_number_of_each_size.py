# -*- coding: UTF-8 -*-
from __future__ import print_function
import pysam
import random
import string
import os
import argparse
import subprocess

#功能
#给定一个bam,对于每个长度集合里的reads,抽取相等数量的reads
parser = argparse.ArgumentParser(description="This tool can extract reads of equal number for each size from a PE bam.", formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('--bam', '-b',type= str,help="Bam file, should be sorted.",required= True)
parser.add_argument('--num','-n',type= int,help="Number of reads, should less than number of shortest fragment. ",required= True)
parser.add_argument('--out','-o',type=str,help="Subsampled bam file name.",required= True)
args = parser.parse_args()

if args.bam:
    print(args.bam)
    allreads=[]
    bamfile = pysam.AlignmentFile(args.bam,'rb')
    rand_str = "".join(random.sample(string.ascii_letters,8))
    with pysam.AlignmentFile("temp_"+rand_str+"_unsorted.bam", "wb", template=bamfile) as outbam:
        allreads={}
        for reads in bamfile.fetch():
            isize = str(abs(reads.template_length))
            #为每个长度集合创建一个字典
            if isize in allreads:
                pass
            else:
                allreads[isize]={}
            #将reads存入长度集合下
            qname=reads.query_name
            if qname in allreads[isize]:
                allreads[isize][qname].append(reads)
            else:
                allreads[isize][qname]=[reads]
        #从每个长度集合中随机抽取read name
        for isize in allreads.keys():
            rand=random.sample(allreads[isize].keys(),int(args.num/2))
            for qname in rand:
                for read in allreads[isize][qname]:
                    outbam.write(read)
    pysam.sort("-o",args.out,"temp_"+rand_str+"_unsorted.bam")
    subprocess.call(["samtools","index",args.out])
    subprocess.call(["rm","temp_"+rand_str+"_unsorted.bam"])