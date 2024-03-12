# -*- coding: UTF-8 -*-
from __future__ import print_function
import pysam
import random
import os
import argparse
import subprocess

#功能
#给定一个bam,随机抽取指定reads数
parser = argparse.ArgumentParser(description="This tool can extract fixed number of reads from a bam.", formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('--bam', '-b',type= str,help="Bam file, should be sorted.",required= True)
parser.add_argument('--type', '-t',type= str,help="Pair end or single end. <PE/SE>",required= True)
parser.add_argument('--num','-n',type= float,help="MILLION of reads, should less than total reads number. "
                    "Note: The number of read PAIRS should be used if it is a paired end bam file. ",required= True)
parser.add_argument('--out','-o',type=str,help="Subsampled bam file name.",required= True)
args = parser.parse_args()

if args.bam:
    print(args.bam)
    allreads=[]
    bamfile = pysam.AlignmentFile(args.bam,'rb')
    count=bamfile.count()
    if args.type == "PE":
        fragmentNumber = count//2
        #约等于一半，实际还有multi alignment
    else:
        fragmentNumber = count
    if args.num*1000000 > fragmentNumber:
        print("The reads number you given are bigger than total reads number "+str(fragmentNumber)+"!")
    else:
        with pysam.AlignmentFile("temp_unsorted.bam", "wb", template=bamfile) as outbam:
            allreads={}
            for reads in bamfile.fetch():
                qname=reads.query_name
                if qname in allreads:
                    allreads[qname].append(reads)
                else:
                    allreads[qname]=[reads]
            rand=random.sample(allreads.keys(),int(args.num*1000000))
            for qname in rand:
                for read in allreads[qname]:
                    outbam.write(read)
        pysam.sort("-o",args.out,"temp_unsorted.bam")
        subprocess.call(["samtools","index",args.out])
        subprocess.call(["rm","temp_unsorted.bam"])