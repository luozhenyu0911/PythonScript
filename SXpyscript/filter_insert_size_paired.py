# -*- coding: UTF-8 -*-
from __future__ import print_function
import pysam
import argparse
import subprocess
import os

example_text = '''example:
 2020.7.3
 python **.py -b input.bam --length-below 100 -o below100.bam
 python **.py -b input.bam --length-above 100 -o above100.bam
 python **.py -b input.bam --length-between 140 180 -o between140_180.bam
'''
parser = argparse.ArgumentParser(description="This script is used to extract specific fragment size from PAIRED alignment bam file.", 
                                 formatter_class= argparse.RawTextHelpFormatter,
                                 prog='base_maker',
                                 epilog=example_text)
parser.add_argument("-b","--bam",help="Bam format file, should be sorted.",required= True)
group=parser.add_mutually_exclusive_group(required= True)
group.add_argument("--length-below",help="Paired reads with fragment (insert size) below this number.",type=int)
group.add_argument("--length-above",help="Paired reads with fragment (insert size) above this number.",type=int)
group.add_argument("--length-between",nargs='+',help="Paired reads with fragment (insert size) between an range.",type=int)
parser.add_argument("-o","--output",help="Output bam file name.",required=True)

args = parser.parse_args()
#print(args)
bam = args.bam
out = args.output

if args.bam:
    bamfile = pysam.AlignmentFile(bam,'rb')
    if args.length_below:
        with pysam.AlignmentFile(out, "wb",template=bamfile) as outf:
            for read in bamfile.fetch():
                #只要比上的reads
                if not read.is_unmapped:
                    if 0 < abs(read.template_length) <= args.length_below:
                        outf.write(read)
                    else:
                        pass
        subprocess.call(["samtools","index",out]) 

    elif args.length_above:
        with pysam.AlignmentFile(out, "wb",template=bamfile) as outf:
            for read in bamfile.fetch():
                #只要比上的reads
                if not read.is_unmapped:
                    if abs(read.template_length) > args.length_above:
                        outf.write(read)
                    else:
                        pass
        subprocess.call(["samtools","index",out])

    elif args.length_between:
        short=int(min(args.length_between))
        long=int(max(args.length_between))
        with pysam.AlignmentFile(out, "wb",template=bamfile) as outf:
            for read in bamfile.fetch():
                #只要比上的reads
                if not read.is_unmapped:
                    if short <= abs(read.template_length) <= long:
                        outf.write(read)
                    else:
                        pass
        subprocess.call(["samtools","index",out])
        
    else:
        pass
