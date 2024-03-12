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
 2021.1.2
 python **.py -b **.bam -G **.bed -L bamlabel -o output_filename
'''
parser = argparse.ArgumentParser(description="This tool caclulate Peak length versus average fragment size within peak. ", 
                                 formatter_class= argparse.RawTextHelpFormatter,
                                 prog='base_maker',
                                 epilog=example_text)

parser.add_argument('--bam','-b',type= str,help="<bam or sam files>"+"\t"+"data files, should be sorted format (Note: Need to filter out low-quality reads first).",required= True)
parser.add_argument('--bamlabel','-L',type=str,help="<bam legend names>"+"\t"+"Labels of Bam file. ",required= True)
parser.add_argument('--bed','-G',type=str,help="<bedfiles>",required= True)
parser.add_argument('--outf','-o',type=str,help="<out filename>",required= True)
args = parser.parse_args()

bamlist=args.bam.split(",")
Bamlabels=args.bamlabel.split(",")
Bed=args.bed
outf=args.outf

bamf=pysam.AlignmentFile(bamlist[0],'rb')
out=open(outf,'w')
print("Chr","start","end","width","median_fragment_size","average_fragment_size","label",sep='\t',file=out)
bedf=open(Bed,'r')
for line in bedf:
    line=line.strip().split()
    chr=line[0]
    start=int(line[1])
    end=int(line[2])
    length=end-start
    insertSize = []
    for read in bamf.fetch(chr,start,end):
        insertSize.append(abs(read.template_length))
    med=np.median(insertSize)
    avg=np.mean(insertSize)
    print(chr,start,end,length,med,avg,Bamlabels[0],sep='\t',file=out)
bedf.close()
out.close()
