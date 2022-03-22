from __future__ import print_function
from __future__ import division
import sys
import argparse
from pybedtools import BedTool
import pybedtools as pybt
import pandas as pd
import numpy as np
import subprocess
import pysam

example_text = '''example:
 python **.py -b **.bam,**.bam... -g **.bed,**.bed... -c 6 -o output.txt
'''
parser = argparse.ArgumentParser(description="Calculate the FPKM of peak",
                                formatter_class= argparse.RawTextHelpFormatter,
                                epilog=example_text)
group=parser.add_mutually_exclusive_group(required= True)
group.add_argument('--bam','-b',type=str,help="<bamfiles>"+"\t"+"data files(Support for multiple files, Separated by commas)",metavar='')
group.add_argument('--bamList','-bl',type=str,help="<bamfiles>"+"\t"+"data files(Support for multiple files in one file, one raw present one file).",metavar='')

group=parser.add_mutually_exclusive_group(required= True)
group.add_argument('--bed','-g',type=str,help="<bedfiles>"+"\t"+"data files(Support for multiple files, Separated by commas)",metavar='')
group.add_argument('--bedList','-gl',type=str,help="<bedfiles>"+"\t"+"data files(Support for multiple files in one file, one raw present one file).",metavar='')

parser.add_argument('--bedcol','-c',type=str,help="How many columns are there per bed file ",required= True)
parser.add_argument('--outfile', '-o',type= str,help="output file col('RPKM','bamlabel','bedlabel')",required= True)
args = parser.parse_args()

bam=args.bam
bamList=args.bamList
bed=args.bed
bedList=args.bedList
bedcol=args.bedcol
outfile=args.outfile

if args.bam:
    bamlist=args.bam.split(",")
else:
    bamlist=[]
    with open(args.bamList,'r') as bedtxt:
        for line in bedtxt:
            line=line.strip().split()
            bamlist.append(line[0])
if args.bed:
    bedlist=args.bed.split(",")
else:
    bedlist=[]
    with open(args.bedList,'r') as bedtxt:
        for line in bedtxt:
            line=line.strip().split()
            bedlist.append(line[0])

colnames1=['Chr','start','end','count']
colnames2=['Chr','start','end','v1','v2','v3','count']
if int(args.bedcol) == 3:
    name=colnames1
elif int(args.bedcol) == 6:
    name=colnames2
df3 = pd.DataFrame(columns=['RPKM','bamlabel','bedlabel'])
df4 = pd.DataFrame(columns=['RPKM','bamlabel','bedlabel'])
for i,bed in enumerate(bedlist):
    if i <1:
        for j,bam in enumerate(bamlist):
            Input1=pysam.AlignmentFile(bam,"rb")
            total_num = int(Input1.count())/1000000000
            if j <1:
                bed_read=open(bed.strip().split('.')[0]+bam.strip().split('.')[0]+'tmp'+'.txt',"w")
                subprocess.call(["bedtools","multicov","-bams",bam,"-bed",bed],stdout=bed_read)
                bed_read.close()
                df1=pd.read_csv(bed.strip().split('.')[0]+bam.strip().split('.')[0]+'tmp'+'.txt',sep='\t',header=None,names=name)
                df1['length']=df1['end']-df1['start']
                df1['RPKM']=(df1['count']/df1['length'])/total_num
                df1['bamlabel']=bam.strip().split('.')[0]
                df1['bedlabel']=bed.strip().split('.')[0]
                df2=df1.loc[:,['RPKM','bamlabel','bedlabel']]
            else:
                bed_read=open(bed.strip().split('.')[0]+bam.strip().split('.')[0]+'tmp'+'.txt',"w")
                subprocess.call(["bedtools","multicov","-bams",bam,"-bed",bed],stdout=bed_read)
                bed_read.close()
                df1=pd.read_csv(bed.strip().split('.')[0]+bam.strip().split('.')[0]+'tmp'+'.txt',sep='\t',header=None,names=name)
                df1['length']=df1['end']-df1['start']
                df1['RPKM']=(df1['count']/df1['length'])/total_num
                df1['bamlabel']=bam.strip().split('.')[0]
                df1['bedlabel']=bed.strip().split('.')[0]
                df3=df1.loc[:,['RPKM','bamlabel','bedlabel']]
                df2=df2.append(df3)
    else:
        for k,bam in enumerate(bamlist):
            Input1=pysam.AlignmentFile(bam,"rb")
            total_num = int(Input1.count())/1000000000
            if k <1:
                bed_read=open(bed.strip().split('.')[0]+bam.strip().split('.')[0]+'tmp'+'.txt',"w")
                subprocess.call(["bedtools","multicov","-bams",bam,"-bed",bed],stdout=bed_read)
                bed_read.close()
                df1=pd.read_csv(bed.strip().split('.')[0]+bam.strip().split('.')[0]+'tmp'+'.txt',sep='\t',header=None,names=name)
                df1['length']=df1['end']-df1['start']
                df1['RPKM']=(df1['count']/df1['length'])/total_num
                df1['bamlabel']=bam.strip().split('.')[0]
                df1['bedlabel']=bed.strip().split('.')[0]
                df4=df1.loc[:,['RPKM','bamlabel','bedlabel']]
            else:
                bed_read=open(bed+bam+'tmp'+'.txt',"w")
                subprocess.call(["bedtools","multicov","-bams",bam,"-bed",bed],stdout=bed_read)
                bed_read.close()
                df1=pd.read_csv(bed+bam+'tmp'+'.txt',sep='\t',header=None,names=name)
                df1['length']=df1['end']-df1['start']
                df1['RPKM']=(df1['count']/df1['length'])/total_num
                df1['bamlabel']=bam.strip().split('.')[0]
                df1['bedlabel']=bed.strip().split('.')[0]
                df3=df1.loc[:,['RPKM','bamlabel','bedlabel']]
                df4=df4.append(df3)
    df2=df2.append(df4)
df2.to_csv(outfile,index=None,sep='\t')
