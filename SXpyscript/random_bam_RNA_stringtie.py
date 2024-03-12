# -*- coding: UTF-8 -*-
from __future__ import print_function
from __future__ import division
import pysam
import sys 
import random
import subprocess
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import argparse

parser = argparse.ArgumentParser(description="This tool can compute gene numbers by step 10% reads from a bam and plot a smooth line to check the RNA-seq saturation level", formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('--bamlist', '-b',type= str,help="Input files contain one bam filenames for each line",required= True)
parser.add_argument('--outf','-o',type= str,help="Figure name",required= True)
parser.add_argument('--outGeneNum','-p',type=str,help="Expresed gene number file")
parser.add_argument('--Gff','-G',type=str,help="Gene annotation GFF file.",required= True)
args = parser.parse_args()

def plotsmooth(values):
    n=len(values)
    num=range(0,n)
    num=np.array(num)
    xnew=np.linspace(num.min(),num.max(),5*n)
    func=interp1d(num,values,kind="cubic")
    ynew=func(xnew)
    plt.plot(xnew,ynew)

def randomReads(bamfiles,bedfile):
    prefixall=[]
    bamprefix=bamfiles.replace(".bam","")
    for i in range(0,20):
        name=bamprefix+"_"+str(i+1)+".bam"
        prefix=bamprefix+"_"+str(i+1)
        prefixall.append(prefix)
    genenum=[]
    for i in range(0,20):
        if i < 19:
            output=prefixall[i]+".bam"
            subprocess.call(["samtools","view","-s",str(float((i+1)*0.05)),"-o",output,bamfiles])
            subprocess.call(["samtools","index",output])
            #计算每个百分比梯度取样的FPKM>1的基因数
            subprocess.call(["stringtie","-e","-p","20","-G",bedfile,"-o",prefixall[i]+".gtf","-l",prefixall[i],"-A",prefixall[i]+"_FPKM.txt",output])
            fpkmfile=open(prefixall[i]+"_FPKM.txt","r")
            next(fpkmfile)
            n=0
            for line2 in fpkmfile:
                row=line2.rstrip().split("\t")
                if float(row[7]) >= 1:
                    n+=1
            genenum.append(n)
            fpkmfile.close()
            bamindex=output+".bai"
            subprocess.call(["rm",output])
            subprocess.call(["rm",bamindex])
            subprocess.call(["rm",prefixall[i]+"_FPKM.txt"])
            subprocess.call(["rm",prefixall[i]+".gtf"])
        else:
            subprocess.call(["stringtie","-e","-p","20","-G",bedfile,"-o",prefixall[i]+".gtf","-l",prefixall[i],"-A",prefixall[i]+"_FPKM.txt",bamfiles])
            fpkmfile=open(prefixall[i]+"_FPKM.txt","r")
            next(fpkmfile)
            n=0
            for line2 in fpkmfile:
                row=line2.rstrip().split("\t")
                if float(row[7]) >= 1:
                    n+=1
            genenum.append(n)
            fpkmfile.close()
            subprocess.call(["rm",prefixall[i]+"_FPKM.txt"])
            subprocess.call(["rm",prefixall[i]+".gtf"])
    return genenum

plt.switch_backend('PDF')

listfile=open(args.bamlist,"r")
label=[]

Genebed=args.Gff
outP=open(args.outGeneNum,"w")

for sample in listfile:
    bamdir=str(sample.rstrip())
    dirall=sample.rstrip().split()
    prefix=dirall[-1].replace(",bam","")
    label.append(prefix)
    valuelist=randomReads(bamdir,Genebed)
    x=range(len(valuelist))
    print(*valuelist,sep="\t",end="\n",file=outP)
    plt.plot(x,valuelist)

xlabs=[]
xlabp=[]

for i in range(0,100,5):
    xlabs.append(i/5)
    xlabp.append(str(i+5))
plt.ylabel("Expressed gene number")
plt.xticks(xlabs,xlabp)
plt.legend(label)
plt.savefig(args.outf)
