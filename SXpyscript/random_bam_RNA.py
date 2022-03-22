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
parser.add_argument('--bedfile','-G',type=str,help="Gene bed file.",required= True)
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
    for i in range(0,10):
        name=bamprefix+"_"+str(i+1)+".bam"
        prefix=bamprefix+"_"+str(i+1)
        prefixall.append(prefix)
    genenum=[]
    for i in range(0,10):
        if i < 9:
            output=prefixall[i]+".bam"
            subprocess.call(["samtools","view","-s",str(float((i+1)*0.1)),"-o",output,bamfiles])
            subprocess.call(["samtools","index",output])
            #计算每个百分比梯度取样的能检测到reads的基因个数
            countout=open("tmp.count","w")
            subprocess.call(["bedtools","intersect","-a",bedfile,"-b",output,"-c"],stdout=countout)
            countout.close()
            countfile=open("tmp.count","r")
            n=0
            for line2 in countfile:
                row=line2.rstrip().split()
                if int(row[6]) > 0:
                    n+=1
            genenum.append(n)

            bamindex=output+".bai"
            subprocess.call(["rm",output])
            subprocess.call(["rm",bamindex])
            subprocess.call(["rm","tmp.count"])
        else:
            countout=open("tmp.count","w")
            subprocess.call(["bedtools","intersect","-a",bedfile,"-b",bamfiles,"-c"],stdout=countout)
            countout.close()
            countfile=open("tmp.count","r")
            n=0
            for line2 in countfile:
                row=line2.rstrip().split()
                if int(row[6]) > 0:
                    n+=1
            genenum.append(n)
            subprocess.call(["rm","tmp.count"])
    return genenum

plt.switch_backend('PDF')

listfile=open(args.bamlist,"r")
label=[]

Genebed=args.bedfile
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

for i in range(0,100,10):
    xlabs.append(i/10)
    xlabp.append(str(i+10))
plt.ylabel("Expressed gene number")
plt.xticks(xlabs,xlabp)
plt.legend(label)
plt.savefig(args.outf)
