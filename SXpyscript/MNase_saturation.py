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

parser = argparse.ArgumentParser(description="This tool can compute peak numbers by step 5% reads from a bam and plot a smooth line for the peak numbers to check the saturation level", formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('--bamlist', '-b',type= str,help="Input files contain one bam filenames for each line",required= True)
parser.add_argument('--outf','-o',type= str,help="Figure name",required= True)
parser.add_argument('--outPeakNum','-p',type=str,help="Peak number file")
parser.add_argument('--tempDir','-d',type=str,help="temp file directory.",required= True)
args = parser.parse_args()

def plotsmooth(values):
    n=len(values)
    num=range(0,n)
    num=np.array(num)
    xnew=np.linspace(num.min(),num.max(),5*n)
    func=interp1d(num,values,kind="cubic")
    ynew=func(xnew)
    plt.plot(xnew,ynew)

def randomPeak(bamfiles):
    #output=[]
    prefixall=[]
    bamprefix=bamfiles.replace(".bam","")
    for i in range(0,20):
        name=bamprefix+"_"+str(i+1)+".bam"
        prefix=bamprefix+"_"+str(i+1)
        #output.append(name)
        prefixall.append(prefix)
    peaknum=[]
    for i in range(0,20):
        if i < 19:
            N=[]
            for j in range(3):
                output=prefixall[i]+"_rand"+str(j+1)+".bam"
                #samtools view -s FLOAT 整数部分表示随机种子，小数部分表示抽取比例
                subprocess.call(["samtools","view","-s",str(j+float((i+1)*0.05)),"-o",output,bamfiles])
                subprocess.call(["samtools","index",output])
                subprocess.call(["python","/gss1/home/lxx20180918/software/danpos2/danpos.py","dpos",output,"-o",args.tempDir,"--paired","1","-p","1e-3","-q","0"])
                inputname="./saturation_check/pooled/"+prefixall[i]+"_rand"+str(j+1)+".smooth.positions.xls"
                hotspot=open(inputname,"r")
                n=len(hotspot.readlines())
                hotspot.close()
                bamindex=output+".bai"
                #subprocess.call(["rm",inputname])
                #subprocess.call(["rm",peakfile])
                subprocess.call(["rm",output])
                subprocess.call(["rm",bamindex])
                N.append(n)
            peaknum.append(int(np.mean(N)))
        else:
            #subprocess.call(["samtools","index",bamfiles])
            subprocess.call(["python","/gss1/home/lxx20180918/software/danpos2/danpos.py","dpos",bamfiles,"-o",args.tempDir,"--paired","1","-p","1e-3","-q","0"])
            inputname="./saturation_check/pooled/"+bamprefix+".smooth.positions.xls"
            hotspot=open(inputname,"r")
            n=len(hotspot.readlines())
            hotspot.close()
            peaknum.append(n)
    return peaknum

plt.switch_backend('PDF')

listfile=open(args.bamlist,"r")
label=[]

outP=open(args.outPeakNum,"w")

for sample in listfile:
    bamdir=str(sample.rstrip())
    dirall=sample.rstrip().split()
    prefix=dirall[-1].replace(",bam","")
    label.append(prefix)
    valuelist=randomPeak(bamdir)
    x=range(len(valuelist))
    print(*valuelist,sep="\t",end="\n",file=outP)
    plt.plot(x,valuelist)

xlabs=[]
xlabp=[]

for i in range(0,100,5):
    xlabs.append(i/5)
    xlabp.append(str(i+5))
plt.ylabel("nucleosome number")
plt.xticks(xlabs,xlabp)
plt.legend(label)
plt.savefig(args.outf)
