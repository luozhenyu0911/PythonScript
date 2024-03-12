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

parser = argparse.ArgumentParser(description="This tool can compute peak numbers by step 10% reads from a bam and plot a smooth line for the peak numbers to check the saturation level", formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('--bamlist', '-b',type= str,help="Input files contain one bam filenames for each line",required= True)
#parser.add_argument('--GenomeSize', '-g',type=float,help="Genome size for peak calling",required=True)
parser.add_argument('--outf','-o',type= str,help="Figure name",required= True)
parser.add_argument('--outPeakNum','-p',type=str,help="Peak number file")
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
    output=[]
    prefixall=[]
    bamprefix=bamfiles.replace(".bam","")
    for i in range(0,10):
        name=bamprefix+"_"+str(i+1)+".bam"
        prefix=bamprefix+"_"+str(i+1)
        output.append(name)
        prefixall.append(prefix)
    peaknum=[]
    for i in range(0,10):
        subprocess.call(["samtools","view","-s",str(float((i+1)*0.1)),"-o",output[i],bamfiles])
        #with pysam.AlignmentFile("tmp.bam", "wb", template = bamfile) as outf:
            #print(int((i+1)/10*N),end="\n")
            #rand=random.sample(allreads,int((i+1)/10*N))
            #for readr in rand:
                #outf.write(readr)
            #outf.close()
        #pysam.sort("-@","20","-o",output[i],"tmp.bam")
        #subprocess.call(["rm","tmp.bam"])
        subprocess.call(["samtools","index",output[i]])
        #subprocess.call(["python","/gss1/home/tst/Popera/Popera.py","-d",output[i],"-n",prefixall[i],"-t","5","-l","50","--threads=20"])
        subprocess.call(["macs2","callpeak","-t",output[i],"-n",prefixall[i],"-f","BAMPE","-g","2.2e+09","--nomodel","-q","0.001","--fe-cutoff","5"])
        #inputname=prefixall[i]+"_hotspots.bed"
        inputname=prefixall[i]+"_peaks.narrowPeak"
        hotspot=open(inputname,"r")
        n=len(hotspot.readlines())
        hotspot.close()
        peakfile=prefixall[i]+"_peaks.xls"
        summit=prefixall[i]+"_summits.bed"
        bamindex=output[i]+".bai"
        subprocess.call(["rm",inputname])
        subprocess.call(["rm",peakfile])
        subprocess.call(["rm",summit])
        subprocess.call(["rm",output[i]])
        subprocess.call(["rm",bamindex])
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

for i in range(0,100,10):
    xlabs.append(i/10)
    xlabp.append(str(i+10))
plt.ylabel("Peak number")
plt.xticks(xlabs,xlabp)
plt.legend(label)
plt.savefig(args.outf)
