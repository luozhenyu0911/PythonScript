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
parser.add_argument('--fastq', '-f',type = str, help="Paired end fastq files, Fq1 and Fq2 seperated by comma.",required= True)
parser.add_argument('--index',type = str, help="bowtie2 index prefix",required= True)
parser.add_argument('--outf','-o',type = str, help="Figure name",required= True)
parser.add_argument('--outPeakNum','-p',type = str, help="Peak number file")
parser.add_argument('--tempDir','-d',type = str, help="temp file directory.",required= True)
args = parser.parse_args()

def plotsmooth(values):
    n=len(values)
    num=range(0,n)
    num=np.array(num)
    xnew=np.linspace(num.min(),num.max(),5*n)
    func=interp1d(num,values,kind="cubic")
    ynew=func(xnew)
    plt.plot(xnew,ynew)

def Replace_all(text):
    dic={'.fq':'','.fastq':'','.gz':'',
         '.R1':'','_R1':'','-R1':'',
         '.R2':'','_R2':'','-R2':''}
    for i,j in dic.iteritems():
        text=text.replace(i,j)
    return text

def randomPeak(fastq,index):
    fastq1=fastq.split(",")[0]
    fastq2=fastq.split(",")[1]
    prefixR1=[]
    prefixR2=[]
    prefix=Replace_all(fastq1)
    for i in range(0,20):
        sampleFq1=prefix+"_random_sample"+str(i+1)+'.R1'
        sampleFq2=prefix+"_random_sample"+str(i+1)+'.R2'
        prefixR1.append(sampleFq1)
        prefixR2.append(sampleFq2)
    peaknum=[]
    for i in range(0,20):
        if i < 19:
            subfq1=prefixR1[i]+'.fq'
            subfq2=prefixR2[i]+'.fq'
            outsam=prefix+"_random_sample"+str(i+1)+".sam"
            outbam=prefix+"_random_sample"+str(i+1)+".bam"
            outuniqunsort=prefix+"_random_sample"+str(i+1)+".unsort.bam"
            outuniqbam=prefix+"_random_sample"+str(i+1)+".uniq.bam"
            fraction=(i+1)/20
            #按比例抽取fastq
            with open(subfq1,'w') as f1,open(subfq2,'w') as f2:
                subprocess.call(["seqtk","sample","-s100",fastq1,str(fraction)],stdout=f1)
                subprocess.call(["seqtk","sample","-s100",fastq2,str(fraction)],stdout=f2)
            #比对, subprocess call wait antomatically.
            subprocess.call(["bowtie2","-p","10","-x",index,"-1",subfq1,"-2",subfq2,"-S",outsam])
            subprocess.call(["samtools","sort","-@","10","-O","BAM","-o",outbam,outsam]) #已经排过序，后面的也一定是有序的
            subprocess.call(["samtools","view","-@","10","-h","-q","20",outbam,'-o',outuniqbam])
#             subprocess.call(["samtools","view","-@","10","-h","-q","20",outbam,'-o',outuniqunsort])
#             subprocess.call(["samtools","sort","-@","10","-O","BAM","-o",outuniqbam,outuniqunsort])
            subprocess.call(["samtools","index",outuniqbam])
            #call 核小体
            subprocess.call(["python","/gss1/home/lxx20180918/software/danpos2/danpos.py","dpos",outuniqbam,"-o",args.tempDir,"--paired","1","-p","1e-3","-q","0"])
            inputname=args.tempDir+"/pooled/"+prefix+"_random_sample"+str(i+1)+".uniq.smooth.positions.xls"
            hotspot=open(inputname,"r")
            n=len(hotspot.readlines())
            hotspot.close()
            peaknum.append(n)
            bamindex=outuniqbam+".bai"
            #删除文件
            subprocess.call(["rm",subfq1])
            subprocess.call(["rm",subfq2])
            subprocess.call(["rm",outsam])
            subprocess.call(["rm",outbam])
            #subprocess.call(["rm",outuniqunsort])
            subprocess.call(["rm",outuniqbam])
            subprocess.call(["rm",bamindex])
        else:
            outsam=prefix+"_random_sample"+str(i+1)+".sam"
            outbam=prefix+"_random_sample"+str(i+1)+".bam"
            outuniqunsort=prefix+"_random_sample"+str(i+1)+".unsort.bam"
            outuniqbam=prefix+"_random_sample"+str(i+1)+".uniq.bam"
            #比对
            subprocess.call(["bowtie2","-p","10","-x",index,"-1",fastq1,"-2",fastq2,"-S",outsam])
            subprocess.call(["samtools","sort","-@","10","-O","BAM","-o",outbam,outsam])
            subprocess.call(["samtools","view","-@","10","-h","-q","20",outbam,'-o',outuniqbam])
#             subprocess.call(["samtools","view","-@","10","-h","-q","20",outbam,'-o',outuniqunsort])
#             subprocess.call(["samtools","sort","-@","10","-O","BAM","-o",outuniqbam,outuniqunsort])
            subprocess.call(["samtools","index",outuniqbam])
            #calling
            subprocess.call(["python","/gss1/home/lxx20180918/software/danpos2/danpos.py","dpos",outuniqbam,"-o",args.tempDir,"--paired","1","-p","1e-3","-q","0"])
            inputname=args.tempDir+"/pooled/"+prefix+"_random_sample"+str(i+1)+".uniq.smooth.positions.xls"
            hotspot=open(inputname,"r")
            n=len(hotspot.readlines())
            peaknum.append(n)
            hotspot.close()
            bamindex=outuniqbam+".bai"
            #删除
            subprocess.call(["rm",outsam])
            subprocess.call(["rm",outbam])
            #subprocess.call(["rm",outuniqunsort])
            subprocess.call(["rm",outuniqbam])
            subprocess.call(["rm",bamindex])
    return peaknum

plt.switch_backend('PDF')

outP=open(args.outPeakNum,"w")
label=[]
fastq=args.fastq
index=args.index

fq1=fastq.split(",")[0]
fq2=fastq.split(",")[1]

prefix=Replace_all(fq1)
label.append(prefix)
valuelist=randomPeak(fastq,index)
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
