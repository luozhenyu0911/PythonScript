from __future__ import print_function
from __future__ import division
import sys
import argparse
import pysam
import numpy as np
import matplotlib.pyplot as plt
from pyfasta import Fasta
from scipy.interpolate import interp1d

parser = argparse.ArgumentParser(description="This tool can profile average footprint metaplot on binding sites", formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('--referencePoint',type=str,help="<referencePoint {start or end}>"+"\t"+"plot the distribution of the reads at reference point from the genome.")
parser.add_argument('--scaleRegion',type=int,help="<Bins for region>"+"\t"+"plot the distribution of the reads over the region from the genome.")
parser.add_argument('--bam', '-b',type= str,help="Bam file for open chromatin assay",required= True)
parser.add_argument('--bedlist', '-r',type= str,help="Bed file for motif binding sites,one file per line",required= True)
parser.add_argument('--bias', '-f',type= str,help="Biased file for nake DNA library or Open chromatin library",required= True)
parser.add_argument('--genome', '-g',type= str,help="fasta file for the genome",required= True)
parser.add_argument('--upstream', '-u',type= int,help="Length of the upstream reigon",required= True)
parser.add_argument('--downstream', '-d',type= int,help="Length of the downstream reigon",required= True)
parser.add_argument('--windowsize', '-w',type= int,help="Length of the windowsize for flank region",required= True)
parser.add_argument('--correct', '-c',type= str,choices=["yes","no"],help="whether to correct seqence bias ",default="yes")
parser.add_argument('--outfileprefix','-o',type=str,help="<OUTPUT file name Prefix>"+"\t"+"It will output full matrix,average data and average profile.",required= True)
parser.add_argument('--outfigure','-O',type=str,help="<OUTPUT figure name>"+"\t"+"It will output all average plot.",required= True)
args = parser.parse_args()

bamIn=args.bam
bedlist=args.bedlist
biasfile=args.bias
genome=args.genome
ustream=args.upstream
dstream=args.upstream
wdsize=args.windowsize
prefix=args.outfileprefix
referencePoint=args.referencePoint
scaleRegion=args.scaleRegion
correct=args.correct
figure=args.outfigure

plt.switch_backend('PDF')


def bamInfo(Input):
    ReadInfo={"+":{},"-":{}}
    for read in Input.fetch():
        chrom=read.reference_name
        if chrom in ReadInfo["+"]:
            pass
        else:
            ReadInfo["+"][chrom]={}
            ReadInfo["-"][chrom]={}
        cutpoint=read.reference_end+1 if read.is_reverse else read.reference_start+1
        strand="-" if read.is_reverse else "+"
        if cutpoint in ReadInfo[strand][chrom]:
            ReadInfo[strand][chrom][int(cutpoint)]+=1
        else:
            ReadInfo[strand][chrom][int(cutpoint)]=1
    return ReadInfo

def partition(lst, n):
    division = len(lst) / float(n)
    slit=[lst[int(round(division * i)): int(round(division * (i + 1)))] for i in xrange(n) ]
    return [sum(slit[i])/len(slit[i]) for i in range(0,len(slit))]    

def Readscount(bed,cutdict):
    [chrom,start,end,strand]=bed
    count=[cutdict[strand][chrom][i] if i in cutdict[strand][chrom] else 0 for i in range(start,end+1)]
    return count

def ReadscountAll(bed,cutdict):
    [chrom,start,end]=bed
    countf=[cutdict["+"][chrom][i] if i in cutdict["+"][chrom] else 0 for i in range(start,end+1)]
    countr=[cutdict["-"][chrom][i] if i in cutdict["-"][chrom] else 0 for i in range(start,end+1)]
    return list(np.sum([countf,countr],axis=0))

def countCorrect(bed,cutdict,fasta,biasfile):
    [chrom,start,end]=bed
    countf=[BaseCorrection(chrom,i,fasta,cutdict,biasfile,"+") if i in cutdict["+"][chrom] else 0 for i in range(start,end+1)]
    countr=[BaseCorrection(chrom,i,fasta,cutdict,biasfile,"-") if i in cutdict["-"][chrom] else 0 for i in range(start,end+1)]      
    return list(np.sum([countf,countr],axis=0))

def Base6mer(chrom,loci,fasta,strand):
    cut6mer=str(fasta.sequence({'chr':chrom,'start':loci-3,'stop':loci+2,'strand':strand}))
    if len(cut6mer) ==0:
        cut6mer = "NNNNNN"
    return cut6mer

def BaseCorrection(chrom,loci,fasta,cutdict,biasfile,strand):
    if cutdict[strand][chrom][loci]:
        cut50bp=sum(Readscount([chrom,loci-25,loci+26,strand],cutdict))
        cut506mer=sum([biasfile[Base6mer(chrom,i,fasta,strand)] if "N" not in Base6mer(chrom,i,fasta,strand) else 0 for i in range(loci-25,loci+26)])
        if "N" not in Base6mer(chrom,loci,fasta,strand):
            correctValue=(cutdict[strand][chrom][loci]+1)/(biasfile[Base6mer(chrom,loci,fasta,strand)]/cut506mer*cut50bp+1)
        else:
            correctValue=cutdict[strand][chrom][loci]
    else:
        if "N" not in Base6mer(chrom,loci,fasta,strand):
            correctValue=1/(biasfile[Base6mer(chrom,loci,fasta,strand)]/cut506mer*cut50bp+1)
        else:
            correctValue=0
    return correctValue

def ScaleCalculate(bedfile,cutdict,up,down,wd,binNum,biasfile,fasta):
    allsite=[]
    for line in bedfile:
        row=line.rstrip().split("\t")
        upregion=[row[0],int(row[1])-up,int(row[1])-1]
        downregion=[row[0],int(row[2])+1,int(row[2])+down]
        upcounts=partition(countCorrect(upregion,cutdict,fasta,biasfile),int(up/wd))
        downcounts=partition(countCorrect(downregion,cutdict,fasta,biasfile),int(down/wd))
        regioncount=partition(countCorrect([row[0],int(row[1]),int(row[2])],cutdict,fasta,biasfile),binNum)
        allcount=upcounts+regioncount+downcounts
        if "+" in row[5]:
            allsite.append(allcount)
        else:
            allsite.append(allcount[::-1])
    return allsite

def ScaleCalculateRaw(bedfile,cutdict,up,down,wd,binNum):
    allsite=[]
    for line in bedfile:
        row=line.rstrip().split("\t")
        upregion=[row[0],int(row[1])-up,int(row[1])-1]
        downregion=[row[0],int(row[2])+1,int(row[2])+down]
        region=[row[0],int(row[1]),int(row[2])]
        upcounts=partition(ReadscountAll(upregion,cutdict),int(up/wd))
        downcounts=partition(ReadscountAll(downregion,cutdict),int(down/wd))
        regioncounts=partition(ReadscountAll(region,cutdict),binNum)
        allcount=upcounts+regioncounts+downcounts
        if "+" in row[5]:
            allsite.append(allcount)
        else:
            allsite.append(allcount[::-1])
    return allsite
         
def plotsmooth(values):
    n=len(values)
    num=range(0,n)
    num=np.array(num)
    xnew=np.linspace(num.min(),num.max(),5*n)
    func=interp1d(num,values,kind="cubic")
    ynew=func(xnew)
    plt.plot(xnew,ynew)

def ScalePlot(mat,up,down,wd,binNum,label):
    upbin=int(up/wd)
    downbin=int(up/wd)
    if binNum > 1:
        labs=["-"+str(up)+"bp","start","end",str(down)+"bp"]
        labp=[0,upbin,upbin+binNum-1,upbin+downbin+binNum-1]
    else:
        labs=["-"+str(up)+"bp","motif",str(down)+"bp"]
        labp=[0,upbin,upbin+downbin]
    value=np.array(mat).mean(axis=0)
    plotsmooth(value)
    plt.title(label)
    plt.ylabel("bias correct cuts")
    plt.xticks(labp,labs)

def DataSave(mat,prefix):
    average=np.array(mat).mean(axis=0)
    outf=open(prefix+"_dens.txt","w")
    print(*average,end="\n",sep="\t",file=outf)
    outf.close()
    outp=open(prefix+"_mat.txt","w")
    header=["h"]+["bin_"+str(i) for i in range(0,len(average))]
    print(*header,end="\n",sep="\t",file=outp)
    for line in mat:
        print("row",*line,end="\n",sep="\t",file=outp)
    outp.close()

def loadBias(biasfile):
    a=open(biasfile,"r")
    bias={}
    for line in a:
        row=line.rstrip().split("\t")
        bias[row[0]]=float(row[1])
    return bias

def main():
    if scaleRegion:
        opbam=pysam.AlignmentFile(bamIn,"rb")
        cuts=bamInfo(opbam)
        biaslist=loadBias(biasfile)
        fa=Fasta(genome)
        beds=[line.rstrip() for line in open(bedlist,"r")]
        plist=[line.rstrip() for line in open(prefix,"r")]
        plotnum=len(beds)
        if plotnum%(int(plotnum ** 0.5)) > 0:
            c,r=plotnum//int(plotnum ** 0.5)+1,int(plotnum ** 0.5)
        else:
            c,r=plotnum//int(plotnum ** 0.5),int(plotnum ** 0.5)
        plt.figure(figsize=(5*c,5*r))
        for i in range(0,len(beds)):
            inputbed=open(beds[i],"r")
            if "yes" in correct:
                cutmat=ScaleCalculate(inputbed,cuts,ustream,dstream,wdsize,scaleRegion,biaslist,fa)
            else:
                cutmat=ScaleCalculateRaw(inputbed,cuts,ustream,dstream,wdsize,scaleRegion)       
            DataSave(cutmat,plist[i])
            plt.subplot(r,c,i+1)
            ScalePlot(cutmat,ustream,dstream,wdsize,scaleRegion,plist[i])
        plt.savefig(figure)

if __name__ == '__main__':
    main()
