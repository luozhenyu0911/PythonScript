# -*- coding: UTF-8 -*-
from __future__ import print_function
from __future__ import division
#from numba import jit
import argparse
#import sys, getopt
import pysam 
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
plt.switch_backend('PDF')
#@jit

example_text = '''example:

 python **.py --referencePoint start/end -b **.bam,**.bam... -G **.bed,**.bed... -g **.txt,**.txt... -u 1000 -d 1000 -w 50 -l **,** -o **_dens,**_matrix -F **.pdf
 python **.py --scaleRegion -b **.bam,**bam... -G **.bed,**.bed... -g **.txt,**.txt... -u 1000 -d 1000 -w 50 -N 100 -l **,** -o **_dens,**_matrix -F **.pdf
'''
parser = argparse.ArgumentParser(description="This tool creates scores and profile plot over sets of genomic regions. Typically, these regions are genes, but any other regions defined in BED will work.", 
                                 formatter_class= argparse.RawTextHelpFormatter,
                                 prog='base_maker',
                                 epilog=example_text)
group=parser.add_mutually_exclusive_group(required= True)
group.add_argument('--referencePoint',type=str,help="<referencePoint {start or end}>"+"\t"+"plot the distribution of the reads at reference point from the genome.")
group.add_argument('--scaleRegion',action='store_true',help="<Scaled Region>"+"\t"+"plot the distribution of the reads over the region from the genome.")

parser.add_argument('--control', '-c',type= str,help="<Input bam files>"+"\t"+"Control data files(Support for multiple files, Separated by commas, the order is consistent with bam files), should be sorted bam format.")
parser.add_argument('--bam','-b',type= str,help="<bamfiles>"+"\t"+"data files(Support for multiple files, Separated by commas), should be sorted bam format.")

group2=parser.add_mutually_exclusive_group(required= True)
group2.add_argument('--bed','-G',type=str,help="<bedfiles>"+"\t"+"data files(Support for multiple files, Separated by commas), should be original bed files.")
group2.add_argument('--bedList',type= str,help="<bedfiles>"+"\t"+"bed file text file.")

parser.add_argument('--gene','-g',type=str,help="<gene list files>"+"\t"+"data files(Support for multiple files, Separated by commas), should be a row genes for every file.")
parser.add_argument('--upstream','-u',type=int,help="<upstream distance>"+"\t"+"Distance upstream of the reference-point selected.",required= True)
parser.add_argument('--downstream','-d',type=int,help="<downstream distance>"+"\t"+"Distance downstream of the reference-point selected.",required= True)
parser.add_argument('--windows','-w',type=int,help="<window size>"+"\t"+"Length, in bases, binSize for averaging the score over the regions length.",required= True)
parser.add_argument('--Nbin','-N',type=int,help="<Number of bins>"+"\t"+"The number of equally divided bins, ONLY works when used with <--scaleRegion>")

group3=parser.add_mutually_exclusive_group(required= True)
group3.add_argument('--legends','-l',type=str,help="<legend names>"+"\t"+"Name of the plot legend. (Support for multiple legends).")
group3.add_argument('--legendList',type= str,help="<legend names>"+"\t"+"legend name text file.")

parser.add_argument('--outfile','-o',type=str,help="<outputfiles average data,final matrix>"+"\t"+"File name to save the average data file and final matrix file. (Separated by commas).",required= True)
parser.add_argument('--figure','-F',type=str,help="<Figure name>"+"\t"+"Name of the plot name. (plot file format is PDF).",required= True)
args = parser.parse_args()

if args.control:
    Inputbamlist=args.control.split(",")
else:
    Inputbamlist=None
if args.gene:
    genelist=args.gene.split(",")
else:
    genelist=None

bamlist=args.bam.split(",")

if args.bed:
    bedlist=args.bed.split(",")
else:
    bedlist=[]
    with open(args.bedList,'r') as bedtxt:
        for line in bedtxt:
            line=line.strip().split()
            bedlist.append(line[0])

up=args.upstream
L=str(up/1000)
down=args.downstream
R=str(down/1000)
referencePoint=args.referencePoint
scaleRegion=args.scaleRegion
wdsize=args.windows
blocknum=args.Nbin

if args.legends:
    legendname=args.legends.split(",")
else:
    legendname=[]
    with open(args.legendList,'r') as legendtxt:
        for line in legendtxt:
            line=line.strip().split()
            legendname.append(line[0])
figurename=args.figure
outf=args.outfile.split(",")


def bamTodict(Input):
    ReadInfo={}
    for read in Input.fetch():
        ID=read.query_name
        if ID in ReadInfo:
            if int(read.reference_start) < int(ReadInfo[ID][1]):
                ReadInfo[ID][1]=int(read.reference_start)
            else:
                ReadInfo[ID][2]=int(read.reference_end)
        else:
            a=[]
            a.append(str(read.reference_name)) #染色体
            a.append(int(read.reference_start)) #左切点
            a.append(int(read.reference_end)) # 右切点
            #a.append(int(a[2])-int(a[1])+1)
            ReadInfo[ID]=a
    region={}
    for value in ReadInfo.values():
        chrom=value[0]
        start=int(value[1])
        end=int(value[2])
        if chrom in region:
            pass
        else:
            region[chrom]={}
        if str(start) in region[chrom]:
            region[chrom][str(start)]+=1
        else:
            region[chrom][str(start)]=1
        if str(end) in region[chrom]:
            region[chrom][str(end)]+=1
        else:
            region[chrom][str(end)]=1
    return region
    
def GenomebedCount(bedf,readdict,count,wdsize,up,down):
    bed6col=open(bedf,"r")
    sumReads=0
    genomeRead={}
    for bed in bed6col:
        row2=bed.rstrip().split()
        chrom=row2[0]
        sublist=[]
        if "+" in row2[5]:
            left=int(row2[1])-up
            right=int(row2[1])+down
            for k in range(int(left),int(right),wdsize):
                for j in range(k,(k+wdsize)):
                    if str(j) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(j)])
                    else:
                        pass
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            genomeRead[row2[3]]=sublist
        else:
            left=int(row2[2])-down
            right=int(row2[2])+up
            for k in range(int(right),int(left),-wdsize):
                for j in range(k,(k-wdsize),-1):
                    if str(j) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(j)])
                    else:
                        pass
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            genomeRead[row2[3]]=sublist
    return genomeRead
    bed6col.close()

def GenomebedCountEnd(bedf,readdict,count,wdsize,up,down):
    bed6col=open(bedf,"r")
    sumReads=0
    genomeRead={}
    for bed in bed6col:
        row2=bed.rstrip().split()
        chrom=row2[0]
        sublist=[]
        if "+" in row2[5]:
            left=int(row2[2])-up
            right=int(row2[2])+down
            for k in range(int(left),int(right),wdsize):
                for j in range(k,(k+wdsize)):
                    if str(j) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(j)])
                    else:
                        pass
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            genomeRead[row2[3]]=sublist
        else:
            left=int(row2[1])-down
            right=int(row2[1])+up
            for k in range(int(right),int(left),-wdsize):
                for j in range(k,(k-wdsize),-1):
                    if str(j) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(j)])
                    else:
                        pass
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            genomeRead[row2[3]]=sublist
    return genomeRead
    bed6col.close()

def bedCount(bedf,readdict,count,wdsize,up,down):
    bed6col=open(bedf,"r")
    Lis=[]
    sumReads=0  
    for bed in bed6col:
        row2=bed.rstrip().split()
        chrom=row2[0]
        sublist=[]
        if "+" in row2[5]:
            left=int(row2[1])-up
            right=int(row2[1])+down
            for k in range(int(left),int(right),wdsize):
                for j in range(k,(k+wdsize)):
                    if str(j) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(j)])
                    else:
                        pass
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            Lis.append(sublist)
        else:
            left=int(row2[2])-down
            right=int(row2[2])+up
            for k in range(int(right),int(left),-wdsize):
                for j in range(k,(k-wdsize),-1):
                    if str(j) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(j)])
                    else:
                        pass
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            Lis.append(sublist)
    return Lis
    bed6col.close()

def bedCountEnd(bedf,readdict,count,wdsize,up,down):
    bed6col=open(bedf,"r")
    Lis=[]
    sumReads=0
    for bed in bed6col:
        row2=bed.rstrip().split()
        chrom=row2[0]
        sublist=[]
        if "+" in row2[5]:
            left=int(row2[2])-up
            right=int(row2[2])+down
            for k in range(int(left),int(right),wdsize):
                for j in range(k,(k+wdsize)):
                    if str(j) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(j)])
                    else:
                        pass
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            Lis.append(sublist)
        else:
            left=int(row2[1])-down
            right=int(row2[1])+up
            for k in range(int(right),int(left),-wdsize):
                for j in range(k,(k-wdsize),-1):
                    if str(j) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(j)])
                    else:
                        pass
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            Lis.append(sublist)
    return Lis
    bed6col.close()

def BodyGenomeCount(bedf,readdict,count,wdsize,up,down,BlockNum):
    bed6col=open(bedf,"r")
    sumReads=0
    genomeRead={}
    for bed in bed6col:
        row2=bed.rstrip().split()
        chrom=row2[0]
        left=int(row2[1])
        right=int(row2[2])
        sublist=[]
        remainder = (right-left)%BlockNum
        #如果区间能被整除
        if remainder == 0:
            blockSize=int((right-left)/BlockNum)
            reside=0
        #如果区间不能被整除
        else:
            reside=BlockNum-int(((right-left)%BlockNum))
            blockSize=int(((right-left)+reside)/BlockNum)
        if "+" in row2[5]:
            for k in range(left-up,left,wdsize):
                for j in range(k,(k+wdsize)):
                    if str(j) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(j)])
                    else:
                        pass
              #  sumReads=int((sumReads/10)*10)
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            for x in range(left,right+reside,blockSize):
                for y in range(x,(x+blockSize)):
                    if str(y) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(y)])
                    else:
                        pass
              #  sumReads=int((sumReads/blockSize)*10)
                sublist.append((int(sumReads)*1000000)/(count*blockSize))
                sumReads=0
            for u in range(right,right+down,wdsize):
                for v in range(u,(u+wdsize)):
                    if str(v) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(v)])
                    else:
                        pass
             #   sumReads=int((sumReads/10)*10)
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            genomeRead[row2[3]]=sublist
        else:
            for k in range(right+up,right,-wdsize):
                for j in range(k,(k-wdsize),-1):
                    if str(j) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(j)])
                    else:
                        pass
            #    sumReads=int((sumReads/10)*10)
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            for x in range(right,left-reside,-blockSize):
                for y in range(x,(x-blockSize),-1):
                    if str(y) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(y)])
                    else:
                        pass
           #     sumReads=int((sumReads/blockSize)*10)
                sublist.append((int(sumReads)*10000000)/(count*blockSize))
                sumReads=0
            for u in range(left,left-down,-wdsize):
                for v in range(u,(u-wdsize),-1):
                    if str(v) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(v)])
                    else:
                        pass
          #      sumReads=int((sumReads/10)*10)
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            genomeRead[row2[3]]=sublist
    return genomeRead
    bed6col.close()

def BodyCount(bedf,readdict,count,wdsize,up,down,BlockNum):
    bed6col=open(bedf,"r")
    Lis=[]
    sumReads=0
    for bed in bed6col:
        row2=bed.rstrip().split()
        chrom=row2[0]
        left=int(row2[1])
        right=int(row2[2])
        sublist=[]
        remainder = (right-left)%BlockNum
        #如果区间能被整除
        if remainder == 0:
            blockSize=int((right-left)/BlockNum)
            reside=0
        #如果区间不能被整除
        else:
            reside=BlockNum-int(((right-left)%BlockNum))
            blockSize=int(((right-left)+reside)/BlockNum)
        if "+" in row2[5]:
            for k in range(left-up,left,wdsize):
                for j in range(k,(k+wdsize)):
                    if str(j) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(j)])
                    else:
                        pass
               # sumReads=int((sumReads/10)*10)
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            for x in range(left,right+reside,blockSize):
                for y in range(x,(x+blockSize)):
                    if str(y) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(y)])
                    else:
                        pass
                #sumReads=int((sumReads/blockSize)*10)
                sublist.append((int(sumReads)*10000000)/(count*blockSize))
                sumReads=0
            for u in range(right,right+down,wdsize):
                for v in range(u,(u+wdsize)):
                    if str(v) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(v)])
                    else:
                        pass
               # sumReads=int((sumReads/10)*10)
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            Lis.append(sublist)
        else:
            for k in range(right+up,right,-wdsize):
                for j in range(k,(k-wdsize),-1):
                    if str(j) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(j)])
                    else:
                        pass
              #  sumReads=int((sumReads/10)*10)
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            for x in range(right,left-reside,-blockSize):
                for y in range(x,(x-blockSize),-1):
                    if str(y) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(y)])
                    else:
                        pass
              #  sumReads=int((sumReads/blockSize)*10)
                sublist.append((int(sumReads)*10000000)/(count*blockSize))
                sumReads=0
            for u in range(left,left-down,-wdsize):
                for v in range(u,(u-wdsize),-1):
                    if str(v) in readdict[chrom]:
                        sumReads+=int(readdict[chrom][str(v)])
                    else:
                        pass
              #  sumReads=int((sumReads/10)*10)
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            Lis.append(sublist)
    return Lis
    bed6col.close()

def bamtocountStart(bed,Input,windowSize,Up,Down):
    bamf=pysam.AlignmentFile(Input,"rb")
    region=bamTodict(bamf)
    num=bamf.count()
    allgene=GenomebedCount(bed,region,num,windowSize,Up,Down)
    return allgene

def bamTocountStart(bed,Input,windowSize,Up,Down):
    bamf=pysam.AlignmentFile(Input,"rb")
    region=bamTodict(bamf)
    num=bamf.count()
    allvalue=bedCount(bed,region,num,windowSize,Up,Down)
    return allvalue

def bamtocountEnd(bed,Input,windowSize,Up,Down):
    bamf=pysam.AlignmentFile(Input,"rb")
    region=bamTodict(bamf)
    num=bamf.count()
    allgene=GenomebedCountEnd(bed,region,num,windowSize,Up,Down)
    return allgene

def bamTocountEnd(bed,Input,windowSize,Up,Down):
    bamf=pysam.AlignmentFile(Input,"rb")
    region=bamTodict(bamf)
    num=bamf.count()
    allvalue=bedCountEnd(bed,region,num,windowSize,Up,Down)
    return allvalue

def bamtobodycount(bed,Input,windowSize,Up,Down,BlockNUM):
    bamf=pysam.AlignmentFile(Input,"rb")
    region=bamTodict(bamf)
    num=bamf.count()
    bodyvalue=BodyCount(bed,region,num,windowSize,Up,Down,BlockNUM)
    return bodyvalue

def bamGenometobodycount(bed,Input,windowSize,Up,Down,BlockNUM):
    bamf=pysam.AlignmentFile(Input,"rb")
    region=bamTodict(bamf)
    num=bamf.count()
    Genomebodyvalue=BodyGenomeCount(bed,region,num,windowSize,Up,Down,BlockNUM)
    return Genomebodyvalue

def Listgene(geneFile,GenomeDict):
    Input=open(geneFile,"r")
    Lismat=[]
    tmp=[]
    for line in Input:
        Id=line.rstrip()
        tmp.append(Id)
    for g in range(len(tmp)):    
        if tmp[g] in GenomeDict:
            Lismat.append(GenomeDict[tmp[g]])
    return Lismat

def plotsmooth(Values):
    n=len(Values)
    num=range(0,n)
    num=np.array(num)
    xnew=np.linspace(num.min(),num.max(),5*n)
    func=interp1d(num,Values,kind="cubic")
    ynew=func(xnew)
    plt.plot(xnew,ynew)

if Inputbamlist:

    if referencePoint:
        Ylength=int(up+down)
        ylength=int(Ylength/wdsize)
        left=int(0)
        medium=int(ylength/2)
        right=int(ylength)
        NumOfBin = int((int(up+down))/wdsize)
        if referencePoint=="start":
            if genelist is None:
                out_profile=open(outf[0],"w")
                out_matrix=open(outf[1],"w")
                out_figure=open(figurename,"w")
                figurelist=[]
                header = ["label","num"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
                print(*header,sep="\t",end="\n",file=out_matrix)
                for Bam in range(len(bamlist)):
                    for Bed in range(len(bedlist)):
                        values=bamTocountStart(bedlist[Bed],bamlist[Bam],wdsize,up,down)
                        Invalues=bamTocountStart(bedlist[Bed],Inputbamlist[Bam],wdsize,up,down)
                        tup1=np.array(values)	
                        Intup=np.array(Invalues)
                        tup=tup1-Intup
                        tup[tup<0]=0
                        tup2=tup.mean(axis=0)
                        figurelist.append(tup2)
                        Values=list(tup)
                        for index in range(len(Values)):
                            print(legendname[Bed],str(index+1),*Values[index],sep="\t",end="\n",file=out_matrix)
                print("region","model","value",sep="\t",end="\n",file=out_profile)
                N = 0
                for lis in figurelist:
                    n = 0
                    for l in lis:
                        print(n*wdsize,legendname[N],l,sep="\t",end="\n",file=out_profile)
                        n += 1
                    N += 1
                for line in figurelist:
                    plotsmooth(line)
                plt.legend(legendname)
                plt.ylabel("Normalized reads count")
                plt.xticks([left,medium,right],["-"+L+"kb","TSS",R+"kb"])
                plt.savefig(out_figure)	
                out_profile.close()
                out_matrix.close()
                out_figure.close()
            else:
                out_figure=open(figurename,"w")
                out_profile=open(outf[0],"w")
                out_matrix=open(outf[1],"w")
                figurelist=[]
                header = ["label","num"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
                print(*header,sep="\t",end="\n",file=out_matrix)
                for Bam in range(len(bamlist)):
                    for Bed in range(len(bedlist)):
                        for Gene in range(len(genelist)):
                            values=bamtocountStart(bedlist[Bed],bamlist[Bam],wdsize,up,down)
                            Invalues=bamtocountStart(bedlist[Bed],Inputbamlist[Bam],wdsize,up,down)
                            finalmat=Listgene(genelist[Gene],values)
                            Infinalmat=Listgene(genelist[Gene],Invalues)
                            tup1=np.array(finalmat)
                            Intup=np.array(Infinalmat)
                            tup=tup1-Intup
                            tup[tup<0]=0
                            tup2=tup.mean(axis=0)
                            figurelist.append(tup2)
                            Values=list(tup)
                            for index in range(len(Values)):
                                print(legendname[Bed],str(index+1),*Values[index],sep="\t",end="\n",file=out_matrix)
                print("region","model","value",sep="\t",end="\n",file=out_profile)
                N = 0
                for lis in figurelist:
                    n = 0
                    for l in lis:
                        print(n*wdsize,legendname[N],l,sep="\t",end="\n",file=out_profile)
                        n += 1
                    N += 1
                for line in figurelist:
                    plotsmooth(line)
                plt.legend(legendname)
                plt.ylabel("Normalized reads count")
                plt.xticks([left,medium,right],["-"+L+"kb","TSS",R+"kb"])
                plt.savefig(out_figure)     
                out_figure.close()
                out_profile.close()
                out_matrix.close()
        elif referencePoint=="end":
            if genelist is None:
                out_profile=open(outf[0],"w")
                out_matrix=open(outf[1],"w")
                out_figure=open(figurename,"w")
                figurelist=[]
                header = ["label","num"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
                print(*header,sep="\t",end="\n",file=out_matrix)
                for Bam in range(len(bamlist)):
                    for Bed in range(len(bedlist)):
                        values=bamTocountEnd(bedlist[Bed],bamlist[Bam],wdsize,up,down)
                        Invalues=bamTocountEnd(bedlist[Bed],Inputbamlist[Bam],wdsize,up,down)
                        tup1=np.array(values)
                        Intup=np.array(Invalues)
                        tup=tup1-Intup
                        tup[tup<0]=0
                        tup2=tup.mean(axis=0)
                        Values=list(tup)
                        figurelist.append(tup2)
                        for index in range(len(Values)):
                            print(legendname[Bed],str(index+1),*Values[index],sep="\t",end="\n",file=out_matrix)
                print("region","model","value",sep="\t",end="\n",file=out_profile)
                N = 0
                for lis in figurelist:
                    n = 0
                    for l in lis:
                        print(n*wdsize,legendname[N],l,sep="\t",end="\n",file=out_profile)
                        n += 1
                    N += 1
                for line in figurelist:
                    plotsmooth(line)
                plt.legend(legendname)
                plt.ylabel("Normalized reads count")
                plt.xticks([left,medium,right],["-"+L+"kb","TSS",R+"kb"])
                plt.savefig(out_figure)     
                out_figure.close()
                out_profile.close()
                out_matrix.close()
            else:
                out_profile=open(outf[0],"w")
                out_matrix=open(outf[1],"w")
                out_figure=open(figurename,"w")
                figurelist=[]
                header = ["label","num"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
                print(*header,sep="\t",end="\n",file=out_matrix)
                for Bam in range(len(bamlist)):
                    for Bed in range(len(bedlist)):
                        for Gene in range(len(genelist)):
                            values=bamtocountEnd(bedlist[Bed],bamlist[Bam],wdsize,up,down)
                            Invalues=bamtocountEnd(bedlist[Bed],Inputbamlist[Bam],wdsize,up,down)
                            finalmat=Listgene(genelist[Gene],values)
                            Infinalmat=Listgene(genelist[Gene],Invalues)
                            tup1=np.array(finalmat)
                            Intup=np.array(Infinalmat)
                            tup=tup1-Intup
                            tup[tup<0]=0
                            tup2=tup.mean(axis=0)
                            Values=list(tup)
                            figurelist.append(tup2)
                            for index in range(len(Values)):
                                print(legendname[Bed],str(index+1),*Values[index],sep="\t",end="\n",file=out_matrix)
                print("region","model","value",sep="\t",end="\n",file=out_profile)
                N = 0
                for lis in figurelist:
                    n = 0
                    for l in lis:
                        print(n*wdsize,legendname[N],l,sep="\t",end="\n",file=out_profile)
                        n += 1
                    N += 1
                for line in figurelist:
                    plotsmooth(line)
                plt.legend(legendname)
                plt.ylabel("Normalized reads count")
                plt.xticks([left,medium,right],["-"+L+"kb","TSS",R+"kb"])
                plt.savefig(out_figure)     
                out_figure.close()
                out_profile.close()
                out_matrix.close()
        else:
            pass
	#print("Usage: python tmp1.py <<--scale-regions or --reference-point>> -b <bamfile> -G <bedfile> -g <genelist> -r <referencePoint {start or end}> -u <upstream_distance> -d <downstream_distance> -w <window_size> -bn <Equal_block_number> -o <outputfile>")

    elif scaleRegion:
        left=int(0)
        mediumA=int(up/wdsize)
        mediumB=mediumA+blocknum
        right=int(down/wdsize)+mediumA+blocknum
        NumOfBin = int(up/wdsize)+blocknum+int(down/wdsize)
        if genelist is None:
            out_profile=open(outf[0],"w")
            out_matrix=open(outf[1],"w")
            out_figure=open(figurename,"w")
            figurelist=[]
            header = ["label","num"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
            print(*header,sep="\t",end="\n",file=out_matrix)
            for Bam in range(len(bamlist)):
                for Bed in range(len(bedlist)):
                    values=bamtobodycount(bedlist[Bed],bamlist[Bam],wdsize,up,down,blocknum)
                    Invalues=bamtobodycount(bedlist[Bed],Inputbamlist[Bam],wdsize,up,down,blocknum)
                    Intup=np.array(Invalues)
                    tup1=np.array(values)
                    tup=tup1-Intup
                    tup[tup<0]=0
                    tup2=tup.mean(axis=0)
                    Values=list(tup)
                    figurelist.append(tup2)
                    for index in range(len(Values)):
                        print(legendname[Bed],str(index+1),*Values[index],sep="\t",end="\n",file=out_matrix)
            print("region","model","value",sep="\t",end="\n",file=out_profile)
            N = 0
            for lis in figurelist:
                n = 0
                for l in lis:
                    print(n*wdsize,legendname[N],l,sep="\t",end="\n",file=out_profile)
                    n += 1
                N += 1
            for line in figurelist:
                plotsmooth(line)
            plt.legend(legendname)
            plt.ylabel("Normalized reads count")
            plt.xticks([left,mediumA,mediumB,right],["-"+L+"kb","start","end",R+"kb"])
            plt.savefig(out_figure)     
            out_figure.close()
            out_profile.close()
            out_matrix.close()
        else:
            out_profile=open(outf[0],"w")
            out_matrix=open(outf[1],"w")
            out_figure=open(figurename,"w")
            figurelist=[]
            header = ["label","num"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
            print(*header,sep="\t",end="\n",file=out_matrix)
            for Bam in range(len(bamlist)):
                for Bed in range(len(bedlist)):
                    for Gene in range(len(genelist)):
                        values=bamGenometobodycount(bedlist[Bed],bamlist[Bam],wdsize,up,down,blocknum)
                        Invalues=bamGenometobodycount(bedlist[Bed],Inputbamlist[Bam],wdsize,up,down,blocknum)
                        finalmat=Listgene(genelist[Gene],values)
                        Infinalmat=Listgene(genelist[Gene],Invalues)
                        tup1=np.array(finalmat)
                        Intup=np.array(Infinalmat)
                        tup=tup1-Intup
                        tup[tup<0]=0
                        tup2=tup.mean(axis=0)
                        Values=list(tup)
                        figurelist.append(tup2)
                        for index in range(len(Values)):
                            print(legendname[Bed],str(index+1),*Values[index],sep="\t",end="\n",file=out_matrix)
            print("region","model","value",sep="\t",end="\n",file=out_profile)
            N = 0
            for lis in figurelist:
                n = 0
                for l in lis:
                    print(n*wdsize,legendname[N],l,sep="\t",end="\n",file=out_profile)
                    n += 1
                N += 1
            for line in figurelist:
                plotsmooth(line)
            plt.legend(legendname)
            plt.ylabel("Normalized reads count")
            plt.xticks([left,mediumA,mediumB,right],["-"+L+"kb","start","end",R+"kb"])
            plt.savefig(out_figure)     
            out_figure.close()
            out_profile.close()
            out_matrix.close()
    else:
        pass
    
else:

    if referencePoint:
        Ylength=int(up+down)
        ylength=int(Ylength/wdsize)
        left=int(0)
        medium=int(ylength/2)
        right=int(ylength)
        NumOfBin = int((int(up+down))/wdsize)
        if referencePoint=="start":
            if genelist is None:
                out_profile=open(outf[0],"w")
                out_matrix=open(outf[1],"w")
                out_figure=open(figurename,"w")
                figurelist=[]
                header = ["label","num"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
                print(*header,sep="\t",end="\n",file=out_matrix)
                for Bam in range(len(bamlist)):
                    for Bed in range(len(bedlist)):
                        values=bamTocountStart(bedlist[Bed],bamlist[Bam],wdsize,up,down)
                        tup1=np.array(values)	
                        tup2=tup1.mean(axis=0)
                        figurelist.append(tup2)
                        for index in range(len(values)):
                            print(legendname[Bed],str(index+1),*values[index],sep="\t",end="\n",file=out_matrix)
                print("region","model","value",sep="\t",end="\n",file=out_profile)
                N = 0
                for lis in figurelist:
                    n = 0
                    for l in lis:
                        print(n*wdsize,legendname[N],l,sep="\t",end="\n",file=out_profile)
                        n += 1
                    N += 1
                for line in figurelist:
                    plotsmooth(line)
                plt.legend(legendname)
                plt.ylabel("Normalized reads count")
                plt.xticks([left,medium,right],["-"+L+"kb","TSS",R+"kb"])
                plt.savefig(out_figure)	
                out_profile.close()
                out_matrix.close()
                out_figure.close()
            else:
                out_figure=open(figurename,"w")
                out_profile=open(outf[0],"w")
                out_matrix=open(outf[1],"w")
                figurelist=[]
                header = ["label","num"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
                print(*header,sep="\t",end="\n",file=out_matrix)
                for Bam in range(len(bamlist)):
                    for Bed in range(len(bedlist)):
                        for Gene in range(len(genelist)):
                            values=bamtocountStart(bedlist[Bed],bamlist[Bam],wdsize,up,down)
                            finalmat=Listgene(genelist[Gene],values)
                            tup1=np.array(finalmat)
                            tup2=tup1.mean(axis=0)
                            figurelist.append(tup2)
                            for index in range(len(values)):
                                print(legendname[Bed],str(index+1),*values[index],sep="\t",end="\n",file=out_matrix)
                print("region","model","value",sep="\t",end="\n",file=out_profile)
                N = 0
                for lis in figurelist:
                    n = 0
                    for l in lis:
                        print(n*wdsize,legendname[N],l,sep="\t",end="\n",file=out_profile)
                        n += 1
                    N += 1
                for line in figurelist:
                    plotsmooth(line)
                plt.legend(legendname)
                plt.ylabel("Normalized reads count")
                plt.xticks([left,medium,right],["-"+L+"kb","TSS",R+"kb"])
                plt.savefig(out_figure)     
                out_figure.close()
                out_profile.close()
                out_matrix.close()
        elif referencePoint=="end":
            if genelist is None:
                out_profile=open(outf[0],"w")
                out_matrix=open(outf[1],"w")
                out_figure=open(figurename,"w")
                figurelist=[]
                header = ["label","num"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
                print(*header,sep="\t",end="\n",file=out_matrix)
                for Bam in range(len(bamlist)):
                    for Bed in range(len(bedlist)):
                        values=bamTocountEnd(bedlist[Bed],bamlist[Bam],wdsize,up,down)
                        tup1=np.array(values)
                        tup2=tup1.mean(axis=0)
                        figurelist.append(tup2)
                        for index in range(len(values)):
                            print(legendname[Bed],str(index+1),*values[index],sep="\t",end="\n",file=out_matrix)
                print("region","model","value",sep="\t",end="\n",file=out_profile)
                N = 0
                for lis in figurelist:
                    n = 0
                    for l in lis:
                        print(n*wdsize,legendname[N],l,sep="\t",end="\n",file=out_profile)
                        n += 1
                    N += 1
                for line in figurelist:
                    plotsmooth(line)
                plt.legend(legendname)
                plt.ylabel("Normalized reads count")
                plt.xticks([left,medium,right],["-"+L+"kb","TSS",R+"kb"])
                plt.savefig(out_figure)     
                out_figure.close()
                out_profile.close()
                out_matrix.close()
            else:
                out_profile=open(outf[0],"w")
                out_matrix=open(outf[1],"w")
                out_figure=open(figurename,"w")
                figurelist=[]
                header = ["label","num"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
                print(*header,sep="\t",end="\n",file=out_matrix)
                for Bam in range(len(bamlist)):
                    for Bed in range(len(bedlist)):
                        for Gene in range(len(genelist)):
                            values=bamtocountEnd(bedlist[Bed],bamlist[Bam],wdsize,up,down)
                            finalmat=Listgene(genelist[Gene],values)
                            tup1=np.array(finalmat)
                            tup2=tup1.mean(axis=0)
                            figurelist.append(tup2)
                            for index in range(len(values)):
                                print(legendname[Bed],str(index+1),*values[index],sep="\t",end="\n",file=out_matrix)
                print("region","model","value",sep="\t",end="\n",file=out_profile)
                N = 0
                for lis in figurelist:
                    n = 0
                    for l in lis:
                        print(n*wdsize,legendname[N],l,sep="\t",end="\n",file=out_profile)
                        n += 1
                    N += 1
                for line in figurelist:
                    plotsmooth(line)
                plt.legend(legendname)
                plt.ylabel("Normalized reads count")
                plt.xticks([left,medium,right],["-"+L+"kb","TSS",R+"kb"])
                plt.savefig(out_figure)     
                out_figure.close()
                out_profile.close()
                out_matrix.close()
        else:
            pass
	#print("Usage: python tmp1.py <<--scale-regions or --reference-point>> -b <bamfile> -G <bedfile> -g <genelist> -r <referencePoint {start or end}> -u <upstream_distance> -d <downstream_distance> -w <window_size> -bn <Equal_block_number> -o <outputfile>")

    elif scaleRegion:
        left=int(0)
        mediumA=int(up/wdsize)
        mediumB=mediumA+blocknum
        right=int(down/wdsize)+mediumA+blocknum
        NumOfBin = int(up/wdsize)+blocknum+int(down/wdsize)
        if genelist is None:
            out_profile=open(outf[0],"w")
            out_matrix=open(outf[1],"w")
            out_figure=open(figurename,"w")
            figurelist=[]
            header = ["label","num"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
            print(*header,sep="\t",end="\n",file=out_matrix)
            for Bam in range(len(bamlist)):
                for Bed in range(len(bedlist)):
                    values=bamtobodycount(bedlist[Bed],bamlist[Bam],wdsize,up,down,blocknum)
                    tup1=np.array(values)
                    tup2=tup1.mean(axis=0)
                    figurelist.append(tup2)
                    for index in range(len(values)):
                        print(legendname[Bed],str(index+1),*values[index],sep="\t",end="\n",file=out_matrix)
            print("region","model","value",sep="\t",end="\n",file=out_profile)
            N = 0
            for lis in figurelist:
                n = 0
                for l in lis:
                    print(n*wdsize,legendname[N],l,sep="\t",end="\n",file=out_profile)
                    n += 1
                N += 1
            for line in figurelist:
                plotsmooth(line)
            plt.legend(legendname)
            plt.ylabel("Normalized reads count")
            plt.xticks([left,mediumA,mediumB,right],["-"+L+"kb","start","end",R+"kb"])
            plt.savefig(out_figure)     
            out_figure.close()
            out_profile.close()
            out_matrix.close()
        else:
            out_profile=open(outf[0],"w")
            out_matrix=open(outf[1],"w")
            out_figure=open(figurename,"w")
            figurelist=[]
            header = ["label","num"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
            print(*header,sep="\t",end="\n",file=out_matrix)
            for Bam in range(len(bamlist)):
                for Bed in range(len(bedlist)):
                    for Gene in range(len(genelist)):
                        values=bamGenometobodycount(bedlist[Bed],bamlist[Bam],wdsize,up,down,blocknum)
                        finalmat=Listgene(genelist[Gene],values)
                        tup1=np.array(finalmat)
                        tup2=tup1.mean(axis=0)
                        figurelist.append(tup2)
                        for index in range(len(values)):
                            print(legendname[Bed],str(index+1),*values[index],sep="\t",end="\n",file=out_matrix)
            print("region","model","value",sep="\t",end="\n",file=out_profile)
            N = 0
            for lis in figurelist:
                n = 0
                for l in lis:
                    print(n*wdsize,legendname[N],l,sep="\t",end="\n",file=out_profile)
                    n += 1
                N += 1
            for line in figurelist:
                plotsmooth(line)
            plt.legend(legendname)
            plt.ylabel("Normalized reads count")
            plt.xticks([left,medium1,medium2,right],["-"+L+"kb","start","end",R+"kb"])
            plt.savefig(out_figure)     
            out_figure.close()
            out_profile.close()
            out_matrix.close()
    else:
        pass
