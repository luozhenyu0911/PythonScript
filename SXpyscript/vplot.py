# -*- coding: UTF-8 -*-
from __future__ import print_function
from __future__ import division
import argparse
import pandas as pd
import numpy as np
import pysam
from math import sqrt,ceil
import matplotlib.pyplot as plt
import seaborn as sns
#from scipy.interpolate import interp1d
from scipy.interpolate import make_interp_spline, BSpline
plt.switch_backend('PDF')
plt.rcParams['pdf.fonttype'] = 42
##@jit

example_text = '''example:
 2020.6.14
 python **.py --referencePoint start/end -b **.bam,**.bam... -G **.bed,**.bed... -L bamlabel1,bamlabel2 -l bedlabel1,bedlabel2 -u 1000 -d 1000 -w 50 -o *vplot.matrix -f figure.pdf
'''
parser = argparse.ArgumentParser(description="This tool creates scores and profile plot over sets of genomic regions. Typically, these regions are genes, but any other regions defined in BED will work.", 
                                 formatter_class= argparse.RawTextHelpFormatter, #用于自定义帮助文档输出格式的类
                                 #prog='python *py',                             #显示的程序名称
                                 usage='%(prog)s [-h]',                          #格式化显示方式(描述程序用途的字符串)
                                 epilog=example_text)                            #在参数帮助文档之后显示的文本（默认值：无）
group=parser.add_mutually_exclusive_group(required= True)
group.add_argument('--referencePoint',type=str,help="<referencePoint {start or end}>"+"\t"+"plot the distribution of the reads at reference point from the genome. This option is mutually exclusive with \"--scaleRegion\".",metavar='')
group.add_argument('--scaleRegion',action='store_true',help="<Scaled Region>"+"\t"+"plot the distribution of the reads over the region from the genome. This option is mutually exclusive with \"--referencePoint\".")

parser.add_argument('--bam','-b',type= str,help="<bam or sam files>"+"\t"+"data files(Support for multiple files, Separated by commas), should be sorted format.(Note: Need to filter out low-quality reads first).",required= True,metavar='')
parser.add_argument('--format',type= str,help="<bam/sam>"+"\t"+"Specify Alignment file format. (Default=bam)",default="bam",metavar='')

group2=parser.add_mutually_exclusive_group(required= True)
group2.add_argument('--bed','-G',type=str,help="<bedfiles>"+"\t"+"data files(Support for multiple files, Separated by commas), should be original bed files.",metavar='')
group2.add_argument('--bedList',type= str,help="<bedfiles>"+"\t"+"bed file text file.",metavar='')

parser.add_argument('--upstream','-u',type=int,help="<upstream distance>"+"\t"+"Distance upstream of the reference-point selected.",required= True,metavar='')
parser.add_argument('--downstream','-d',type=int,help="<downstream distance>"+"\t"+"Distance downstream of the reference-point selected.",required= True,metavar='')
parser.add_argument('--windows','-w',type=int,help="<window size>"+"\t"+"Length, in bases, binSize for averaging the score over the regions length.",required= True,metavar='')
parser.add_argument('--Maxfragment','-M',type=int,help="<Max fragment size>"+"\t"+"Maximum fragment size to be plot.",required= True,metavar='')

parser.add_argument('--bamlabel','-L',type=str,help="<bam legend names>"+"\t"+"Labels of Bam file. (Support for multiple legends).",required= True,metavar='')

group3=parser.add_mutually_exclusive_group(required= True)
group3.add_argument('--bedlabel','-l',type=str,help="<bed legend names>"+"\t"+"Labels of bed file. (Support for multiple legends).",metavar='')
group3.add_argument('--bedlabeltxt',type=str,help="<bed legend names>"+"\t"+"A text file containg all labels of bed files. Usually used when you have a large group of bed files.",metavar='')

#parser.add_argument('--Facet',type=str,help="<Bam/Bed>"+"\t"+"If facet by bam file, each bam file as a subplot to draw all bedfiles. If facet by bed file, each bed site as a subplot to draw all bamfiles reads.",required= True,metavar='')
#parser.add_argument('--scaleY',action='store_true',help="<Optional argument>"+"\t"+"Whether or not to scale Y limits to the same. (Default=False)",default=False)

parser.add_argument('--outfile','-o',type=str,help="<outputfiles final matrix>"+"\t"+"File name to save the average data file and final matrix file. (Separated by commas).",required= True,metavar='')
parser.add_argument('--figure','-f',type=str,help="<Figure name>"+"\t"+"Name of the plot name. (plot file format is PDF).",required= True,metavar='')
args = parser.parse_args()

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
maxf=args.Maxfragment

Bamlabels=args.bamlabel.split(",")
#Bedlabels=args.bedlabel.split(",")
if args.bedlabel:
    Bedlabels=args.bedlabel.split(",")
else:
    Bedlabels=[]
    with open(args.bedlabeltxt,'r') as f:
        for line in f:
            line=line.strip().split()
            Bedlabels.append(line[0])

#Facet=args.Facet

figurename=args.figure
outf=args.outfile

format={'bam':'rb','sam':'r'}

#in: bam
#out: [read1:[chr,start,end,length],read2:[chr,start,end,length]]
def bamTobed(Input):
    ReadInfo={}
    for read in Input.fetch():
        ID=read.query_name
        if ID in ReadInfo:
            if int(read.reference_start) < int(ReadInfo[ID][1]):
                ReadInfo[ID][1]=int(read.reference_start)
                ReadInfo[ID][3]=int(ReadInfo[ID][2])-int(ReadInfo[ID][1])+1
            else:
                ReadInfo[ID][2]=int(read.reference_end)
                ReadInfo[ID][3]=int(ReadInfo[ID][2])-int(ReadInfo[ID][1])+1
        else:
            readInfo=[]
            readInfo.append(str(read.reference_name))
            readInfo.append(int(read.reference_start))
            readInfo.append(int(read.reference_end))
            readInfo.append(int(readInfo[2])-int(readInfo[1])+1)
            ReadInfo[ID]=readInfo
    return ReadInfo

#in: see below out
#out: {50:{chr:{'1000':2,'2000':3}},
#      51:{chr:{'1000':2,'2000':3}}}
def bedTodict(bed):
    region={}
    for value in bed.values():
        start=value[1]
        end=value[2]
        length=value[3]
        chrom=value[0]
        remain=length%2
        if length not in region:
            region[length]={}
        if chrom not in region[length]:
            region[length][chrom]={}
        #偶数长度的中点有两个，长度n-1的中点和长度n+1的中点
        if int(remain) == 0:
            centerPoint1=int((start+end-1)/2)
            centerPoint2=int((start+end+1)/2)
            if str(centerPoint1) in region[length][chrom]:
                region[length][chrom][str(centerPoint1)]+=0.5
            else:
                region[length][chrom][str(centerPoint1)]=0.5
            if str(centerPoint2) in region[length][chrom]:
                region[length][chrom][str(centerPoint2)]+=0.5
            else:
                region[length][chrom][str(centerPoint2)]=0.5
        #奇数长度的中点等于首尾之和的一半
        else:
            centerPoint=int((start+end)/2)
            if str(centerPoint) in region[length][chrom]:
                region[length][chrom][str(centerPoint)]+=1
            else:
                region[length][chrom][str(centerPoint)]=1
    return region

#Lis [[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],...]
#regard coordinate start as centerpoint
def bedCount(bedf,readdict,count,wdsize,up,down):
    bed6col=open(bedf,"r")
    Lis=[]
    sumReads=0
    for bed in bed6col:
        row2=bed.rstrip().split()
        chr=row2[0]
        sublist=[]
        if "+" in row2[5]:
            left=int(row2[1])-up
            right=int(row2[1])+down
            #如果字典中有该染色体的reads,则计算
            if chr in readdict:
                for binStart in range(int(left),int(right),wdsize):
                    for j in range(binStart,(binStart+wdsize)):
                        if str(j) in readdict[chr]:
                            sumReads+=readdict[chr][str(j)]
                        else:
                            pass
                    normalizedreads = (sumReads*1000000000)/(count*wdsize)
                    sublist.append(normalizedreads)
                    sumReads=0
                Lis.append(sublist)
            else:
                misvalue=[0 for _ in range(int(left),int(right),wdsize)]
                Lis.append(misvalue)
        else:
            left=int(row2[2])-down
            right=int(row2[2])+up
            #如果字典中有该染色体的reads,则计算
            if chr in readdict:
                for binStart in range(int(right),int(left),-wdsize):
                    for j in range(binStart,(binStart-wdsize),-1):
                        if str(j) in readdict[chr]:
                            sumReads+=readdict[chr][str(j)]
                        else:
                            pass
                    sublist.append((sumReads*1000000000)/(count*wdsize))
                    sumReads=0
                Lis.append(sublist)
            else:
                misvalue=[0 for _ in range(int(right),int(left),-wdsize)]
                Lis.append(misvalue)
    return Lis
    bed6col.close()

#Lis [[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],...]
#regard coordinate end as centerpoint
def bedCountEnd(bedf,readdict,count,wdsize,up,down):
    bed6col=open(bedf,"r")
    Lis=[]
    sumReads=0
    for bed in bed6col:
        row2=bed.rstrip().split()
        chr=row2[0]
        sublist=[]
        if "+" in row2[5]:
            left=int(row2[2])-up
            right=int(row2[2])+down
            for k in range(int(left),int(right),wdsize):
                for j in range(k,(k+wdsize)):
                    #如果字典中有该染色体的reads,则计算
                    if chr in readdict:
                        if str(j) in readdict[chr]:
                            sumReads+=readdict[chr][str(j)]
                        else:
                            pass
                    else:#如果字典中没有该染色体上的reads，则跳过
                        pass
                sublist.append((sumReads*1000000000)/(count*wdsize))
                sumReads=0
            Lis.append(sublist)
        else:
            left=int(row2[1])-down
            right=int(row2[1])+up
            for k in range(int(right),int(left),-wdsize):
                for j in range(k,(k-wdsize),-1):
                    #如果字典中有该染色体的reads,则计算
                    if chr in readdict:
                        if str(j) in readdict[chr]:
                            sumReads+=readdict[chr][str(j)]
                        else:
                            pass
                    else:#如果字典中没有该染色体上的reads，则跳过
                        pass
                sublist.append((sumReads*1000000000)/(count*wdsize))
                sumReads=0
            Lis.append(sublist)
    return Lis
    bed6col.close()

def vplot(df,xtick,xticklabel,ylabel,outfigname):
    mat = df.iloc[:,3:]
    array=np.array(mat)
    q=np.quantile(array,0.99999)
    array[array>q]=q
    mat=pd.DataFrame(array)
    mat.columns = np.array(range(mat.shape[1]))+1
    ax = sns.heatmap(mat)
    ax.invert_yaxis()
    ax.set_xticks(xtick)
    ax.set_xticklabels(xticklabel,rotation=0)
    ax.set_ylabel(ylabel)
    plt.yticks(rotation=0)
    plt.xticks(rotation=0)
    #ax.set(ylim=ylim)
    ax.figure.set_size_inches(5,4)
    plt.savefig(outfigname)

if referencePoint:
    ylength=int(int(up+down)/wdsize)
    left=int(1)
    medium=int(ylength/2)+0.5
    right=int(ylength)
    NumOfBin = int((int(up+down))/wdsize)
    if referencePoint=="start":
        out_matrix=open(outf,"w")
        out_figure=open(figurename,"w")
        #header
        header = ["Bamlabel","Bedlabel","length"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
        print(*header,sep="\t",end="\n",file=out_matrix)
        for Bam in range(len(bamlist)):
            bamf=pysam.AlignmentFile(bamlist[Bam],format[args.format])
            bamTotalReads=bamf.count()
            bamReadsbed=bamTobed(bamf)
            bamrCoverageDict=bedTodict(bamReadsbed)
            for Bed in range(len(bedlist)):
                for length in range(0,maxf):
                    if length in bamrCoverageDict:
                        values=bedCount(bedlist[Bed],bamrCoverageDict[length],bamTotalReads,wdsize,up,down)
                        tup1=np.array(values)
                        tup2=tup1.mean(axis=0)
                        #print average bin value
                        print(Bamlabels[Bam],Bedlabels[Bed],str(length),*list(tup2),sep="\t",end="\n",file=out_matrix)
                    else:
                        meanvalue=[0 for _ in range(NumOfBin)]
                        print(Bamlabels[Bam],Bedlabels[Bed],str(length),*meanvalue,sep="\t",end="\n",file=out_matrix)
        out_matrix.close()
        out_figure.close()
        dataframe=pd.read_csv(outf,sep='\t')
        Xticks=[left,medium,right]
        Xticklabel=["-"+L+"kb","TSS",R+"kb"]
        Ylabel="Fragment length"
        vplot(dataframe,Xticks,Xticklabel,Ylabel,figurename)
    elif referencePoint=="end":
        out_matrix=open(outf,"w")
        out_figure=open(figurename,"w")
        header = ["Bamlabel","Bedlabel","length"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
        print(*header,sep="\t",end="\n",file=out_matrix)
        for Bam in range(len(bamlist)):
            bamf=pysam.AlignmentFile(bamlist[Bam],format[args.format])
            bamTotalReads=bamf.count()
            bamReadsbed=bamTobed(bamf)
            bamrCoverageDict=bedTodict(bamReadsbed)
            for Bed in range(len(bedlist)):
                for length in range(0,maxf):
                    if length in bamrCoverageDict:
                        values=bedCountEnd(bedlist[Bed],bamrCoverageDict[length],bamTotalReads,wdsize,up,down)
                        tup1=np.array(values)
                        tup2=tup1.mean(axis=0)
                        #print average bin value
                        print(Bamlabels[Bam],Bedlabels[Bed],str(length),*list(tup2),sep="\t",end="\n",file=out_matrix)
                    else:
                        meanvalue=[0 for _ in range(NumOfBin)]
                        print(Bamlabels[Bam],Bedlabels[Bed],str(length),*meanvalue,sep="\t",end="\n",file=out_matrix)
        out_matrix.close()
        out_figure.close()
        dataframe=pd.read_csv(outf,sep='\t')
        Xticks=[left,medium,right]
        Xticklabel=["-"+L+"kb","TSS",R+"kb"]
        Ylabel="Fragment length"
        vplot(dataframe,Xticks,Xticklabel,Ylabel,figurename)
    else:
        pass
