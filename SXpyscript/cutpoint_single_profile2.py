# -*- coding: UTF-8 -*-
from __future__ import print_function
from __future__ import division
#from numba import jit
#import sys, getopt
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
 python **.py --referencePoint start/end/center -b **.bam,**.bam... -G **.bed,**.bed... -L bamlabel1,bamlabel2 -l bedlabel1,bedlabel2 --Facet Bam/Bed -u 1000 -d 1000 -w 50 -o average.data,matrix.data -f figure.pdf
 python **.py --scaleRegion -b **.bam,**bam... -G **.bed,**.bed... -L bamlabel1,bamlabel2 -l bedlabel1,bedlabel2 --Facet Bam/Bed -u 1000 -d 1000 -w 50 -N 100 -o average.data,matrix.data -f figure.pdf
'''
parser = argparse.ArgumentParser(description="This tool creates scores and profile plot over sets of genomic regions. Typically, these regions are genes, but any other regions defined in BED will work.", 
                                 formatter_class= argparse.RawTextHelpFormatter,
                                 prog='base_maker',
                                 epilog=example_text)
group=parser.add_mutually_exclusive_group(required= True)
group.add_argument('--referencePoint',type=str,help="<referencePoint {start/end/center}>"+"\t"+"plot the distribution of the reads at reference point from the genome.")
group.add_argument('--scaleRegion',action='store_true',help="<Scaled Region>"+"\t"+"plot the distribution of the reads over the region from the genome.")

parser.add_argument('--control', '-c',type= str,help="<Input bam files>"+"\t"+"Control data files(Support for multiple files, Separated by commas, the order is consistent with bam files), should be sorted bam format.")
parser.add_argument('--bam','-b',type= str,help="<bam or sam files>"+"\t"+"data files(Support for multiple files, Separated by commas), should be sorted format (Note: Need to filter out low-quality reads first).")
parser.add_argument('--format',type= str,help="<bam/sam>"+"\t"+"Specify Alignment file format. (Default=bam)",default="bam")

group2=parser.add_mutually_exclusive_group(required= True)
group2.add_argument('--bed','-G',type=str,help="<bedfiles>"+"\t"+"data files(Support for multiple files, Separated by commas), should be original bed files.")
group2.add_argument('--bedList',type= str,help="<bedfiles>"+"\t"+"bed file text file.")

parser.add_argument('--upstream','-u',type=int,help="<upstream distance>"+"\t"+"Distance upstream of the reference-point selected.",required= True)
parser.add_argument('--downstream','-d',type=int,help="<downstream distance>"+"\t"+"Distance downstream of the reference-point selected.",required= True)
parser.add_argument('--windows','-w',type=int,help="<window size>"+"\t"+"Length, in bases, binSize for averaging the score over the regions length.",required= True)
parser.add_argument('--Nbin','-N',type=int,help="<Number of bins>"+"\t"+"The number of equally divided bins, ONLY works when used with <--scaleRegion>")

parser.add_argument('--bamlabel','-L',type=str,help="<bam legend names>"+"\t"+"Labels of Bam file. (Support for multiple legends).",required= True)

group3=parser.add_mutually_exclusive_group(required= True)
group3.add_argument('--bedlabel','-l',type=str,help="<bed legend names>"+"\t"+"Labels of bed file. (Support for multiple legends).")
group3.add_argument('--bedlabeltxt',type=str,help="<bed legend names>"+"\t"+"A text file containg all labels of bed files. Usually used when you have a large group of bed files.")

parser.add_argument('--Facet',type=str,help="<Bam/Bed>"+"\t"+"If facet by bam file, each bam file as a subplot to draw all bedfiles. If facet by bed file, each bed site as a subplot to draw all bamfiles reads.",required= True)
parser.add_argument('--scaleY',action='store_true',help="<Optional argument>"+"\t"+"Whether or not to scale Y limits to the same. (Default=False)",default=False)

parser.add_argument('--outfile','-o',type=str,help="<outputfiles average data,final matrix>"+"\t"+"File name to save the average data file and final matrix file. (Separated by commas).",required= True)
parser.add_argument('--figure','-f',type=str,help="<Figure name>"+"\t"+"Name of the plot name. (plot file format is PDF).",required= True)
args = parser.parse_args()

if args.control:
    Inputbamlist=args.control.split(",")
else:
    Inputbamlist=None

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

Facet=args.Facet

figurename=args.figure
outf=args.outfile.split(",")

format={'bam':'rb','sam':'r'}

def bamTodict(bf):
    region={}
    for read in bf.fetch():
        a=[]
        a.append(read.reference_name)
        a.append(read.reference_start)
        a.append(read.reference_end)
        chrom=str(a[0])
        if str(read.is_reverse) == "False":
            if chrom in region:
                pass
            else:
                region[chrom]={}
            if str(a[1]) in region[chrom]:
                region[chrom][str(a[1])]+=1
            else:
                region[chrom][str(a[1])]=1
        elif str(read.is_reverse) == "True":
            if chrom in region:
                pass
            else:
                region[chrom]={}
            if str(a[2]) in region[chrom]:
                region[chrom][str(a[2])]+=1
            else:
                region[chrom][str(a[2])]=1
        else:
            pass
    return region

#Lis [[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],...]
#regard coordinate start as centerpoint
def bedCount(bedf,readdict,count,wdsize,up,down):
    bed6col=open(bedf,"r")
    Lis=[]
    sumReads=0  
    for bed in bed6col:
        row2=bed.rstrip().split()
        sublist=[]
        if "+" in row2[5]:
            left=int(row2[1])-up
            right=int(row2[1])+down
            for k in range(int(left),int(right),wdsize):
                for j in range(k,(k+wdsize)):
                    #此处默认了row2[0]染色体已经在字典中了
                    if str(j) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(j)])
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
                    if str(j) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(j)])
                    else:
                        pass
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            Lis.append(sublist)
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
        sublist=[]
        if "+" in row2[5]:
            left=int(row2[2])-up
            right=int(row2[2])+down
            for k in range(int(left),int(right),wdsize):
                for j in range(k,(k+wdsize)):
                    if str(j) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(j)])
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
                    if str(j) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(j)])
                    else:
                        pass
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            Lis.append(sublist)
    return Lis
    bed6col.close()

#2020.12.23
def bedCountCenter(bedf,readdict,count,wdsize,up,down):
    bed6col=open(bedf,"r")
    Lis=[]
    sumReads=0
    for bed in bed6col:
        row2=bed.rstrip().split()
        chr=row2[0]
        sublist=[]
        if "+" in row2[5]:
            mid=(int(row2[1])+int(row2[2]))//2
            left=mid-up
            right=mid+down
            for binStart in range(int(left),int(right),wdsize):
                for j in range(binStart,(binStart+wdsize)):
                    #此处默认了chr染色体已经在字典中了
                    if str(j) in readdict[chr]:
                        sumReads+=readdict[chr][str(j)]
                    else:
                        pass
                sublist.append((sumReads*10000000)/(count*wdsize))
                sumReads=0
            Lis.append(sublist)
        else:
            mid=(int(row2[1])+int(row2[2]))//2
            left=mid-down
            right=mid+up
            for binStart in range(int(right),int(left),-wdsize):
                for j in range(binStart,(binStart-wdsize),-1):
                    if str(j) in readdict[chr]:
                        sumReads+=readdict[chr][str(j)]
                    else:
                        pass
                sublist.append((sumReads*10000000)/(count*wdsize))
                sumReads=0
            Lis.append(sublist)
    return Lis
    bed6col.close()

#Lis [[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8]]
# up start end down
def BodyCount(bedf,readdict,count,wdsize,up,down,BlockNum):
    bed6col=open(bedf,"r")
    Lis=[]
    sumReads=0
    for bed in bed6col:
        row2=bed.rstrip().split()
        sublist=[]
       # blockSize=((int(row2[2])-int(row2[1]))/BlockNum)+1
       # reside=BlockNum-((int(row2[2])-int(row2[1]))%BlockNum)
        if int(((int(row2[2])-int(row2[1]))%BlockNum)) == 0:
            blockSize = int((int(row2[2])-int(row2[1]))/BlockNum)
            reside=0
        else:
            reside=BlockNum-int(((int(row2[2])-int(row2[1]))%BlockNum))
            blockSize=int(((int(row2[2])-int(row2[1]))+reside)/BlockNum)
        if "+" in row2[5]:
            for k in range(int(row2[1])-up,int(row2[1]),wdsize):
                for j in range(k,(k+wdsize)):
                    if str(j) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(j)])
                    else:
                        pass
               # sumReads=int((sumReads/10)*10)
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            for x in range(int(row2[1]),int(row2[2])+reside,blockSize):
                for y in range(x,(x+blockSize)):
                    if str(y) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(y)])
                    else:
                        pass
                #sumReads=int((sumReads/blockSize)*10)
                sublist.append((int(sumReads)*10000000)/(count*blockSize))
                sumReads=0
            for u in range(int(row2[2]),int(row2[2])+down,wdsize):
                for v in range(u,(u+wdsize)):
                    if str(v) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(v)])
                    else:
                        pass
               # sumReads=int((sumReads/10)*10)
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            Lis.append(sublist)
        else:
            for k in range(int(row2[2])+up,int(row2[2]),-wdsize):
                for j in range(k,(k-wdsize),-1):
                    if str(j) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(j)])
                    else:
                        pass
              #  sumReads=int((sumReads/10)*10)
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            for x in range(int(row2[2]),int(row2[1])-reside,-blockSize):
                for y in range(x,(x-blockSize),-1):
                    if str(y) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(y)])
                    else:
                        pass
              #  sumReads=int((sumReads/blockSize)*10)
                sublist.append((int(sumReads)*10000000)/(count*blockSize))
                sumReads=0
            for u in range(int(row2[1]),int(row2[1])-down,-wdsize):
                for v in range(u,(u-wdsize),-1):
                    if str(v) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(v)])
                    else:
                        pass
              #  sumReads=int((sumReads/10)*10)
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            Lis.append(sublist)
    return Lis
    bed6col.close()

# def plotsmooth(Values):
#     n=len(Values)
#     num=range(0,n)
#     num=np.array(num)
#     xnew=np.linspace(num.min(),num.max(),5*n)
#     func=interp1d(num,Values,kind="cubic")
#     ynew=func(xnew)
#     plt.plot(xnew,ynew)

def plotsmooth(Values):
    n=len(Values)
    num=range(0,n)
    x=np.array(num)
    xnew=np.linspace(x.min(),x.max(),5*n)
    func=make_interp_spline(x, Values, k=3)
    ynew=func(xnew)
    plt.plot(xnew,ynew)

def plotfacet(facet,df,xtick,xticklabel,ylabel,outfigname,scaleY):
    Pal="Set1"
    if facet=="Bam":
        u = df.Bamlabel.unique()
        #设定分面布局行列数
        nfacet=len(u)
        if nfacet >= 4:
            ncol=int(ceil(sqrt(nfacet)))
            nrow=int(round(sqrt(nfacet)))
        else:
            ncol=nfacet
            nrow=1
        #判断有多个分面还是1个
        if len(u)>1:
            fig, axes = plt.subplots(nrows=nrow,ncols=ncol)
            axes=axes.flatten()
            if scaleY:
                max_value=max(df.value)
                min_value=min(df.value)
                max_y=max_value+(max_value-min_value)*0.05
                min_y=max(0,min_value-(max_value-min_value)*0.05)
                custom_ylim = (min_y,max_y)
                # Setting the values for all axes.
                plt.setp(axes,ylim=custom_ylim)
            #遍历每个分面
            for Bamlabel, ax in zip(u, axes):
                subdf=df[df.Bamlabel==Bamlabel]
                subdf_new=pd.DataFrame()
                #列表储存插值后的新标签和值
                newbin=[]
                newvalue=[]
                newBedlabel=[]
                #遍历每条线
                for Bedlabel in subdf.Bedlabel.unique():
                    n=len(subdf[subdf.Bedlabel==Bedlabel])
                    x=subdf.loc[subdf.Bedlabel==Bedlabel,'bin']
                    y=subdf.loc[subdf.Bedlabel==Bedlabel,'value']
                    xnew=np.linspace(x.min(),x.max(),5*n)
                    func=make_interp_spline(x, y, k=3)
                    ynew=func(xnew)
                    newBedlabel+=[Bedlabel]*n*5
                    newbin=newbin+list(xnew)
                    newvalue=newvalue+list(ynew)
                #插值后赋值给一个新Dataframe
                subdf_new['Bedlabel']=newBedlabel
                subdf_new['bin']=newbin
                subdf_new['value']=newvalue
                #开始画图
                pal=sns.color_palette(Pal,len(subdf_new['Bedlabel'].unique()))
                g=sns.lineplot(x="bin",y='value',hue='Bedlabel',ax =ax ,data=subdf_new,palette=pal)
                #保留左边和底部坐标轴
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.yaxis.set_ticks_position('left')
                ax.xaxis.set_ticks_position('bottom')
                #移除图例标题
                handles, labels = ax.get_legend_handles_labels()
                ax.legend(handles=handles[1:], labels=labels[1:])
                g.set(title=Bamlabel)
                g.set_xticks(xtick)
                g.set_xticklabels(xticklabel)
                g.set_xlabel("")
                g.set_ylabel(ylabel)
                g.figure.set_size_inches(ncol*4.5,nrow*3.5)
            #plt.tight_layout()
            fig.set_tight_layout(True)
        else:
            fig, axes = plt.subplots(ncols=1)
            if scaleY:
                max_value=max(df.value)
                min_value=min(df.value)
                max_y=max_value+(max_value-min_value)*0.05
                min_y=max(0,min_value-(max_value-min_value)*0.05)
                custom_ylim = (min_y,max_y)
                # Setting the values for all axes.
                plt.setp(axes,ylim=custom_ylim)
            Bamlabel=u[0]
            df_new=pd.DataFrame()
            newbin=[]
            newvalue=[]
            newBedlabel=[]
            for Bedlabel in df.Bedlabel.unique():
                n=len(df[df.Bedlabel==Bedlabel])
                x=df.loc[df.Bedlabel==Bedlabel,'bin']
                y=df.loc[df.Bedlabel==Bedlabel,'value']
                xnew=np.linspace(x.min(),x.max(),5*n)
                func=make_interp_spline(x, y, k=3)
                ynew=func(xnew)
                newBedlabel+=[Bedlabel]*n*5
                newbin=newbin+list(xnew)
                newvalue=newvalue+list(ynew)
            #插值后赋值给一个新Dataframe
            df_new['Bedlabel']=newBedlabel
            df_new['bin']=newbin
            df_new['value']=newvalue
            #开始画图
            pal=sns.color_palette(Pal,len(df_new['Bedlabel'].unique()))
            g=sns.lineplot(x="bin",y='value',hue='Bedlabel',ax = axes,data=df_new,palette=pal)
            #保留左边和底部坐标轴
            axes.spines['right'].set_visible(False)
            axes.spines['top'].set_visible(False)
            axes.yaxis.set_ticks_position('left')
            axes.xaxis.set_ticks_position('bottom')
            #移除图例标题
            handles, labels = axes.get_legend_handles_labels()
            axes.legend(handles=handles[1:], labels=labels[1:])
            #设置标签
            g.set(title=Bamlabel)
            g.set_xticks(xtick)
            g.set_xticklabels(xticklabel)
            g.set_xlabel("")
            g.set_ylabel(ylabel)
            g.figure.set_size_inches(4.5,3.5)
            #plt.tight_layout()
            fig.set_tight_layout(True)
        plt.savefig(outfigname)
    elif facet=="Bed":
        #按照bed分面，每个图是一个bed,1个bed就是一个子图，多个bed就是多个子图; 
        #bam文件可以是1个或多个，1个bam一条线，多个bam多条线
        u = df.Bedlabel.unique()
        #设定分面布局行列数
        nfacet=len(u)
        if nfacet >= 4:
            ncol=int(ceil(sqrt(nfacet)))
            nrow=int(round(sqrt(nfacet)))
        else:
            ncol=nfacet
            nrow=1
        #判断有多个分面还是1个
        if len(u)>1:
            fig, axes = plt.subplots(nrows=nrow,ncols=ncol)
            #这一行很重要，当subplot不是一维的时候，axes是一个二维数组
            axes=axes.flatten()
            if scaleY:
                max_value=max(df.value)
                min_value=min(df.value)
                max_y=max_value+(max_value-min_value)*0.05
                min_y=max(0,min_value-(max_value-min_value)*0.05)
                custom_ylim = (min_y,max_y)
                # Setting the values for all axes.
                plt.setp(axes,ylim=custom_ylim)
            #遍历每个分面
            for Bedlabel, ax in zip(u, axes):
                subdf=df[df.Bedlabel==Bedlabel]
                subdf_new=pd.DataFrame()
                #列表储存插值后的新标签和值
                newbin=[]
                newvalue=[]
                newBamlabel=[]
                #遍历每条线
                for Bamlabel in subdf.Bamlabel.unique():
                    n=len(subdf[subdf.Bamlabel==Bamlabel])
                    x=subdf.loc[subdf.Bamlabel==Bamlabel,'bin']
                    y=subdf.loc[subdf.Bamlabel==Bamlabel,'value']
                    xnew=np.linspace(x.min(),x.max(),5*n)
                    func=make_interp_spline(x, y, k=3)
                    ynew=func(xnew)
                    newBamlabel+=[Bamlabel]*n*5
                    newbin=newbin+list(xnew)
                    newvalue=newvalue+list(ynew)
                #插值后赋值给一个新Dataframe
                subdf_new['Bamlabel']=newBamlabel
                subdf_new['bin']=newbin
                subdf_new['value']=newvalue
                #开始画图
                pal=sns.color_palette(Pal,len(subdf_new['Bamlabel'].unique()))
                g=sns.lineplot(x="bin",y='value',hue='Bamlabel',ax =ax ,data=subdf_new,palette=pal)
                #保留左边和底部坐标轴
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.yaxis.set_ticks_position('left')
                ax.xaxis.set_ticks_position('bottom')
                #移除图例标题
                handles, labels = ax.get_legend_handles_labels()
                ax.legend(handles=handles[1:], labels=labels[1:])
                #设置标签
                g.set(title=Bedlabel)
                g.set_xticks(xtick)
                g.set_xticklabels(xticklabel)
                g.set_xlabel("")
                g.set_ylabel(ylabel)
                g.figure.set_size_inches(ncol*4.5,nrow*3.5)
            #plt.tight_layout()
            fig.set_tight_layout(True)
        else:
            fig, axes = plt.subplots(ncols=1)
            if scaleY:
                max_value=max(df.value)
                min_value=min(df.value)
                max_y=max_value+(max_value-min_value)*0.05
                min_y=max(0,min_value-(max_value-min_value)*0.05)
                custom_ylim = (min_y,max_y)
                # Setting the values for all axes.
                plt.setp(axes,ylim=custom_ylim)
            Bedlabel=u[0]
            df_new=pd.DataFrame()
            newbin=[]
            newvalue=[]
            newBamlabel=[]
            for Bamlabel in df.Bamlabel.unique():
                n=len(df[df.Bamlabel==Bamlabel])
                x=df.loc[df.Bamlabel==Bamlabel,'bin']
                y=df.loc[df.Bamlabel==Bamlabel,'value']
                xnew=np.linspace(x.min(),x.max(),5*n)
                func=make_interp_spline(x, y, k=3)
                ynew=func(xnew)
                newBamlabel+=[Bamlabel]*n*5
                newbin=newbin+list(xnew)
                newvalue=newvalue+list(ynew)
            #插值后赋值给一个新Dataframe
            df_new['Bamlabel']=newBamlabel
            df_new['bin']=newbin
            df_new['value']=newvalue
            #开始画图
            pal=sns.color_palette(Pal,len(df_new['Bamlabel'].unique()))
            g=sns.lineplot(x="bin",y='value',hue='Bamlabel',ax = axes,data=df_new,palette=pal)
            #保留左边和底部坐标轴
            axes.spines['right'].set_visible(False)
            axes.spines['top'].set_visible(False)
            axes.yaxis.set_ticks_position('left')
            axes.xaxis.set_ticks_position('bottom')
            #移除图例标题
            handles, labels = axes.get_legend_handles_labels()
            axes.legend(handles=handles[1:], labels=labels[1:])
            g.set(title=Bedlabel)
            g.set_xticks(xtick)
            g.set_xticklabels(xticklabel)
            g.set_xlabel("")
            g.set_ylabel(ylabel)
            g.figure.set_size_inches(4.5,3.5)
            #plt.tight_layout()
            fig.set_tight_layout(True)
        plt.savefig(outfigname)
                    
if Inputbamlist:

    if referencePoint:
        Ylength=int(up+down)
        ylength=int(Ylength/wdsize)
        left=int(1)
        medium=int(ylength/2)
        right=int(ylength)
        NumOfBin = int((int(up+down))/wdsize)
        if referencePoint=="start":
            out_profile=open(outf[0],"w")
            out_matrix=open(outf[1],"w")
            out_figure=open(figurename,"w")
            figurelist=[]
            #fullmat header
            header = ["Bamlabel","Bedlabel"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
            print(*header,sep="\t",end="\n",file=out_matrix)
            #average header
            print("Bamlabel","Bedlabel","bin","value",sep="\t",end="\n",file=out_profile)
            for Bam in range(len(bamlist)):
                bamf=pysam.AlignmentFile(bamlist[Bam],format[args.format])
                controlf=pysam.AlignmentFile(Inputbamlist[Bam],format[args.format])
                bamTotalReads=bamf.count()
                controlTotalReads=controlf.count()
                bamrCoverageDict=bamTodict(bamf)
                controlCoverageDict=bamTodict(controlf)
                for Bed in range(len(bedlist)):
                    values=bedCount(bedlist[Bed],bamrCoverageDict,bamTotalReads,wdsize,up,down)
                    Invalues=bedCount(bedlist[Bed],controlCoverageDict,controlTotalReads,wdsize,up,down)
                    #values=bamTocountStart(bedlist[Bed],bamlist[Bam],wdsize,up,down)
                    #Invalues=bamTocountStart(bedlist[Bed],Inputbamlist[Bam],wdsize,up,down)
                    tup1=np.array(values)	
                    Intup=np.array(Invalues)
                    tup=tup1-Intup
                    tup[tup<0]=0
                    #去掉极值
                    q=np.quantile(tup,0.998)
                    tup[tup>q]=q
                    tup2=tup.mean(axis=0)
                    figurelist.append(tup2)
                    Values=list(tup)
                    #print matrix
                    for index in range(len(values)):
                        print(Bamlabels[Bam],Bedlabels[Bed],*Values[index],sep="\t",end="\n",file=out_matrix)
                    #print average bin value
                    i=1
                    for binValue in tup2:
                        print(Bamlabels[Bam],Bedlabels[Bed],i,binValue,sep="\t",end="\n",file=out_profile)
                        i+=1
            out_matrix.close()
            out_profile.close()
            dataframe=pd.read_csv(outf[0],sep='\t')
            Xticks=[left,medium,right]
            Xticklabel=["-"+L+"kb","TSS",R+"kb"]
            Ylabel="Normalized cuts"
            plotfacet(Facet,dataframe,Xticks,Xticklabel,Ylabel,figurename,args.scaleY)
        elif referencePoint=="end":
            out_profile=open(outf[0],"w")
            out_matrix=open(outf[1],"w")
            out_figure=open(figurename,"w")
            figurelist=[]
            #fullmat header
            header = ["Bamlabel","Bedlabel"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
            print(*header,sep="\t",end="\n",file=out_matrix)
            #average header
            print("Bamlabel","Bedlabel","bin","value",sep="\t",end="\n",file=out_profile)
            for Bam in range(len(bamlist)):
                bamf=pysam.AlignmentFile(bamlist[Bam],format[args.format])
                controlf=pysam.AlignmentFile(Inputbamlist[Bam],format[args.format])
                bamTotalReads=bamf.count()
                controlTotalReads=controlf.count()
                bamrCoverageDict=bamTodict(bamf)
                controlCoverageDict=bamTodict(controlf)
                for Bed in range(len(bedlist)):
                    values=bedCountEnd(bedlist[Bed],bamrCoverageDict,bamTotalReads,wdsize,up,down)
                    Invalues=bedCountEnd(bedlist[Bed],controlCoverageDict,controlTotalReads,wdsize,up,down)
                    #values=bamTocountEnd(bedlist[Bed],bamlist[Bam],wdsize,up,down)
                    #Invalues=bamTocountEnd(bedlist[Bed],Inputbamlist[Bam],wdsize,up,down)
                    tup1=np.array(values)
                    Intup=np.array(Invalues)
                    tup=tup1-Intup
                    tup[tup<0]=0
                    #去掉极值
                    q=np.quantile(tup,0.998)
                    tup[tup>q]=q
                    tup2=tup.mean(axis=0)
                    Values=list(tup)
                    figurelist.append(tup2)
                    #print matrix
                    for index in range(len(values)):
                        print(Bamlabels[Bam],Bedlabels[Bed],*Values[index],sep="\t",end="\n",file=out_matrix)
                    #print average bin value
                    i=1
                    for binValue in tup2:
                        print(Bamlabels[Bam],Bedlabels[Bed],i,binValue,sep="\t",end="\n",file=out_profile)
                        i+=1
            out_matrix.close()
            out_profile.close()
            dataframe=pd.read_csv(outf[0],sep='\t')
            Xticks=[left,medium,right]
            Xticklabel=["-"+L+"kb","TTS",R+"kb"]
            Ylabel="Normalized cuts"
            plotfacet(Facet,dataframe,Xticks,Xticklabel,Ylabel,figurename,args.scaleY)
        else:
            pass

    elif scaleRegion:
        left=int(1)
        mediumA=int(up/wdsize)
        mediumB=mediumA+blocknum
        right=int(down/wdsize)+mediumA+blocknum
        NumOfBin = mediumA+blocknum+int(down/wdsize)
        out_profile=open(outf[0],"w")
        out_matrix=open(outf[1],"w")
        out_figure=open(figurename,"w")
        figurelist=[]
        #fullmat header
        header = ["Bamlabel","Bedlabel"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
        print(*header,sep="\t",end="\n",file=out_matrix)
        #average header
        print("Bamlabel","Bedlabel","bin","value",sep="\t",end="\n",file=out_profile)
        for Bam in range(len(bamlist)):
            bamf=pysam.AlignmentFile(bamlist[Bam],format[args.format])
            controlf=pysam.AlignmentFile(Inputbamlist[Bam],format[args.format])
            bamTotalReads=bamf.count()
            controlTotalReads=controlf.count()
            bamrCoverageDict=bamTodict(bamf)
            controlCoverageDict=bamTodict(controlf)
            for Bed in range(len(bedlist)):
                values=BodyCount(bedlist[Bed],bamrCoverageDict,bamTotalReads,wdsize,up,down,blocknum)
                Invalues=BodyCount(bedlist[Bed],controlCoverageDict,controlTotalReads,wdsize,up,down,blocknum)
#                 values=bamtobodycount(bedlist[Bed],bamlist[Bam],wdsize,up,down,blocknum)
#                 Invalues=bamtobodycount(bedlist[Bed],Inputbamlist[Bam],wdsize,up,down,blocknum)
                Intup=np.array(Invalues)
                tup1=np.array(values)
                tup=tup1-Intup
                tup[tup<0]=0
                #去掉极值
                q=np.quantile(tup,0.998)
                tup[tup>q]=q
                tup2=tup.mean(axis=0)
                Values=list(tup)
                figurelist.append(tup2)
                #print matrix
                for index in range(len(values)):
                    print(Bamlabels[Bam],Bedlabels[Bed],*values[index],sep="\t",end="\n",file=out_matrix)
                #print average bin value
                i=1
                for binValue in tup2:
                    print(Bamlabels[Bam],Bedlabels[Bed],i,binValue,sep="\t",end="\n",file=out_profile)
                    i+=1
        out_matrix.close()
        out_profile.close()
        dataframe=pd.read_csv(outf[0],sep='\t')
        Xticks=[left,mediumA,mediumB,right]
        Xticklabel=["-"+L+"kb","TSS","TTS",R+"kb"]
        Ylabel="Normalized cuts"
        plotfacet(Facet,dataframe,Xticks,Xticklabel,Ylabel,figurename,args.scaleY)
    else:
        pass
    
else:

    if referencePoint:
        ylength=int(int(up+down)/wdsize)
        left=int(1)
        medium=int(ylength/2)+0.5
        right=int(ylength)
        NumOfBin = int((int(up+down))/wdsize)
        if referencePoint=="start":
            out_profile=open(outf[0],"w")
            out_matrix=open(outf[1],"w")
            out_figure=open(figurename,"w")
            figurelist=[]
            #fullmat header
            header = ["Bamlabel","Bedlabel"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
            print(*header,sep="\t",end="\n",file=out_matrix)
            #average header
            print("Bamlabel","Bedlabel","bin","value",sep="\t",end="\n",file=out_profile)
            for Bam in range(len(bamlist)):
                bamf=pysam.AlignmentFile(bamlist[Bam],format[args.format])
                bamTotalReads=bamf.count()
                bamrCoverageDict=bamTodict(bamf)
                for Bed in range(len(bedlist)):
                    values=bedCount(bedlist[Bed],bamrCoverageDict,bamTotalReads,wdsize,up,down)
                    #values=bamTocountStart(bedlist[Bed],bamlist[Bam],wdsize,up,down)
                    tup1=np.array(values)
                    #去掉极值
                    q=np.quantile(tup1,0.998)
                    tup1[tup1>q]=q
                    tup2=tup1.mean(axis=0)
                    figurelist.append(tup2)
                    #print matrix
                    for index in range(len(values)):
                        print(Bamlabels[Bam],Bedlabels[Bed],*values[index],sep="\t",end="\n",file=out_matrix)
                    #print average bin value
                    i=1
                    for binValue in tup2:
                        print(Bamlabels[Bam],Bedlabels[Bed],i,binValue,sep="\t",end="\n",file=out_profile)
                        i+=1
            out_matrix.close()
            out_profile.close()
            dataframe=pd.read_csv(outf[0],sep='\t')
            Xticks=[left,medium,right]
            Xticklabel=["-"+L+"kb","TSS",R+"kb"]
            Ylabel="Normalized cuts"
            plotfacet(Facet,dataframe,Xticks,Xticklabel,Ylabel,figurename,args.scaleY)
        elif referencePoint=="end":
            out_profile=open(outf[0],"w")
            out_matrix=open(outf[1],"w")
            out_figure=open(figurename,"w")
            figurelist=[]
            #fullmat header
            header = ["Bamlabel","Bedlabel"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
            print(*header,sep="\t",end="\n",file=out_matrix)
            #average header
            print("Bamlabel","Bedlabel","bin","value",sep="\t",end="\n",file=out_profile)
            for Bam in range(len(bamlist)):
                bamf=pysam.AlignmentFile(bamlist[Bam],format[args.format])
                bamTotalReads=bamf.count()
                bamrCoverageDict=bamTodict(bamf)
                for Bed in range(len(bedlist)):
                    values=bedCountEnd(bedlist[Bed],bamrCoverageDict,bamTotalReads,wdsize,up,down)
                    #values=bamTocountEnd(bedlist[Bed],bamlist[Bam],wdsize,up,down)
                    tup1=np.array(values)
                    #去掉极值
                    q=np.quantile(tup1,0.998)
                    tup1[tup1>q]=q
                    tup2=tup1.mean(axis=0)
                    figurelist.append(tup2)
                    #print matrix
                    for index in range(len(values)):
                        print(Bamlabels[Bam],Bedlabels[Bed],*values[index],sep="\t",end="\n",file=out_matrix)
                    #print average bin value
                    i=1
                    for binValue in tup2:
                        print(Bamlabels[Bam],Bedlabels[Bed],i,binValue,sep="\t",end="\n",file=out_profile)
                        i+=1
            out_matrix.close()
            out_profile.close()
            dataframe=pd.read_csv(outf[0],sep='\t')
            Xticks=[left,medium,right]
            Xticklabel=["-"+L+"kb","TTS",R+"kb"]
            Ylabel="Normalized cuts"
            plotfacet(Facet,dataframe,Xticks,Xticklabel,Ylabel,figurename,args.scaleY)
        elif referencePoint=="center":
            out_profile=open(outf[0],"w")
            out_matrix=open(outf[1],"w")
            out_figure=open(figurename,"w")
            #fullmat header
            header = ["Bamlabel","Bedlabel"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
            print(*header,sep="\t",end="\n",file=out_matrix)
            #average header
            print("Bamlabel","Bedlabel","bin","value",sep="\t",end="\n",file=out_profile)
            for Bam in range(len(bamlist)):
                bamf=pysam.AlignmentFile(bamlist[Bam],format[args.format])
                bamTotalReads=bamf.count()
                bamrCoverageDict=bamTodict(bamf)
                for Bed in range(len(bedlist)):
                    values=bedCountCenter(bedlist[Bed],bamrCoverageDict,bamTotalReads,wdsize,up,down)
                    tup1=np.array(values)
                    #去掉极值
                    q=np.quantile(tup1,0.998)
                    tup1[tup1>q]=q
                    tup2=tup1.mean(axis=0)
                    #print matrix
                    for index in range(len(values)):
                        print(Bamlabels[Bam],Bedlabels[Bed],*values[index],sep="\t",end="\n",file=out_matrix)
                    #print average bin value
                    i=1
                    for binValue in tup2:
                        print(Bamlabels[Bam],Bedlabels[Bed],i,binValue,sep="\t",end="\n",file=out_profile)
                        i+=1
            out_matrix.close()
            out_profile.close()
            dataframe=pd.read_csv(outf[0],sep='\t')
            Xticks=[left,medium,right]
            Xticklabel=["-"+L+"kb","center",R+"kb"]
            Ylabel="Normalized read counts"
            plotfacet(Facet,dataframe,Xticks,Xticklabel,Ylabel,figurename,args.scaleY)
        else:
            pass
#print("Usage: python tmp1.py <<--scale-regions or --reference-point>> -b <bamfile> -G <bedfile> -g <genelist> -r <referencePoint {start or end}> -u <upstream_distance> -d <downstream_distance> -w <window_size> -bn <Equal_block_number> -o <outputfile>")

    elif scaleRegion:
        left=int(1)
        mediumA=int(up/wdsize)
        mediumB=mediumA+blocknum
        right=int(down/wdsize)+mediumA+blocknum
        NumOfBin = mediumA+blocknum+int(down/wdsize)
        out_profile=open(outf[0],"w")
        out_matrix=open(outf[1],"w")
        out_figure=open(figurename,"w")
        figurelist=[]
        #fullmat header
        header = ["Bamlabel","Bedlabel"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
        print(*header,sep="\t",end="\n",file=out_matrix)
        #average header
        print("Bamlabel","Bedlabel","bin","value",sep="\t",end="\n",file=out_profile)
        for Bam in range(len(bamlist)):
            bamf=pysam.AlignmentFile(bamlist[Bam],format[args.format])
            bamTotalReads=bamf.count()
            bamrCoverageDict=bamTodict(bamf)
            for Bed in range(len(bedlist)):
                values=BodyCount(bedlist[Bed],bamrCoverageDict,bamTotalReads,wdsize,up,down,blocknum)
                #values=bamtobodycount(bedlist[Bed],bamlist[Bam],wdsize,up,down,blocknum)
                tup1=np.array(values)
                #去掉极值
                q=np.quantile(tup1,0.998)
                tup1[tup1>q]=q
                tup2=tup1.mean(axis=0)
                figurelist.append(tup2)
                #print matrix
                for index in range(len(values)):
                    print(Bamlabels[Bam],Bedlabels[Bed],*values[index],sep="\t",end="\n",file=out_matrix)
                #print average bin value
                i=1
                for binValue in tup2:
                    print(Bamlabels[Bam],Bedlabels[Bed],i,binValue,sep="\t",end="\n",file=out_profile)
                    i+=1
        out_matrix.close()
        out_profile.close()
        dataframe=pd.read_csv(outf[0],sep='\t')
        Xticks=[left,mediumA,mediumB,right]
        Xticklabel=["-"+L+"kb","TSS","TTS",R+"kb"]
        Ylabel="Normalized cuts"
        plotfacet(Facet,dataframe,Xticks,Xticklabel,Ylabel,figurename,args.scaleY)
    else:
        pass
