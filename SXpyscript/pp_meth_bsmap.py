# -*- coding: UTF-8 -*-
from __future__ import print_function
from __future__ import division
#from numba import jit
#import sys, getopt
import argparse
import pandas as pd
import numpy as np
import gzip
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
 python **.py --referencePoint start/end -b bsmap.methratio.gz,... -G **.bed,**.bed... -L bamlabel1,bamlabel2 -l bedlabel1,bedlabel2 --Facet Bam/Bed -u 1000 -d 1000 -w 50 -o prefix
 python **.py --scaleRegion -b bsmap.methratio.gz,... -G **.bed,**.bed... -L bamlabel1,bamlabel2 -l bedlabel1,bedlabel2 --Facet Bam/Bed -u 1000 -d 1000 -w 50 -N 100 -o prefix
'''
parser = argparse.ArgumentParser(description="This tool creates methylation level profile over sets of genomic regions. Typically, these regions are genes, but any other regions defined in BED will work.", 
                                 formatter_class= argparse.RawTextHelpFormatter,
                                 prog='base_maker',
                                 epilog=example_text)
group=parser.add_mutually_exclusive_group(required= True)
group.add_argument('--referencePoint',type=str,help="<referencePoint {start or end}>"+"\t"+"plot the distribution of the reads at reference point from the genome.")
group.add_argument('--scaleRegion',action='store_true',help="<Scaled Region>"+"\t"+"plot the distribution of the reads over the region from the genome.")

parser.add_argument('--bam','-b',type= str,help="<bam or sam files>"+"\t"+"data files(Support for multiple files, Separated by commas), should be sorted format (Note: Need to filter out low-quality reads first).")
#parser.add_argument('--format',type= str,help="<gz/text>"+"\t"+"Specify Alignment file format. (Default=gz)",default="gz")
parser.add_argument('--minCov','-c',type= int,help="<minimal C coverage>"+"\t"+"The minimal Cytosine coverage per base to be count. (Default=5)",default=5)

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

parser.add_argument('--Facet',type=str,help="<Bam/Bed>"+"\t"+"Note: If you only have one bam(CGmap), all 3 contexts of methylation will be drawn into ONE plot. If facet by bam file, each bam file as a subplot to draw all bedfiles. If facet by bed file, each bed site as a subplot to draw all bamfiles reads. Each methylation context will save as single pdf respectively.",required= True)
parser.add_argument('--scaleY',action='store_true',help="<Optional argument>"+"\t"+"Whether or not to scale Y limits to the same. (Default=False)",default=False)

parser.add_argument('--outfileprefix','-o',type=str,help="<OUTPUT file name Prefix>"+"\t"+"It will output full matrix,average data and average profile.",required= True)
#parser.add_argument('--figure','-f',type=str,help="<Figure name>"+"\t"+"Name of the plot name. (plot file format is PDF).",required= True)
args = parser.parse_args()


bamlist=args.bam.split(",")
minCov=args.minCov

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

#figurename=args.figure
outf=args.outfileprefix
figurename=outf

format={'gz':'rb','text':'r'}

def CGmapTodict(Input,minimalCoverage):
    CG={}
    CHG={}
    CHH={}
    #跳过表头
    next(Input)
    for line in Input:
        col=line.rstrip().split()
        chr=str(col[0])
        site=str(col[1])
        context=str(col[3])
        allC=int(float(col[5]))
        methyC=int(float(col[6]))
        if allC >= minimalCoverage:
            if context == 'CG':
                if chr in CG:
                    CG[chr][site]=(methyC,allC)
                else:
                    CG[chr]={}
                    CG[chr][site]=(methyC,allC)
            elif context == 'CHG':
                if chr in CHG:
                    CHG[chr][site]=(methyC,allC)
                else:
                    CHG[chr]={}
                    CHG[chr][site]=(methyC,allC)
            elif context == 'CHH':
                if chr in CHH:
                    CHH[chr][site]=(methyC,allC)
                else:
                    CHH[chr]={}
                    CHH[chr][site]=(methyC,allC)
            else:
                pass
    region={'CG':CG,'CHG':CHG,'CHH':CHH}
    return region

def bedMethyStart(bedfile,readdict,windows,up,down):
    metC=0
    sumC=0
    methyList=[]
    bedf=open(bedfile,'r')
    for line in bedf:
        line=line.strip().split()
        sublist=[]
        if line[5] is "+":
            left=int(line[1])-up
            right=int(line[1])+down
            if line[0] in readdict:
                for k in range(left,right,windows):
                    for j in range(k,(k+windows)):
                        if str(j) in readdict[line[0]]:
                            sumC+=readdict[line[0]][str(j)][1]
                            metC+=readdict[line[0]][str(j)][0]
                        else:
                            pass
                    if sumC == 0:
                        methyratio=np.nan
                    else:
                        methyratio=metC/sumC
                    sublist.append(methyratio)
                    metC=0
                    sumC=0
            else:
                for k in range(left,right,windows):
                    sublist.append(np.nan)
            methyList.append(sublist)
        else:
            left=int(line[2])-down
            right=int(line[2])+up
            if line[0] in readdict:
                for k in range(right,left,-windows):
                    for j in range(k,(k-windows),-1):
                        if str(j) in readdict[line[0]]:
                            sumC+=readdict[line[0]][str(j)][1]
                            metC+=readdict[line[0]][str(j)][0]
                        else:
                            pass
                    if sumC == 0:
                        methyratio=np.nan
                    else:
                        methyratio=metC/sumC
                    sublist.append(methyratio)
                    metC=0
                    sumC=0
            else:
                for k in range(right,left,-windows):
                    sublist.append(np.nan)
            methyList.append(sublist)
    return methyList
    bedfile.close()

def bedMethyEnd(bedfile,readdict,windows,up,down):
    metC=0
    sumC=0
    methyList=[]
    bedf=open(bedfile,'r')
    for line in bedf:
        line=line.strip().split()
        sublist=[]
        if line[5] is "+":
            left=int(line[2])-up
            right=int(line[2])+down
            if line[0] in readdict:
                for k in range(left,right,windows):
                    for j in range(k,(k+windows)):
                        if str(j) in readdict[line[0]]:
                            sumC+=readdict[line[0]][str(j)][1]
                            metC+=readdict[line[0]][str(j)][0]
                        else:
                            pass
                    if sumC == 0:
                        methyratio=np.nan
                    else:
                        methyratio=metC/sumC
                    sublist.append(methyratio)
                    metC=0
                    sumC=0
            else:
                for k in range(left,right,windows):
                    sublist.append(np.nan)
            methyList.append(sublist)
        else:
            left=int(line[1])-down
            right=int(line[1])+up
            if line[0] in readdict:
                for k in range(right,left,-windows):
                    for j in range(k,(k-windows),-1):
                        if str(j) in readdict[line[0]]:
                            sumC+=readdict[line[0]][str(j)][1]
                            metC+=readdict[line[0]][str(j)][0]
                        else:
                            pass
                    if sumC == 0:
                        methyratio=np.nan
                    else:
                        methyratio=metC/sumC
                    sublist.append(methyratio)
                    metC=0
                    sumC=0
            else:
                for k in range(right,left,-windows):
                    sublist.append(np.nan)
            methyList.append(sublist)
    return methyList
    bedfile.close()


def bedMethyBody(bedfile,readdict,windows,up,down,Blocknum):
    metC=0
    sumC=0
    methyList=[]
    bedf=open(bedfile,'r')
    for line in bedf:
        line=line.strip().split()
        sublist=[]
        reside=Blocknum-int(((int(line[2])-int(line[1]))%Blocknum))
        blockSize=int(((int(line[2])-int(line[1]))+reside)/Blocknum)
        if line[5] is "+":
            if line[0] in readdict:
                for k in range(int(line[1])-up,int(line[1]),windows):
                    for j in range(k,(k+windows)):
                        if str(j) in readdict[line[0]]:
                            sumC+=readdict[line[0]][str(j)][1]
                            metC+=readdict[line[0]][str(j)][0]
                        else:
                            pass
                    if sumC == 0:
                        methyratio=np.nan
                    else:
                        methyratio=metC/sumC
                    sublist.append(methyratio)
                    metC=0
                    sumC=0
            else:
                for k in range(int(line[1])-up,int(line[1]),windows):
                    sublist.append(np.nan)
            if line[0] in readdict:
                for x in range(int(line[1]),int(line[2])+reside,blockSize):
                    for y in range(x,(x+blockSize)):
                        if str(y) in readdict[line[0]]:
                            sumC+=readdict[line[0]][str(y)][1]
                            metC+=readdict[line[0]][str(y)][0]
                        else:
                            pass
                    if sumC == 0:
                        methyratio=np.nan
                    else:
                        methyratio=metC/sumC
                    sublist.append(methyratio)
                    metC=0
                    sumC=0
            else:
                for x in range(int(line[1]),int(line[2])+reside,blockSize):
                    sublist.append(np.nan)
            if line[0] in readdict:
                for u in range(int(line[2]),int(line[2])+down,windows):
                    for v in range(u,(u+windows)):
                        if str(v) in readdict[line[0]]:
                            sumC+=readdict[line[0]][str(v)][1]
                            metC+=readdict[line[0]][str(v)][0]
                        else:
                            pass
                    if sumC == 0:
                        methyratio=np.nan
                    else:
                        methyratio=metC/sumC
                    sublist.append(methyratio)
                    metC=0
                    sumC=0
            else:
                for u in range(int(line[2]),int(line[2])+down,windows):
                    sublist.append(np.nan)
            methyList.append(sublist)
        else:
            if line[0] in readdict:
                for k in range(int(line[2])+up,int(line[2]),-windows):
                    for j in range(k,(k-windows),-1):
                        if str(j) in readdict[line[0]]:
                            sumC+=readdict[line[0]][str(j)][1]
                            metC+=readdict[line[0]][str(j)][0]
                        else:
                            pass
                    if sumC == 0:
                        methyratio=np.nan
                    else:
                        methyratio=metC/sumC
                    sublist.append(methyratio)
                    metC=0
                    sumC=0
            else:
                for k in range(int(line[2])+up,int(line[2]),-windows):
                    sublist.append(np.nan)
            if line[0] in readdict:
                for x in range(int(line[2]),int(line[1])-reside,-blockSize):
                    for y in range(x,(x-blockSize),-1):
                        if str(y) in readdict[line[0]]:
                            sumC+=readdict[line[0]][str(y)][1]
                            metC+=readdict[line[0]][str(y)][0]
                        else:
                            pass
                    if sumC == 0:
                        methyratio=np.nan
                    else:
                        methyratio=metC/sumC
                    sublist.append(methyratio)
                    metC=0
                    sumC=0
            else:
                for x in range(int(line[2]),int(line[1])-reside,-blockSize):
                    sublist.append(np.nan)
            if line[0] in readdict:
                for u in range(int(line[1]),int(line[1])-down,-windows):
                    for v in range(u,(u-windows),-1):
                        if str(v) in readdict[line[0]]:
                            sumC+=readdict[line[0]][str(v)][1]
                            metC+=readdict[line[0]][str(v)][0]
                        else:
                            pass
                    if sumC == 0:
                        methyratio=np.nan
                    else:
                        methyratio=metC/sumC
                    sublist.append(methyratio)
                    metC=0
                    sumC=0
            else:
                for u in range(int(line[1]),int(line[1])-down,-windows):
                    sublist.append(np.nan)
            methyList.append(sublist)
    return methyList
    bedfile.close()

def plotsmooth(Values):
    n=len(Values)
    num=range(0,n)
    x=np.array(num)
    xnew=np.linspace(x.min(),x.max(),5*n)
    func=make_interp_spline(x, Values, k=3)
    ynew=func(xnew)
    plt.plot(xnew,ynew)

def dataframe_Interpolation(df,byWhichlabel):
    df_new=pd.DataFrame()
    newbin=[]
    newvalue=[]
    newlabel=[]
    #遍历每条线
    for label in df[byWhichlabel].unique():
        n=len(df[df[byWhichlabel]==label])
        x=df.loc[df[byWhichlabel]==label,'bin']
        y=df.loc[df[byWhichlabel]==label,'value']
        xnew=np.linspace(x.min(),x.max(),5*n)
        func=make_interp_spline(x, y, k=3)
        ynew=func(xnew)
        newlabel+=[label]*n*5
        newbin=newbin+list(xnew)
        newvalue=newvalue+list(ynew)
    #插值后赋值给一个新Dataframe
    df_new[byWhichlabel]=newlabel
    df_new['bin']=newbin
    df_new['value']=newvalue
    return df_new

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
            Context=df.Context.unique()
            #三种context分开画，各输出一张图
            for context in Context:
                contextdf=df[df.Context==context]
                fig, axes = plt.subplots(nrows=nrow,ncols=ncol)
                axes=axes.flatten()
                if scaleY:
                    max_value=max(contextdf.value)
                    min_value=min(contextdf.value)
                    max_y=max_value+(max_value-min_value)*0.05
                    min_y=max(0,min_value-(max_value-min_value)*0.05)
                    custom_ylim = (min_y,max_y)
                    # Setting the values for all axes.
                    plt.setp(axes,ylim=custom_ylim)
                #遍历每个分面
                for Bamlabel, ax in zip(u, axes):
                    subdf=contextdf[contextdf.Bamlabel==Bamlabel]
                    subdf_new=dataframe_Interpolation(subdf,'Bedlabel')
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
                plt.savefig(outfigname+"_"+context+'.pdf')
        else:
            Context=df.Context.unique()
            fig, axes = plt.subplots(ncols=len(Context))
            if scaleY:
                max_value=max(df.value)
                min_value=min(df.value)
                max_y=max_value+(max_value-min_value)*0.05
                min_y=max(0,min_value-(max_value-min_value)*0.05)
                custom_ylim = (min_y,max_y)
                # Setting the values for all axes.
                plt.setp(axes,ylim=custom_ylim)
            for context, ax in zip(Context, axes):
                subdf=df[df.Context==context]
                subdf_new=dataframe_Interpolation(subdf,'Bedlabel')
                #开始画图
                pal=sns.color_palette(Pal,len(subdf_new['Bedlabel'].unique()))
                g=sns.lineplot(x="bin",y='value',hue='Bedlabel',ax = ax,data=subdf_new,palette=pal)
                #保留左边和底部坐标轴
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.yaxis.set_ticks_position('left')
                ax.xaxis.set_ticks_position('bottom')
                #移除图例标题
                handles, labels = ax.get_legend_handles_labels()
                ax.legend(handles=handles[1:], labels=labels[1:])
                #设置标签
                g.set(title=context)
                g.set_xticks(xtick)
                g.set_xticklabels(xticklabel)
                g.set_xlabel("")
                g.set_ylabel(ylabel)
                g.figure.set_size_inches(len(Context)*4.5,3.5)
                #plt.tight_layout()
            fig.set_tight_layout(True)
            plt.savefig(outfigname+'.pdf')
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
            Context=df.Context.unique()
            #三种context分开画，各输出一张图
            for context in Context:
                contextdf=df[df.Context==context]
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
                    subdf=contextdf[contextdf.Bedlabel==Bedlabel]
                    subdf_new=dataframe_Interpolation(subdf,'Bamlabel')
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
                plt.savefig(outfigname+"_"+context+'.pdf')
        else:
            Context=df.Context.unique()
            fig, axes = plt.subplots(ncols=len(Context))
            if scaleY:
                max_value=max(df.value)
                min_value=min(df.value)
                max_y=max_value+(max_value-min_value)*0.05
                min_y=max(0,min_value-(max_value-min_value)*0.05)
                custom_ylim = (min_y,max_y)
                # Setting the values for all axes.
                plt.setp(axes,ylim=custom_ylim)
            for context, ax in zip(Context, axes):
                subdf=df[df.Context==context]
                subdf_new=dataframe_Interpolation(subdf,'Bamlabel')
                #开始画图
                pal=sns.color_palette(Pal,len(subdf_new['Bamlabel'].unique()))
                g=sns.lineplot(x="bin",y='value',hue='Bamlabel',ax = ax,data=subdf_new,palette=pal)
                #保留左边和底部坐标轴
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.yaxis.set_ticks_position('left')
                ax.xaxis.set_ticks_position('bottom')
                #移除图例标题
                handles, labels = ax.get_legend_handles_labels()
                ax.legend(handles=handles[1:], labels=labels[1:])
                #设置标签
                g.set(title=context)
                g.set_xticks(xtick)
                g.set_xticklabels(xticklabel)
                g.set_xlabel("")
                g.set_ylabel(ylabel)
                g.figure.set_size_inches(len(Context)*4.5,3.5)
                #plt.tight_layout()
            fig.set_tight_layout(True)
            plt.savefig(outfigname+'.pdf')

if referencePoint:
    ylength=int(int(up+down)/wdsize)
    left=int(1)
    medium=int(ylength/2)+0.5
    right=int(ylength)
    NumOfBin = int((int(up+down))/wdsize)
    if referencePoint=="start":
        out_profile=open(outf+'.mean',"w")
        out_matrix=open(outf+'.mat',"w")
        #out_figure=open(figurename,"w")
        #fullmat header
        header = ["Bamlabel","Bedlabel","Context"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
        print(*header,sep="\t",end="\n",file=out_matrix)
        #average header
        print("Bamlabel","Bedlabel","Context","bin","value",sep="\t",end="\n",file=out_profile)
        for Bam in range(len(bamlist)):
            bamf=gzip.open(bamlist[Bam],'rb')
            MethyDict = CGmapTodict(bamf,minCov)
            for Bed in range(len(bedlist)):
                for context in ['CG','CHG','CHH']:
                    values=bedMethyStart(bedlist[Bed],MethyDict[context],wdsize,up,down)
                    tup1=np.array(values)
                    tup2=np.nanmean(tup1,axis=0)
                    #print matrix
                    for index in range(len(values)):
                        print(Bamlabels[Bam],Bedlabels[Bed],context,*values[index],sep="\t",end="\n",file=out_matrix)
                    #print average bin value
                    i=1
                    for binValue in tup2:
                        print(Bamlabels[Bam],Bedlabels[Bed],context,i,binValue,sep="\t",end="\n",file=out_profile)
                        i+=1
        out_matrix.close()
        out_profile.close()
        dataframe=pd.read_csv(outf+'.mean',sep='\t')
        Xticks=[left,medium,right]
        Xticklabel=["-"+L+"kb","TSS",R+"kb"]
        Ylabel="Methylation level"
        plotfacet(Facet,dataframe,Xticks,Xticklabel,Ylabel,figurename,args.scaleY)
    elif referencePoint=="end":
        out_profile=open(outf+'.mean',"w")
        out_matrix=open(outf+'.mat',"w")
        #out_figure=open(figurename,"w")
        #fullmat header
        header = ["Bamlabel","Bedlabel","Context"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
        print(*header,sep="\t",end="\n",file=out_matrix)
        #average header
        print("Bamlabel","Bedlabel","Context","bin","value",sep="\t",end="\n",file=out_profile)
        for Bam in range(len(bamlist)):
            bamf=gzip.open(bamlist[Bam],'rb')
            MethyDict = CGmapTodict(bamf,minCov)
            for Bed in range(len(bedlist)):
                for context in ['CG','CHG','CHH']:
                    values=bedMethyEnd(bedlist[Bed],MethyDict[context],wdsize,up,down)
                    tup1=np.array(values)
                    tup2=np.nanmean(tup1,axis=0)
                    #print matrix
                    for index in range(len(values)):
                        print(Bamlabels[Bam],Bedlabels[Bed],context,*values[index],sep="\t",end="\n",file=out_matrix)
                    #print average bin value
                    i=1
                    for binValue in tup2:
                        print(Bamlabels[Bam],Bedlabels[Bed],context,i,binValue,sep="\t",end="\n",file=out_profile)
                        i+=1
        out_matrix.close()
        out_profile.close()
        dataframe=pd.read_csv(outf+'.mean',sep='\t')
        Xticks=[left,medium,right]
        Xticklabel=["-"+L+"kb","TTS",R+"kb"]
        Ylabel="Methylation level"
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
    out_profile=open(outf+'.mean',"w")
    out_matrix=open(outf+'.mat',"w")
    #out_figure=open(figurename,"w")
    #fullmat header
    header = ["Bamlabel","Bedlabel","Context"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
    print(*header,sep="\t",end="\n",file=out_matrix)
    #average header
    print("Bamlabel","Bedlabel","Context","bin","value",sep="\t",end="\n",file=out_profile)
    for Bam in range(len(bamlist)):
        bamf=gzip.open(bamlist[Bam],'rb')
        MethyDict = CGmapTodict(bamf,minCov)
        for Bed in range(len(bedlist)):
            for context in ['CG','CHG','CHH']:
                values=bedMethyBody(bedlist[Bed],MethyDict[context],wdsize,up,down,blocknum)
                tup1=np.array(values)
                tup2=np.nanmean(tup1,axis=0)
                #print matrix
                for index in range(len(values)):
                    print(Bamlabels[Bam],Bedlabels[Bed],context,*values[index],sep="\t",end="\n",file=out_matrix)
                #print average bin value
                i=1
                for binValue in tup2:
                    print(Bamlabels[Bam],Bedlabels[Bed],context,i,binValue,sep="\t",end="\n",file=out_profile)
                    i+=1
    out_matrix.close()
    out_profile.close()
    dataframe=pd.read_csv(outf+'.mean',sep='\t')
    Xticks=[left,mediumA,mediumB,right]
    Xticklabel=["-"+L+"kb","TSS","TTS",R+"kb"]
    Ylabel="Methylation level"
    plotfacet(Facet,dataframe,Xticks,Xticklabel,Ylabel,figurename,args.scaleY)
