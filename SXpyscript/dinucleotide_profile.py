# -*- coding: UTF-8 -*-
from __future__ import print_function
from __future__ import division
from pyfasta import Fasta
import sys, getopt
import gzip
import numpy as np
import pandas as pd
from math import sqrt,ceil
import matplotlib.pyplot as plt
import seaborn as sns
#from scipy.interpolate import interp1d
from scipy.interpolate import make_interp_spline, BSpline
plt.switch_backend('PDF')
plt.rcParams['pdf.fonttype'] = 42
##@jit

def message():
    print("\n"+"Usage: python **.py --reference-point or --scale-regions [options]"+"\n")
    print("Description:"+"\n"+"\t"+"\t"+"This tool creates Dinucleotide density profile over sets of genomic regions."+"\n")
    print("Options:"+"\n")
    print("\t"+"-v, "+"--version"+"\t"+"show program's version number and exit")
    print("\t"+"-h, "+"--help"+"\t"+"show this help message and exit.")
    print("\t"+"--reference-point or --scale-regions"+"\t"+"Reference-point refers to a position within a BED region(e.g., the starting point). In this mode, only those genomicpositions before (upstream) and/or after (downstream) of the reference point will be plotted."+"In the scale-regions mode, all regions in the BED file are stretched or shrunken to the length (in bases) indicated by the user.")
    print("\t"+"-b"+"    <Reference>"+"\t"+"Reference genome fasta files")
    print("\t"+"-G"+"    <bedfiles>"+"\t"+"data files(Support for multiple files, Separated by commas), should be original bed files.")
    #print("\t"+"-p"+"    <pattern {Any Dinucleotide, for example [CG] or [AT] or [WW] or [SS] and so on.}>")
    #print("\t"+"-s"+"    <strand {forward or reverse}>"+"\t"+"which DNA strand.")
    print("\t"+"-r"+"    <referencePoint {start or end}>"+"\t"+"The reference point for the plotting could be either the region start (TSS), the region end (TES) of the region.")
    print("\t"+"-u"+"    <upstream distance>"+"\t"+"Distance upstream of the reference-point selected.")
    print("\t"+"-d"+"    <downstream distance>"+"\t"+"Distance downstream of the reference-point selected.")
    print("\t"+"-w"+"    <window size>"+"\t"+"Length, in bases, binSize for averaging the score over the regions length.")
    print("\t"+"-N"+"    <Equally fragment number>"+"\t"+"Equally divided into the region into several fragments.")
    print("\t"+"-l"+"    <legend names>"+"\t"+"Bed file label")
    print("\t"+"-o"+"    <outputfiles average data,final matrix>"+"\t"+"File name to save the average data file and final matrix file. (Separated by commas)")
    print("\t"+"-F"+"    <Figure name>"+"\t"+"Name of the plot name. (plot file format is PDF)"+"\n")
    print("An example usage is:"+"\n")
    print("\t"+"python **.py --reference-point -b **.fasta -G **.bed,**.bed... -g **.txt,**.txt... -r start/end -u 1000 -d 1000 -w 50 -l **,** -o **_dens,**_matrix -F **.pdf"+"\n")
    print("\t"+"python **.py --scale-regions -b **.fasta -G **.bed,**.bed... -g **.txt,**.txt... -u 1000 -d 1000 -w 50 -N 100 -l **,** -o **_dens,**_matrix -F **.pdf"+"\n")
    return

opts,args = getopt.getopt(sys.argv[1:],'-h-v-b:-G:-s:-p:-u:-r:-d:-w:-N:-o:-l:-F:',['help','version','scale-regions','reference-point'])
for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        message()
    if opt_name in ('-v','--version'):
        print("**.py 1.0 2020-05-28")
    if opt_name in ('-b'):
        fasta = opt_value
    if opt_name in ('-G'):
        bedfile = opt_value
        bedlist=bedfile.split(",")
    if opt_name in ('-l'):
        legend = opt_value
        bedlabel=legend.split(",")
    if opt_name in ('-s'):
        strand = opt_value
#     if opt_name in ('-p'):
#         Pattern = opt_value
    if opt_name in ('-u'):
        upstream = opt_value
        up=int(upstream)
        L=str(up/1000)
    if opt_name in ('-d'):
        downstream = opt_value
        down=int(downstream)
        R=str(down/1000)
    if opt_name in ('-r'):
        referencePoint = opt_value
    if opt_name in ('-w'):
        window = opt_value
        windows=int(window)
    if opt_name in ('-N'):
        BN = opt_value
        BlockNum=int(BN)
    if opt_name in ('-F'):
        figurename = opt_value
    if opt_name in ('-o'):
        outputfile = opt_value
        outf=outputfile.split(",")
    if opt_name not in ('-h,--help,-v,--version,-b,-G,-r,-u,-d,-w,-N,-o,-l,-F,-S,-c,--scale-regions,--reference-point'):
        message()
    #print("Usage: python tmp1.py <--scale-regions or --reference-point> -b <CGmapfiles> -G <bedfiles> -r <referencePoint {start or end}> -u <upstream_distance> -d <downstream_distance> -N <Equal_block_number> -l <legend_names> -o <outputfiles average.matrix> -F <Figure_name>")
Option={}
for option in opts:
    Option[option[0]]=option[1]

Dinuc=['AA','AT','AC','AG',
       'TA','TT','TC','TG',
       'CA','CT','CC','CG',
       'GA','GT','GG','GC']

Pattern=[['AA'],['AT'],['AC'],['AG'],
         ['TA'],['TT'],['TC'],['TG'],
         ['CA'],['CT'],['CC'],['CG'],
         ['GA'],['GT'],['GG'],['GC']]

def GenomeFrequency(Fa,Pattern):
    gsize=0
    for i in Fa.keys():
        gsize+=len(str(Fa[i]))
    Count=0
    for pattern in Pattern:
        countlist=[str(Fa[chr]).count(pattern) for chr in Fa.keys()]
        countarray=np.array(countlist)
        Sum=countarray.sum()
        Count+=Sum
    frequency=Count/gsize
    return(frequency)

"""
region: {'chr1':{'1643':(13,15),'1644':(12,15),...},
         'chr2':{'...'},...}
"""
def Revcomplement(Seq):
    out = ''
    rule={'A':'T','T':'A','G':'C','C':'G','N':'N','a':'t','t':'a','g':'c','c':'g','n':'n'}
    for i in Seq:
        out+=rule[i]
    return out[::-1]

#正负链如何考虑;基因组有无N，是否大小写  编辑到此处2020.5.28 14:05
def bedDinucStart(bedfile,fasta,windows,up,down,pattern):
    dinucList=[]
    bedf=open(bedfile,'r')
    dinucCount=0
    for line in bedf:
        line=line.strip().split()
        sublist=[]
        if line[5] is "+":
            left=int(line[1])-up
            right=int(line[1])+down
            numofbin = (right-left)//windows #整除
            seq=fasta.sequence({'chr':line[0],'start':left,'stop':right+1},one_based=False)
            #先分bin 取每个bin的起始坐标
            for k in range(0,len(seq)-1,windows):
                #取出每个bin的序列
                subseq=seq[k:k+windows+1]
                subrevseq=Revcomplement(subseq)
                #单碱基遍历数二核苷酸
                for j in range(len(subseq)-1):
                    dinucFor=subseq[j:j+2]
                    dinucRev=subrevseq[j:j+2]
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(windows*2)
                sublist.append(dinucDensity)
                dinucCount=0
            #若位点的长度超出了基因组的end,用0补齐
            if len(sublist)<numofbin:
                for _ in range(numofbin-len(sublist)):sublist.append(0)
            dinucList.append(sublist)
        else:
            left=int(line[2])-down
            right=int(line[2])+up
            numofbin = (right-left)//windows #整除
            seq=fasta.sequence({'chr':line[0],'start':left,'stop':right+1},one_based=False)
            #先分bin 取每个bin的起始坐标
            for k in range(len(seq)-1,0,-windows):
                #取出每个bin的序列
                subseq=seq[k-windows:k+1]
                subrevseq=Revcomplement(subseq)
                #单碱基遍历所有二核苷酸
                for j in range(len(subseq)-1):
                    dinucFor=subseq[j:j+2]
                    dinucRev=subrevseq[j:j+2]
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(windows*2)
                sublist.append(dinucDensity)
                dinucCount=0
            #若位点的长度超出了基因组的end,用0补齐
            if len(sublist)<numofbin:
                for _ in range(numofbin-len(sublist)):sublist.append(0)
            dinucList.append(sublist)
    return dinucList
    bedfile.close()

def bedDinucEnd(bedfile,fasta,windows,up,down,pattern):
    dinucList=[]
    bedf=open(bedfile,'r')
    dinucCount=0
    for line in bedf:
        line=line.strip().split()
        sublist=[]
        if line[5] is "+":
            left=int(line[2])-up
            right=int(line[2])+down
            numofbin = (right-left)//windows #整除
            seq=fasta.sequence({'chr':line[0],'start':left,'stop':right+1},one_based=False)
            #先分bin 取每个bin的起始坐标
            for k in range(0,len(seq)-1,windows):
                #取出每个bin的序列
                subseq=seq[k:k+windows+1]
                subrevseq=Revcomplement(subseq)
                for j in range(len(subseq)-1):
                    dinucFor=subseq[j:j+2]
                    dinucRev=subrevseq[j:j+2]
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(windows*2)
                sublist.append(dinucDensity)
                dinucCount=0
            #若位点的长度超出了基因组的end,用0补齐
            if len(sublist)<numofbin:
                for _ in range(numofbin-len(sublist)):sublist.append(0)
            dinucList.append(sublist)
        else:
            left=int(line[1])-down
            right=int(line[1])+up
            numofbin = (right-left)//windows #整除
            seq=fasta.sequence({'chr':line[0],'start':left,'stop':right+1},one_based=False)
            #先分bin 取每个bin的起始坐标
            for k in range(len(seq)-1,0,-windows):
                #取出每个bin的序列
                subseq=seq[k-windows:k+1]
                subrevseq=Revcomplement(subseq)
                #单碱基遍历所有二核苷酸
                for j in range(len(subseq)-1):
                    dinucFor=subseq[j:j+2]
                    dinucRev=subrevseq[j:j+2]
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(windows*2)
                sublist.append(dinucDensity)
                dinucCount=0
            #若位点的长度超出了基因组的end,用0补齐
            if len(sublist)<numofbin:
                for _ in range(numofbin-len(sublist)):sublist.append(0)
            dinucList.append(sublist)
    return dinucList
    bedfile.close()

def bedDinucBody(bedfile,fasta,windows,up,down,Blocknum,pattern):
    dinucList=[]
    bedf=open(bedfile,'r')
    dinucCount=0
    for line in bedf:
        line=line.strip().split()
        start=int(line[1])
        end=int(line[2])
        sublist=[]
        reside=Blocknum-int((end-start)%Blocknum )
        blockSize=int((end-start+reside)/Blocknum)
        if line[5] is "+":
            seq1=fasta.sequence({'chr':line[0],'start':start-up,'stop':start+1},one_based=False)
            seq2=fasta.sequence({'chr':line[0],'start':start,'stop':start+reside+1},one_based=False)
            seq3=fasta.sequence({'chr':line[0],'start':end,'stop':end+down+1},one_based=False)
            numofbin = ((start-up+end+down)//windows)+Blocknum #整除
            #先分bin 取每个bin的起始坐标
            for k in range(0,len(seq1)-1,windows):
                #取出每个bin的序列
                subseq=seq1[k:k+windows+1]
                subrevseq=Revcomplement(subseq)
                for j in range(len(subseq)-1):
                    dinucFor=subseq[j:j+2]
                    dinucRev=subrevseq[j:j+2]
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(windows*2)
                sublist.append(dinucDensity)
                dinucCount=0
            for x in range(0,len(seq2)-1,blockSize):
                #取出每个bin的序列
                subseq=seq2[x:x+blockSize+1]
                subrevseq=Revcomplement(subseq)
                for y in range(len(subseq)-1):
                    dinucFor=subseq[y:y+2]
                    dinucRev=subrevseq[y:y+2]
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(blockSize*2)
                sublist.append(dinucDensity)
                dinucCount=0
            for u in range(0,len(seq3)-1,windows):
                #取出每个bin的序列
                subseq=seq3[u:u+windows+1]
                subrevseq=Revcomplement(subseq)
                for v in range(len(subseq)-1):
                    dinucFor=subseq[v:v+2]
                    dinucRev=subrevseq[v:v+2]
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(windows*2)
                sublist.append(dinucDensity)
                dinucCount=0
            if len(sublist)<numofbin:
                for _ in range(numofbin-len(sublist)):sublist.append(0)
            dinucList.append(sublist)
        else:
            seq1=fasta.sequence({'chr':line[0],'start':end,'stop':end+up+1},one_based=False)
            seq2=fasta.sequence({'chr':line[0],'start':start-reside,'stop':end+1},one_based=False)
            seq3=fasta.sequence({'chr':line[0],'start':start-down,'stop':start+1},one_based=False)
            numofbin = ((start-down+end+up)//windows)+Blocknum #整除
            for k in range(len(seq1)-1,0,-windows):
                #取出每个bin的序列
                subseq=seq1[k-windows:k+1]
                subrevseq=Revcomplement(subseq)
                for j in range(len(subseq)-1):
                    dinucFor=subseq[j:j+2]
                    dinucRev=subrevseq[j:j+2]
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(windows*2)
                sublist.append(dinucDensity)
                dinucCount=0
            for x in range(len(seq2)-1,0,-blockSize):
                #取出每个bin的序列
                subseq=seq2[x-blockSize:x+1]
                subrevseq=Revcomplement(subseq)
                for y in range(len(subseq)-1):
                    dinucFor=subseq[y:y+2]
                    dinucRev=subrevseq[y:y+2]
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(blockSize*2)
                sublist.append(dinucDensity)
                dinucCount=0
            for u in range(len(seq3)-1,0,-windows):
                #取出每个bin的序列
                subseq=seq3[u-windows:u+1]
                subrevseq=Revcomplement(subseq)
                for v in range(len(subseq)-1):
                    dinucFor=subseq[v:v+2]
                    dinucRev=subrevseq[v:v+2]
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(windows*2)
                sublist.append(dinucDensity)
                dinucCount=0
            if len(sublist)<numofbin:
                for _ in range(numofbin-len(sublist)):sublist.append(0)
            dinucList.append(sublist)
    return dinucList
    bedfile.close()
    
#readdict {'chr1':{'1045':546,''1067:'453',...},...}
#sublist [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8]

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

def plotfacet(df,xtick,xticklabel,ylabel,outfigname):
    Pal="Set1"
    u = df.Dinucletide.unique()
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
        #遍历每个分面
        for Dinucletide, ax in zip(u, axes):
            subdf=df[df.Dinucletide==Dinucletide]
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
            g.set(title=Dinucletide)
            g.set_xticks(xtick)
            g.set_xticklabels(xticklabel)
            g.set_xlabel("")
            g.set_ylabel(ylabel)
            g.figure.set_size_inches(ncol*4.5,nrow*3.5)
        #plt.tight_layout()
        fig.set_tight_layout(True)
    else:
        fig, axes = plt.subplots(ncols=1)
        Dinucletide=u[0]
        df_new=dataframe_Interpolation(df,'Bedlabel')
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
        g.set(title=Dinucletide)
        g.set_xticks(xtick)
        g.set_xticklabels(xticklabel)
        g.set_xlabel("")
        g.set_ylabel(ylabel)
        g.figure.set_size_inches(4.5,3.5)
        #plt.tight_layout()
        fig.set_tight_layout(True)
    plt.savefig(outfigname)
    
if "--reference-point" in Option:
    #ylength=int(up/windows)
    left=int(0)
    medium=int(up/windows)
    right=int((up+down)/windows)
    NumOfBin = int((int(up+down))/windows)
    if "-r" in Option and str(referencePoint)=="start":
        out_profile=open(outf[0],"w")
#         out_matrix=open(outf[1],"w")
        out_figure=open(figurename,"w")
        figurelist=[]
        fa=Fasta(fasta,key_fn=lambda key: key.split()[0])
        print("Bedlabel","Dinucletide","bin","value","norm_value",sep="\t",end="\n",file=out_profile)
        for Bed in range(len(bedlist)):
            for dinuc,pattern in zip(Dinuc,Pattern):
                values=bedDinucStart(bedlist[Bed],fa,windows,up,down,pattern)
                genomefrequency=GenomeFrequency(fa,pattern)
                tup1=np.array(values)
                tup1=tup1.astype(float)
                tup2=tup1.mean(axis=0)
                tup3=tup2/genomefrequency  
                i=1
                for binValue,normValue in zip(tup2,tup3):
                    print(bedlabel[Bed],dinuc,i,binValue,normValue,sep="\t",end="\n",file=out_profile)
                    i+=1
        #out_matrix.close()
        out_profile.close()
        dataframe=pd.read_csv(outf[0],sep='\t')
        Xticks=[left,medium,right]
        Xticklabel=["-"+L+"kb","TSS",R+"kb"]
        Ylabel="Relative dinucleotide density"
        plotfacet(dataframe,Xticks,Xticklabel,Ylabel,figurename)
    elif "-r" in Option and str(referencePoint)=="end":
        out_profile=open(outf[0],"w")
#         out_matrix=open(outf[1],"w")
        out_figure=open(figurename,"w")
        figurelist=[]
        fa=Fasta(fasta,key_fn=lambda key: key.split()[0])
        print("Bedlabel","Dinucletide","bin","value","norm_value",sep="\t",end="\n",file=out_profile)
        for Bed in range(len(bedlist)):
            for dinuc,pattern in zip(Dinuc,Pattern):
                values=bedDinucEnd(bedlist[Bed],fa,windows,up,down,pattern)
                genomefrequency=GenomeFrequency(fa,pattern)
                tup1=np.array(values)
                tup1=tup1.astype(float)
                tup2=tup1.mean(axis=0)
                tup3=tup2/genomefrequency
                i=1
                for binValue,normValue in zip(tup2,tup3):
                    print(bedlabel[Bed],dinuc,i,binValue,normValue,sep="\t",end="\n",file=out_profile)
                    i+=1
        #out_matrix.close()
        out_profile.close()
        dataframe=pd.read_csv(outf[0],sep='\t')
        Xticks=[left,medium,right]
        Xticklabel=["-"+L+"kb","TSS",R+"kb"]
        Ylabel="Relative dinucleotide density"
        plotfacet(dataframe,Xticks,Xticklabel,Ylabel,figurename)
    else:
        pass

elif "--scale-regions" in Option:
    left=int(0)
    mediumA=int(up/windows)
    mediumB=mediumA+BlockNum
    right=int(down/windows)+mediumA+BlockNum
    out_profile=open(outf[0],"w")
    #out_matrix=open(outf[1],"w")
    out_figure=open(figurename,"w")
    fa=Fasta(fasta,key_fn=lambda key: key.split()[0])
    print("Bedlabel","Dinucletide","bin","value","norm_value",sep="\t",end="\n",file=out_profile)
    for Bed in range(len(bedlist)):
        for dinuc,pattern in zip(Dinuc,Pattern):
            values=bedDinucBody(bedlist[Bed],fa,windows,up,down,BlockNum,pattern)
            genomefrequency=GenomeFrequency(fa,pattern)
            tup1=np.array(values)
            tup1=tup1.astype(float)
            tup2=tup1.mean(axis=0)
            tup3=tup2/genomefrequency
            i=1
            for binValue,normValue in zip(tup2,tup3):
                print(bedlabel[Bed],dinuc,i,binValue,normValue,sep="\t",end="\n",file=out_profile)
                i+=1
    #out_matrix.close()
    out_profile.close()
    dataframe=pd.read_csv(outf[0],sep='\t')
    Xticks=[left,mediumA,mediumB,right]
    Xticklabel=["-"+L+"kb","TSS","TTS",R+"kb"]
    Ylabel="Relative dinucleotide density"
    plotfacet(dataframe,Xticks,Xticklabel,Ylabel,figurename)
else:
    pass
