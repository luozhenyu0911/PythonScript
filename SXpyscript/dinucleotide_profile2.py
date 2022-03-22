# -*- coding: UTF-8 -*-
from __future__ import print_function
from __future__ import division
from pyfasta import Fasta
import sys, getopt
import gzip
import numpy as np
import matplotlib.pyplot as plt
#from scipy.interpolate import interp1d
from scipy.interpolate import make_interp_spline, BSpline
plt.switch_backend('PDF')
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
    print("\t"+"-p"+"    <pattern {Any Dinucleotide, for example [CG] or [AT] or [WW] or [SS] and so on.}>")
    #print("\t"+"-s"+"    <strand {forward or reverse}>"+"\t"+"which DNA strand.")
    print("\t"+"-r"+"    <referencePoint {start or end}>"+"\t"+"The reference point for the plotting could be either the region start (TSS), the region end (TES) of the region.")
    print("\t"+"-u"+"    <upstream distance>"+"\t"+"Distance upstream of the reference-point selected.")
    print("\t"+"-d"+"    <downstream distance>"+"\t"+"Distance downstream of the reference-point selected.")
    print("\t"+"-w"+"    <window size>"+"\t"+"Length, in bases, binSize for averaging the score over the regions length.")
    print("\t"+"-N"+"    <Equally fragment number>"+"\t"+"Equally divided into the region into several fragments.")
    print("\t"+"-l"+"    <legend names>"+"\t"+"Name of the plot legend. (Support for multiple legends)")
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
    if opt_name in ('-s'):
        strand = opt_value
    if opt_name in ('-p'):
        Pattern = opt_value  
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
    if opt_name in ('-l'):
        legend = opt_value
        legendname=legend.split(",")
    if opt_name in ('-o'):
        outputfile = opt_value
        outf=outputfile.split(",")
    if opt_name not in ('-h,--help,-v,--version,-b,-G,-r,-u,-d,-w,-N,-o,-l,-F,-S,-c,--scale-regions,--reference-point'):
        message()
    #print("Usage: python tmp1.py <--scale-regions or --reference-point> -b <CGmapfiles> -G <bedfiles> -r <referencePoint {start or end}> -u <upstream_distance> -d <downstream_distance> -N <Equal_block_number> -l <legend_names> -o <outputfiles average.matrix> -F <Figure_name>")
Option={}
for option in opts:
    Option[option[0]]=option[1]

Dinuc={'AA':['AA'],
       'AT':['AT'],
       'AC':['AC'],
       'AG':['AG'],
       'TA':['TA'],
       'TT':['TT'],
       'TC':['TC'],
       'TG':['TG'],
       'CA':['CA'],
       'CT':['CT'],
       'CC':['CC'],
       'CG':['CG'],
       'GA':['GA'],
       'GT':['GT'],
       'GG':['GG'],
       'GC':['GC'],
       'WW':['AA','AT','TA','TT'],
       'SS':['GG','GC','CG','CC']
      }

def GenomeFrequency(Fa,Dinucdict,pattern):
    Dinuclist=Dinucdict[pattern]
    gsize=0
    for i in Fa.keys():
        gsize+=len(str(Fa[i]))
    Count=0
    for dinuc in Dinuclist:
        countlist=[str(Fa[chr]).count(dinuc) for chr in Fa.keys()]
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
            for k in range(left,right,windows):
                for j in range(k,k+windows):
                    dinucFor=fasta.sequence({'chr':line[0],'start':j,'stop':j+2},one_based=False)
                    dinucRev=fasta.sequence({'chr':line[0],'start':j,'stop':j+2,'strand':'-'},one_based=False)
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(windows*2)
                sublist.append(dinucDensity)
                dinucCount=0
            dinucList.append(sublist)
        else:
            left=int(line[2])-down
            right=int(line[2])+up
            for k in range(right,left,-windows):
                for j in range(k-windows,k):
                    dinucFor=fasta.sequence({'chr':line[0],'start':j,'stop':j+2},one_based=False)
                    dinucRev=fasta.sequence({'chr':line[0],'start':j,'stop':j+2,'strand':'-'},one_based=False)
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(windows*2)
                sublist.append(dinucDensity)
                dinucCount=0
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
            for k in range(left,right,windows):
                for j in range(k,k+windows):
                    dinucFor=fasta.sequence({'chr':line[0],'start':j,'stop':j+2},one_based=False)
                    dinucRev=fasta.sequence({'chr':line[0],'start':j,'stop':j+2,'strand':'-'},one_based=False)
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(windows*2)
                sublist.append(dinucDensity)
                dinucCount=0
            dinucList.append(sublist)
        else:
            left=int(line[1])-down
            right=int(line[1])+up
            for k in range(right,left,-windows):
                for j in range(k-windows,k):
                    dinucFor=fasta.sequence({'chr':line[0],'start':j,'stop':j+2},one_based=False)
                    dinucRev=fasta.sequence({'chr':line[0],'start':j,'stop':j+2,'strand':'-'},one_based=False)
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(windows*2)
                sublist.append(dinucDensity)
                dinucCount=0
            dinucList.append(sublist)
    return dinucList
    bedfile.close()

def bedDinucBody(bedfile,fasta,windows,up,down,Blocknum,pattern):
    dinucList=[]
    bedf=open(bedfile,'r')
    dinucCount=0
    for line in bedf:
        line=line.strip().split()
        sublist=[]
        reside=Blocknum-int(((int(line[2])-int(line[1]))%Blocknum))
        blockSize=int(((int(line[2])-int(line[1]))+reside)/Blocknum)
        if line[5] is "+":
            for k in range(int(line[1])-up,int(line[1]),windows):
                for j in range(k,k+windows):
                    dinucFor=fasta.sequence({'chr':line[0],'start':j,'stop':j+2},one_based=False)
                    dinucRev=fasta.sequence({'chr':line[0],'start':j,'stop':j+2,'strand':'-'},one_based=False)
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(windows*2)
                sublist.append(dinucDensity)
                dinucCount=0
            for x in range(int(line[1]),int(line[2])+reside,blockSize):
                for y in range(x,x+blockSize):
                    dinucFor=fasta.sequence({'chr':line[0],'start':y,'stop':y+2},one_based=False)
                    dinucRev=fasta.sequence({'chr':line[0],'start':y,'stop':y+2,'strand':'-'},one_based=False)
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(blockSize*2)
                sublist.append(dinucDensity)
                dinucCount=0
            for u in range(int(line[2]),int(line[2])+down,windows):
                for v in range(u,u+windows):
                    dinucFor=fasta.sequence({'chr':line[0],'start':v,'stop':v+2},one_based=False)
                    dinucRev=fasta.sequence({'chr':line[0],'start':v,'stop':v+2,'strand':'-'},one_based=False)
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(windows*2)
                sublist.append(dinucDensity)
                dinucCount=0
            dinucList.append(sublist)
        else:
            for k in range(int(line[2])+up,int(line[2]),-windows):
                for j in range(k-windows,k):
                    dinucFor=fasta.sequence({'chr':line[0],'start':j,'stop':j+2},one_based=False)
                    dinucRev=fasta.sequence({'chr':line[0],'start':j,'stop':j+2,'strand':'-'},one_based=False)
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(windows*2)
                sublist.append(dinucDensity)
                dinucCount=0
            for x in range(int(line[2]),int(line[1])-reside,-blockSize):
                for y in range(x-blockSize,x):
                    dinucFor=fasta.sequence({'chr':line[0],'start':y,'stop':y+2},one_based=False)
                    dinucRev=fasta.sequence({'chr':line[0],'start':y,'stop':y+2,'strand':'-'},one_based=False)
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(blockSize*2)
                sublist.append(dinucDensity)
                dinucCount=0
            for u in range(int(line[1]),int(line[1])-down,-windows):
                for v in range(u-windows,u):
                    dinucFor=fasta.sequence({'chr':line[0],'start':v,'stop':v+2},one_based=False)
                    dinucRev=fasta.sequence({'chr':line[0],'start':v,'stop':v+2,'strand':'-'},one_based=False)
                    if dinucFor in pattern:
                        dinucCount+=1
                    if dinucRev in pattern:
                        dinucCount+=1
                dinucDensity=dinucCount/(windows*2)
                sublist.append(dinucDensity)
                dinucCount=0
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
#         header = ["label","num"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
#         print(*header,sep="\t",end="\n",file=out_matrix)
        fa=Fasta(fasta,key_fn=lambda key: key.split()[0])
        genomefrequency=GenomeFrequency(fa,Dinuc,Pattern)
        for Bed in range(len(bedlist)):
            values=bedDinucStart(bedlist[Bed],fa,windows,up,down,Dinuc[Pattern])
            tup1=np.array(values)
            tup1=tup1.astype(float)
            tup2=tup1.mean(axis=0)
            tup3=tup2/genomefrequency
            figurelist.append(tup3)     
#             for index in range(len(values)):
#                 print(legendname[Bed],str(index+1),*values[index],sep="\t",end="\n",file=out_matrix)
        print("region","model","value",sep="\t",end="\n",file=out_profile)
        N = 0
        for lis in figurelist:
            n = 0
            for l in lis:
                print(n*windows,legendname[N],l,sep="\t",end="\n",file=out_profile)
                n += 1
            N += 1
        for line in figurelist:
            plotsmooth(line)
        plt.legend(legendname)
        plt.ylabel("Relative "+Pattern+" dinucleotide density")
        plt.xticks([left,medium,right],["-"+L+"kb","TSS",R+"kb"])
        plt.savefig(out_figure)    
        out_profile.close()
        #out_matrix.close()
        out_figure.close()
    elif "-r" in Option and str(referencePoint)=="end":
        out_profile=open(outf[0],"w")
#         out_matrix=open(outf[1],"w")
        out_figure=open(figurename,"w")
        figurelist=[]
#         header = ["label","num"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
#         print(*header,sep="\t",end="\n",file=out_matrix)
        fa=Fasta(fasta,key_fn=lambda key: key.split()[0])
        genomefrequency=GenomeFrequency(fa,Dinuc,Pattern)
        for Bed in range(len(bedlist)):
            values=bedDinucEnd(bedlist[Bed],fa,windows,up,down,Dinuc[Pattern])
            tup1=np.array(values)
            tup1=tup1.astype(float)
            tup2=tup1.mean(axis=0)
            tup3=tup2/genomefrequency
            figurelist.append(tup3)
#             for index in range(len(values)):
#                 print(legendname[Bed],str(index+1),*values[index],sep="\t",end="\n",file=out_matrix)
        print("region","model","value",sep="\t",end="\n",file=out_profile)
        N = 0
        for lis in figurelist:
            n = 0
            for l in lis:
                print(n*windows,legendname[N],l,sep="\t",end="\n",file=out_profile)
                n += 1
            N += 1
        for line in figurelist:
            plotsmooth(line)
        plt.legend(legendname)
        plt.ylabel("Relative "+Pattern+" dinucleotide density")
        plt.xticks([left,medium,right],["-"+L+"kb","TTS",R+"kb"])
        plt.savefig(out_figure)     
        out_figure.close()
        out_profile.close()
        #out_matrix.close()
    else:
        pass
#print("Usage: python tmp1.py <<--scale-regions or --reference-point>> -b <CGmapfile> -G <bedfile> -g <genelist> -r <referencePoint {start or end}> -u <upstream_distance> -d <downstream_distance> -w <window_size> -bn <Equal_block_number> -o <outputfile>")

elif "--scale-regions" in Option:
    left=int(0)
    mediumA=int(up/windows)
    mediumB=mediumA+BlockNum
    right=int(down/windows)+mediumA+BlockNum
    out_profile=open(outf[0],"w")
    out_matrix=open(outf[1],"w")
    out_figure=open(figurename,"w")
    figurelist=[]
#     header = ["label","num"]+["bin"+str(bin+1) for bin in range(NumOfBin)]
#     print(*header,sep="\t",end="\n",file=out_matrix)
    fa=Fasta(fasta,key_fn=lambda key: key.split()[0])
    genomefrequency=GenomeFrequency(fa,Dinuc,Pattern)
    for Bed in range(len(bedlist)):
        values=bedDinucBody(bedlist[Bed],fa,windows,up,down,BlockNum,Dinuc[Pattern])
        tup1=np.array(values)
        tup1=tup1.astype(float)
        tup2=tup1.mean(axis=0)
        tup3=tup2/genomefrequency
        figurelist.append(tup3)
#         for index in range(len(values)):
#             print(legendname[Bed],str(index+1),*values[index],sep="\t",end="\n",file=out_matrix)
    print("region","model","value",sep="\t",end="\n",file=out_profile)
    N = 0
    for lis in figurelist:
        n = 0
        for l in lis:
            print(n*windows,legendname[N],l,sep="\t",end="\n",file=out_profile)
            n += 1
        N += 1
    for line in figurelist:
        plotsmooth(line)
    plt.legend(legendname)
    plt.ylabel("Relative "+Pattern+" dinucleotide density")
    plt.xticks([left,mediumA,mediumB,right],["-"+L+"kb","TSS","TTS",R+"kb"])
    plt.savefig(out_figure)
    out_figure.close()
    out_profile.close()
    #out_matrix.close()
else:
    pass
