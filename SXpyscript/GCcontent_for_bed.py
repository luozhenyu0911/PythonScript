from __future__ import print_function
from __future__ import division
from pyfasta import Fasta
import sys, getopt
import gzip
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
plt.switch_backend('PDF')
##@jit

def message():
    print("\n"+"Usage: python **.py --reference-point or --scale-regions [options]"+"\n")
    print("Description:"+"\n"+"\t"+"\t"+"This tool creates CG content ratio profile over sets of genomic regions."+"\n")
    print("Options:"+"\n")
    print("\t"+"-v, "+"--version"+"\t"+"show program's version number and exit")
    print("\t"+"-h, "+"--help"+"\t"+"show this help message and exit.")
    print("\t"+"--reference-point or --scale-regions"+"\t"+"Reference-point refers to a position within a BED region(e.g., the starting point). In this mode, only those genomicpositions before (upstream) and/or after (downstream) of the reference point will be plotted."+"In the scale-regions mode, all regions in the BED file are stretched or shrunken to the length (in bases) indicated by the user.")
    print("\t"+"-b"+"    <Reference>"+"\t"+"Reference genome fasta files")
    print("\t"+"-G"+"    <bedfiles>"+"\t"+"data files(Support for multiple files, Separated by commas), should be original bed files.")
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

opts,args = getopt.getopt(sys.argv[1:],'-h-v-b:-G:-u:-r:-d:-w:-N:-o:-l:-F:',['help','version','scale-regions','reference-point'])
for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        message()
    if opt_name in ('-v','--version'):
        print("**.py 1.0 2019-10-16")
    if opt_name in ('-b'):
        fasta = opt_value
        #CGmapfile = opt_value
        #CGmaplist=CGmapfile.split(",")
    if opt_name in ('-G'):
        bedfile = opt_value
        bedlist=bedfile.split(",")
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

def bedGCStart(bedfile,fasta,windows,up,down):
    GCList=[]
    bedf=open(bedfile,'r')
    for line in bedf:
        line=line.strip().split()
        sublist=[]
        if line[5] is "+":
            left=int(line[1])-up
            right=int(line[1])+down
            for k in range(left,right,windows):
                seq=fasta.sequence({'chr':line[0],'start':k,'stop':k+windows},one_based=False)
                G=seq.count('G')+seq.count('g')
                C=seq.count('C')+seq.count('c')
                GC=(G+C)/windows
                sublist.append(GC)
            GCList.append(sublist)
        else:
            left=int(line[2])-up
            right=int(line[2])+down
            for k in range(right,left,-windows):
                seq=fasta.sequence({'chr':line[0],'start':k-windows,'stop':k},one_based=False)
                G=seq.count('G')+seq.count('g')
                C=seq.count('C')+seq.count('c')
                GC=(G+C)/windows
                sublist.append(GC)
            GCList.append(sublist)
    return GCList
    bedfile.close()

def bedGCEnd(bedfile,fasta,windows,up,down):
    GCList=[]
    bedf=open(bedfile,'r')
    for line in bedf:
        line=line.strip().split()
        sublist=[]
        if line[5] is "+":
            left=int(line[2])-up
            right=int(line[2])+down
            for k in range(left,right,windows):
                seq=fasta.sequence({'chr':line[0],'start':k,'stop':k+windows},one_based=False)
                G=seq.count('G')+seq.count('g')
                C=seq.count('C')+seq.count('c')
                GC=(G+C)/windows
                sublist.append(GC)
            GCList.append(sublist)
        else:
            left=int(line[1])-up
            right=int(line[1])+down
            for k in range(right,left,-windows):
                seq=fasta.sequence({'chr':line[0],'start':k-windows,'stop':k},one_based=False)
                G=seq.count('G')+seq.count('g')
                C=seq.count('C')+seq.count('c')
                GC=(G+C)/windows
                sublist.append(GC)
            GCList.append(sublist)
    return GCList
    bedfile.close()

def bedGCBody(bedfile,fasta,windows,up,down,Blocknum):
    CGList=[]
    bedf=open(bedfile,'r')
    for line in bedf:
        line=line.strip().split()
        sublist=[]
        reside=Blocknum-int(((int(line[2])-int(line[1]))%Blocknum))
        blockSize=int(((int(line[2])-int(line[1]))+reside)/Blocknum)
        if line[5] is "+":
            for k in range(int(line[1])-up,int(line[1]),windows):
                seq=fasta.sequence({'chr':line[0],'start':k,'stop':k+windows},one_based=False)
                G=seq.count('G')+seq.count('g')
                C=seq.count('C')+seq.count('c')
                GC=(G+C)/windows
                sublist.append(GC)
            for x in range(int(line[1]),int(line[2])+reside,blockSize):
                seq=fasta.sequence({'chr':line[0],'start':x,'stop':x+blockSize},one_based=False)
                G=seq.count('G')+seq.count('g')
                C=seq.count('C')+seq.count('c')
                GC=(G+C)/blockSize
                sublist.append(GC)
            for u in range(int(line[2]),int(line[2])+down,windows):
                seq=fasta.sequence({'chr':line[0],'start':u,'stop':u+windows},one_based=False)
                G=seq.count('G')+seq.count('g')
                C=seq.count('C')+seq.count('c')
                GC=(G+C)/windows
                sublist.append(GC)
            CGList.append(sublist)
        else:
            for k in range(int(line[2])+up,int(line[2]),-windows):
                seq=fasta.sequence({'chr':line[0],'start':k-windows,'stop':k},one_based=False)
                G=seq.count('G')+seq.count('g')
                C=seq.count('C')+seq.count('c')
                GC=(G+C)/windows
                sublist.append(GC)
            for x in range(int(line[2]),int(line[1])-reside,-blockSize):
                seq=fasta.sequence({'chr':line[0],'start':x-blockSize,'stop':x},one_based=False)
                G=seq.count('G')+seq.count('g')
                C=seq.count('C')+seq.count('c')
                GC=(G+C)/blockSize
                sublist.append(GC)
            for u in range(int(line[1]),int(line[1])-down,-windows):
                seq=fasta.sequence({'chr':line[0],'start':u-windows,'stop':u},one_based=False)
                G=seq.count('G')+seq.count('g')
                C=seq.count('C')+seq.count('c')
                GC=(G+C)/windows
                sublist.append(GC)
            CGList.append(sublist)
    return CGList
    bedfile.close()
    
#readdict {'chr1':{'1045':546,''1067:'453',...},...}
#sublist [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8]


def plotsmooth(Values):
    n=len(Values)
    num=range(0,n)
    num=np.array(num)
    xnew=np.linspace(num.min(),num.max(),5*n)
    func=interp1d(num,Values,kind="cubic")
    ynew=func(xnew)
    plt.plot(xnew,ynew)


if "--reference-point" in Option:
    Ylength=int(up+down)
    ylength=int(Ylength/windows)
    left=int(0)
    medium=int(ylength/2)
    right=int(ylength)
    if "-r" in Option and str(referencePoint)=="start":
        out_profile=open(outf[0],"w")
#         out_matrix=open(outf[1],"w")
        out_figure=open(figurename,"w")
        figurelist=[]
        fa=Fasta(fasta,key_fn=lambda key: key.split()[0])
        for Bed in range(len(bedlist)):
            values=bedGCStart(bedlist[Bed],fa,windows,up,down)
            tup1=np.array(values)
            tup1=tup1.astype(float)
            tup2=tup1.mean(axis=0)
            figurelist.append(tup2)
            tmpnumber = len(values[0])
#             header = ["a" for _ in range(tmpnumber+2)]
#             print(*header,sep="\t",end="\n",file=out_matrix)
#             for tmp1 in values:
#                 print("a","a",*tmp1,sep="\t",end="\n",file=out_matrix)
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
        plt.ylabel("GC content")
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
        fa=Fasta(fasta,key_fn=lambda key: key.split()[0])
        for Bed in range(len(bedlist)):
            values=bedGCEnd(bedlist[Bed],fa,windows,up,down)
            tup1=np.array(values)
            tup1=tup1.astype(float)
            tup2=tup1.mean(axis=0)
            figurelist.append(tup2)
            tmpnumber = len(values[0])
#             header = ["a" for _ in range(tmpnumber+2)]
#             print(*header,sep="\t",end="\n",file=out_matrix)
#             for tmp1 in values:
#                 print("a","a",*tmp1,sep="\t",end="\n",file=out_matrix)
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
        plt.ylabel("GC content")
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
#     out_matrix=open(outf[1],"w")
    out_figure=open(figurename,"w")
    figurelist=[]
    fa=Fasta(fasta,key_fn=lambda key: key.split()[0])
    for Bed in range(len(bedlist)):
        values=bedGCBody(bedlist[Bed],fa,windows,up,down,BlockNum)
        tup1=np.array(values)
        tup1=tup1.astype(float)
        tup2=tup1.mean(axis=0)
        figurelist.append(tup2)
        tmpnumber = len(values[0])
#         header = ["a" for _ in range(tmpnumber+2)]
#         print(*header,sep="\t",end="\n",file=out_matrix)
#         for tmp1 in values:
#             print("a","a",*tmp1,sep="\t",end="\n",file=out_matrix)
    print("region","model","value",sep="\t",end="\n",file=out_profile)
    N = 0
    for line in figurelist:
        for bin in range(len(line)):
            print(bin+1,legendname[N],line[bin],sep="\t",end="\n",file=out_profile)
        N += 1
    for line in figurelist:
        plotsmooth(line)
    plt.legend(legendname)
    plt.ylabel("GC content")
    plt.xticks([left,mediumA,mediumB,right],["-"+L+"kb","TSS","TTS",R+"kb"])
    plt.savefig(out_figure)
    out_figure.close()
    out_profile.close()
    #out_matrix.close()
else:
    pass
