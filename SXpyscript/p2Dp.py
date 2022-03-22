# -*- coding: UTF-8 -*-
from __future__ import print_function
from __future__ import division
#from numba import jit
import sys, getopt
import pysam 
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
plt.switch_backend('PDF')
##@jit

def message():
    print("\n"+"Usage: python **.py --reference-point[options]"+"\n")
    print("Description:"+"\n"+"\t"+"\t"+"This tool creates scores and profile plot over sets of genomic regions. Typically, these regions are genes, but any other regions defined in BED will work."+"\n")
    print("Options:"+"\n")
    print("\t"+"-v, "+"--version"+"\t"+"show program's version number and exit")
    print("\t"+"-h, "+"--help"+"\t"+"show this help message and exit.")
    print("\t"+"-b"+"    <bamfiles>"+"\t"+"data files, should be sorted bam format.")
    print("\t"+"-G"+"    <bedfiles>"+"\t"+"data files, should be original bed files.")
    print("\t"+"-r"+"    <referencePoint {start or end}>"+"\t"+"The reference point for the plotting could be either the region start (TSS), the region end (TES) of the region.")
    print("\t"+"-u"+"    <upstream distance>"+"\t"+"Distance upstream of the reference-point selected.")
    print("\t"+"-d"+"    <downstream distance>"+"\t"+"Distance downstream of the reference-point selected.")
    print("\t"+"-w"+"    <window size>"+"\t"+"Length, in bases, binSize for averaging the score over the regions length.")
    #print("\t"+"-l"+"    <legend names>"+"\t"+"Name of the plot legend. (Support for multiple legends)")
    print("\t"+"-o"+"    <outputfiles average matrix>"+"\t"+"File name to save the average matrix file.")
    #print("\t"+"-F"+"    <Figure name>"+"\t"+"Name of the plot name. (plot file format is PDF)"+"\n")
    print("An example usage is:"+"\n")
    print("\t"+"python **.py --reference-point -b **.bam,**.bam... -G **.bed,**.bed... -r start/end -u 1000 -d 1000 -w 50 -l **,** -o **_matrix -F **.pdf"+"\n")
    return

opts,args = getopt.getopt(sys.argv[1:],'-h-v-b:-G:-g:-u:-r:-d:-w:-N:-o:-l:-F:-c:',['help','version','scale-regions','reference-point'])
for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        message()
    if opt_name in ('-v','--version'):
        print("**.py 1.0 2020-05-09")
    if opt_name in ('-b'):
        bamfile = opt_value
        #bamlist=bamfile.split(",")
    if opt_name in ('-G'):
        bedfile = opt_value
        #bedlist=bedfile.split(",")
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
        wdsize=int(window)
    if opt_name in ('-F'):
        figurename = opt_value
    if opt_name in ('-l'):
        legend = opt_value
        legendname=legend.split(",")
    if opt_name in ('-o'):
        outf=opt_value
    if opt_name not in ('-h,--help,-v,--version,-b,-G,-g,-r,-u,-d,-w,-N,-o,-l,-F,--reference-point'):
        message()
    #print("Usage: python tmp1.py <--scale-regions or --reference-point> -b <bamfiles> -G <bedfiles> -g <genefiles> -r <referencePoint {start or end}> -u <upstream_distance> -d <downstream_distance> -N <Equal_block_number> -l <legend_names> -o <outputfiles average.matrix> -F <Figure_name>")
Option={}
for option in opts:
    Option[option[0]]=option[1]

#将不同片段大小的read存入对应长度的列表中,生成{50:[read1,read2,...],51:[read1,read2,...],...,200:[read1,read2,...]}
def SizeClassify(bam,minLength,maxLength):
    SizeDict={}
    for read in bam.fetch():
        size=abs(read.template_length)
        if minLength <= size <= maxLength:
            if size in SizeDict:
                SizeDict[size].append(read)
            else:
                SizeDict[size]=[]
                SizeDict[size].append(read)
    return SizeDict

#ReadInfo {query_name:['chr1',ref_start,ref_end,length],...}
#region {'chr1':{'1045':546,''1067:'453',...},...}
#直接输入一个包含pysam reads对象的列表，就可以返回一个region字典
def bamTodict(Input):
    ReadInfo={}
    for read in Input:
        ID=read.query_name
        if ID in ReadInfo:
            if int(read.reference_start) < int(ReadInfo[ID][1]):
                ReadInfo[ID][1]=int(read.reference_start)
                ReadInfo[ID][3]=int(ReadInfo[ID][2])-int(ReadInfo[ID][1])+1
            else:
                ReadInfo[ID][2]=int(read.reference_end)
                ReadInfo[ID][3]=int(ReadInfo[ID][2])-int(ReadInfo[ID][1])+1
        else:
            a=[]
            a.append(str(read.reference_name))
            a.append(int(read.reference_start))
            a.append(int(read.reference_end))
            a.append(int(a[2])-int(a[1])+1)
            ReadInfo[ID]=a
    region={}
    for value in ReadInfo.values():
        length=value[3]
        chrom=value[0]
        remain=length%2
        if int(remain) == 0:
            centerPoint=int((int(value[1]) + int(value[2]))/2)
            if chrom in region:
                pass
            else:
                region[chrom]={}
            if str(centerPoint) in region[chrom]:
                region[chrom][str(centerPoint)]+=2
            else:
                region[chrom][str(centerPoint)]=2
        else:
            centerPoint1=int((int(value[1]) + int(value[2])+1)/2)
            centerPoint2=int((int(value[1]) + int(value[2])-1)/2)
            if chrom in region:
                pass
            else:
                region[chrom]={}
            if str(centerPoint1) in region[chrom]:
                region[chrom][str(centerPoint1)]+=1
            else:
                region[chrom][str(centerPoint1)]=1
            if str(centerPoint2) in region[chrom]:
                region[chrom][str(centerPoint2)]+=1
            else:
                region[chrom][str(centerPoint2)]=1
    return region
    #bf.close()

#Lis [[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],...]
#regard coordinate start as centerpoint
def bedCount(bedf,Region,count,wdsize,up,down):
    bed6col=open(bedf,"r")
    Size_X_Bin_X_Howmanycoors=[]
    sumReads=0
    #遍历bed文件，生成一个3维数组，每一张网是bed中的一行
    for bed in bed6col:
        row2=bed.rstrip().split()
        Sizes=Region.keys()#把片段大小取出来
        Sizes.sort() #按升序排列
        MutiSize=[]
        #----------------
        if "+" in row2[5]:
            left=int(row2[1])-up
            right=int(row2[1])+down
            #遍历片段大小，每个片段大小生成一行coverage，注意这里需要按顺序来
            for size in Sizes:
                #分bin算完一行之后，把每一行添加到MutiSize中
                OneSize=[]
                for k in range(int(left),int(right),wdsize):
                    for j in range(k,(k+wdsize)):
                        #此处默认了row2[0]染色体已经在字典中了
                        if str(j) in Region[size][row2[0]]:
                            sumReads+=int(Region[size][row2[0]][str(j)])
                        else:
                            pass
                    OneSize.append((int(sumReads)*10000000)/(count*wdsize))
                    sumReads=0
                MutiSize.append(OneSize)
        else:
            left=int(row2[2])-down
            right=int(row2[2])+up
            for size in Sizes:
                #分bin算完一行之后，把每一行添加到MutiSize中
                OneSize=[]
                for k in range(int(right),int(left),-wdsize):
                    for j in range(k,(k-wdsize),-1):
                        if str(j) in Region[size][row2[0]]:
                            sumReads+=int(Region[size][row2[0]][str(j)])
                        else:
                            pass
                    OneSize.append((int(sumReads)*10000000)/(count*wdsize))
                    sumReads=0
                MutiSize.append(OneSize)
        #----------------
        Size_X_Bin_X_Howmanycoors.append(MutiSize)
    output=np.array(Size_X_Bin_X_Howmanycoors)
    return output
    bed6col.close()

#Lis [[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],...]
#regard coordinate end as centerpoint
"""
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
"""

#allvalue [[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8]]
#calculate coordinate start-centered reads
"""
def bamTocountStart(bed,Input,windowSize,Up,Down):
    bamf=pysam.AlignmentFile(Input,"rb")
    region=bamTodict(bamf)
    num=bamf.count()
    allvalue=bedCount(bed,region,num,windowSize,Up,Down)
    return allvalue
"""
#allvalue [[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8]]
#calculate coordinate end-centered reads
"""
def bamTocountEnd(bed,Input,windowSize,Up,Down):
    bamf=pysam.AlignmentFile(Input,"rb")
    region=bamTodict(bamf)
    num=bamf.count()
    allvalue=bedCountEnd(bed,region,num,windowSize,Up,Down)
    return allvalue
"""

if "--reference-point" in Option:
    ylength=int(up/wdsize)
    left=int(0)
    medium=ylength
    right=int((up+down)/wdsize)
    bins=right
    if "-r" in Option and str(referencePoint)=="start":
        out_matrix=open(outf,"w")
        #out_figure=open(figurename,"w")
        #biaotou = []
        #for i in range(0,right+2):
        #    biaotou.append("a")
        #print(*biaotou,sep="\t",end="\n",file=out_matrix)
        bamf=pysam.AlignmentFile(bamfile,"rb")
        total=bamf.count()
        sizedict=SizeClassify(bamf,50,200)
        sizes=sizedict.keys()
        sizes.sort()
        #sizeratio=np.array([[float(len(sizedict[i]))/total for _ in range(bins)] for i in sizes])
        #将不同片段大小的reads分别生成一个dict
        region={}
        for size,rawread in sizedict.items():
            region[size]=bamTodict(rawread)
        values=bedCount(bedfile,region,total,wdsize,up,down)
        print("dimension:",values.shape)
        avgValues=np.mean(values,axis=0)
        #outmat=avgValues/sizeratio
        for line in avgValues:
            print(*line,sep="\t",end="\n",file=out_matrix)

       # plt.legend(legendname)
       # plt.ylabel("Normalized reads count")
       # plt.xticks([left,medium,right],["-"+L+"kb","TSS",R+"kb"])
        out_matrix.close()
"""
    elif "-r" in Option and str(referencePoint)=="end":
        out_profile=open(outf[0],"w")
        out_matrix=open(outf[1],"w")
        out_figure=open(figurename,"w")
        figurelist=[]
        biaotou = []
        for i in range(0,right+2):
            biaotou.append("a")
        print(*biaotou,sep="\t",end="\n",file=out_matrix)
        for Bam in range(len(bamlist)):
            bamf=pysam.AlignmentFile(bamlist[Bam],"rb")
            region=bamTodict(bamf)
            num=bamf.count()
            for Bed in range(len(bedlist)):
                #values=bamTocountEnd(bedlist[Bed],bamlist[Bam],wdsize,up,down)
                values=bedCountEnd(bedlist[Bed],region,num,wdsize,up,down)
                tup1=np.array(values)
                tup2=tup1.mean(axis=0)
                figurelist.append(tup2)
                #tmpnumber = len(values[0])
                for tmp1 in values:
                    print("a","a",*tmp1,sep="\t",end="\n",file=out_matrix)
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
        plt.xticks([left,medium,right],["-"+L+"kb","TTS",R+"kb"])
        plt.savefig(out_figure)     
        out_figure.close()
        out_profile.close()
        out_matrix.close()
"""