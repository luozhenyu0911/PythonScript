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
    print("\n"+"Usage: python **.py --reference-point or --scale-regions [options]"+"\n")
    print("Description:"+"\n"+"\t"+"\t"+"This tool creates scores and profile plot over sets of genomic regions. Typically, these regions are genes, but any other regions defined in BED will work."+"\n")
    print("Options:"+"\n")
    print("\t"+"-v, "+"--version"+"\t"+"show program's version number and exit")
    print("\t"+"-h, "+"--help"+"\t"+"show this help message and exit.")
    print("\t"+"--reference-point or --scale-regions"+"\t"+"Reference-point refers to a position within a BED region(e.g., the starting point). In this mode, only those genomicpositions before (upstream) and/or after (downstream) of the reference point will be plotted."+"In the scale-regions mode, all regions in the BED file arestretched or shrunken to the length (in bases) indicated by the user.")
    print("\t"+"-c"+"    <Input bam files>"+"\t"+"Control data files(Support for multiple files, Separated by commas, the order is consistent with bam files), should be sorted bam format.")
    print("\t"+"-b"+"    <bamfiles>"+"\t"+"data files(Support for multiple files, Separated by commas), should be sorted bam format.")
    print("\t"+"-G"+"    <bedfiles>"+"\t"+"data files(Support for multiple files, Separated by commas), should be original bed files.")
    print("\t"+"-g"+"    <gene list files>"+"\t"+"data files(Support for multiple files, Separated by commas), should be a row genes for every file.")
    print("\t"+"-r"+"    <referencePoint {start or end}>"+"\t"+"The reference point for the plotting could be either the region start (TSS), the region end (TES) of the region.")
    print("\t"+"-u"+"    <upstream distance>"+"\t"+"Distance upstream of the reference-point selected.")
    print("\t"+"-d"+"    <downstream distance>"+"\t"+"Distance downstream of the reference-point selected.")
    print("\t"+"-w"+"    <window size>"+"\t"+"Length, in bases, binSize for averaging the score over the regions length.")
    print("\t"+"-N"+"    <Equally fragment number>"+"\t"+"Equally divided into the region into several fragments.")
    print("\t"+"-l"+"    <legend names>"+"\t"+"Name of the plot legend. (Support for multiple legends)")
    print("\t"+"-o"+"    <outputfiles average data,final matrix>"+"\t"+"File name to save the average data file and final matrix file. (Separated by commas)")
    print("\t"+"-F"+"    <Figure name>"+"\t"+"Name of the plot name. (plot file format is PDF)"+"\n")
    print("An example usage is:"+"\n")
    print("\t"+"python **.py --reference-point -b **.bam,**.bam... -G **.bed,**.bed... -g **.txt,**.txt... -r start/end -u 1000 -d 1000 -w 50 -l **,** -o **_dens,**_matrix -F **.pdf"+"\n")
    print("\t"+"python **.py --scale-regions -b **.bam,**bam... -G **.bed,**.bed... -g **.txt,**.txt... -u 1000 -d 1000 -w 50 -N 100 -l **,** -o **_dens,**_matrix -F **.pdf"+"\n")
    return

opts,args = getopt.getopt(sys.argv[1:],'-h-v-b:-G:-g:-u:-r:-d:-w:-N:-o:-l:-F:-c:',['help','version','scale-regions','reference-point'])
for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
	message()
    if opt_name in ('-v','--version'):
	print("**.py 1.0 2018-07-22")
    if opt_name in ('-b'):
        bamfile = opt_value
	bamlist=bamfile.split(",")
    if opt_name in ('-G'):
        bedfile = opt_value
	bedlist=bedfile.split(",")
    if opt_name in ('-g'):
        genefile = opt_value
	genelist=genefile.split(",")
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
    if opt_name in ('-N'):
        Blocknum = opt_value
	blocknum=int(Blocknum)
    if opt_name in ('-F'):
        figurename = opt_value
    if opt_name in ('-c'):
        Inputbam = opt_value
        Inputbamlist=Inputbam.split(",")
    if opt_name in ('-l'):
	legend = opt_value
	legendname=legend.split(",")
    if opt_name in ('-o'):
        outputfile = opt_value
	outf=outputfile.split(",")
    if opt_name not in ('-h,--help,-v,--version,-b,-G,-g,-r,-u,-d,-w,-N,-o,-l,-c,-F,--scale-regions,--reference-point'):
	message()
	#print("Usage: python tmp1.py <--scale-regions or --reference-point> -b <bamfiles> -G <bedfiles> -g <genefiles> -r <referencePoint {start or end}> -u <upstream_distance> -d <downstream_distance> -N <Equal_block_number> -l <legend_names> -o <outputfiles average.matrix> -F <Figure_name>")
Option={}
for option in opts:
    Option[option[0]]=option[1]

#ReadInfo {query_name:['chr1',ref_start,ref_end,length],...}
#region {'chr1':{'1045':546,''1067:'453',...},...}
def bamTodict(Input):
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
    bf.close()

#readdict {'chr1':{'1045':546,''1067:'453',...},...}
#sublist [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8]
#genomeRead {'Loc****':[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],'Loc':...}
#regard coordinate start as centerpoint
def GenomebedCount(bedf,readdict,count,wdsize,up,down):
    bed6col=open(bedf,"r")
    sumReads=0
    genomeRead={}
    for bed in bed6col:
        row2=bed.rstrip().split()
        sublist=[]
        if "+" in row2[5]:
            left=int(row2[1])-up
            right=int(row2[1])+down
            for k in range(int(left),int(right),wdsize):
                for j in range(k,(k+wdsize)):
                    if str(j) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(j)])
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
                    if str(j) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(j)])
                    else:
                        pass
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            genomeRead[row2[3]]=sublist
    return genomeRead
    bed6col.close()

#genomeRead {'Loc****':[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],'Loc':...}
#regard coordinate end as centerpoint
def GenomebedCountEnd(bedf,readdict,count,wdsize,up,down):
    bed6col=open(bedf,"r")
    sumReads=0
    genomeRead={}
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
            genomeRead[row2[3]]=sublist
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
            genomeRead[row2[3]]=sublist
    return genomeRead
    bed6col.close()

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

    
    
def equalDivide(x,y,part):
    width=int(y-x)
    L = []
    if width % part == 0:
        for _ in range(part):
            L.append(int(width/part))
    else:
        reside = width%part
        main = int(width/part)
        for _ in range(part):
            if reside != 0:
                L.append(main+1)
                reside += -1
            else:
                L.append(main)
    J=[]
    start=x
    for i in L:
        J.append([start,start+i])
        start+=i
    return(J)

#genomeRead {'Loc****':[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],'Loc':...}
# up TSS TES down
def BodyGenomeCount(bedf,readdict,count,wdsize,up,down,BlockNum):
    bed6col=open(bedf,"r")
    sumReads=0
    genomeRead={}
    for bed in bed6col:
        row2=bed.rstrip().split()
        sublist=[]
        block = equalDivide(int(row2[1]),int(row2[2]),BlockNum)
        #reside=BlockNum-int(((int(row2[2])-int(row2[1]))%BlockNum))
        #blockSize=int(((int(row2[2])-int(row2[1]))+reside)/BlockNum)
        if "+" in row2[5]:
            for k in range(int(row2[1])-up,int(row2[1]),wdsize):
                for j in range(k,(k+wdsize)):
                    if str(j) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(j)])
                    else:
                        pass
              #  sumReads=int((sumReads/10)*10)
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            #for x in range(int(row2[1]),int(row2[2])+reside,blockSize):
            #    for y in range(x,(x+blockSize)):
            for x in block:
                for y in range(x[0],x[1]):
                    if str(y) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(y)])
                    else:
                        pass
              #  sumReads=int((sumReads/blockSize)*10)
                blockSize=int(x[1]-x[0])
                sublist.append((int(sumReads)*10000000)/(count*blockSize))
                sumReads=0
            for u in range(int(row2[2]),int(row2[2])+down,wdsize):
                for v in range(u,(u+wdsize)):
                    if str(v) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(v)])
                    else:
                        pass
             #   sumReads=int((sumReads/10)*10)
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            genomeRead[row2[3]]=sublist
        else:
            for k in range(int(row2[2])+up,int(row2[2]),-wdsize):
                for j in range(k,(k-wdsize),-1):
                    if str(j) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(j)])
                    else:
                        pass
            #    sumReads=int((sumReads/10)*10)
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            for x in reversed(block):
                for y in range(x[1],x[0],-1):
                    if str(y) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(y)])
                    else:
                        pass
           #     sumReads=int((sumReads/blockSize)*10)
                blockSize=int(x[1]-x[0])
                sublist.append((int(sumReads)*10000000)/(count*blockSize))
                sumReads=0
            for u in range(int(row2[1]),int(row2[1])-down,-wdsize):
                for v in range(u,(u-wdsize),-1):
                    if str(v) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(v)])
                    else:
                        pass
          #      sumReads=int((sumReads/10)*10)
                sublist.append((int(sumReads)*10000000)/(count*wdsize))
                sumReads=0
            genomeRead[row2[3]]=sublist
    return genomeRead
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
        block = equalDivide(int(row2[1]),int(row2[2]),BlockNum)
        #reside=BlockNum-int(((int(row2[2])-int(row2[1]))%BlockNum))
        #blockSize=int(((int(row2[2])-int(row2[1]))+reside)/BlockNum)
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
            for x in block:
                for y in range(x[0],x[1]):
                    if str(y) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(y)])
                    else:
                        pass
              #  sumReads=int((sumReads/blockSize)*10)
                blockSize=int(x[1]-x[0])
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
            for x in reversed(block):
                for y in range(x[1],x[0],-1):
                    if str(y) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(y)])
                    else:
                        pass
           #     sumReads=int((sumReads/blockSize)*10)
                blockSize=int(x[1]-x[0])
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

#{'Loc****':[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],'Loc':...}
#calculate gene TSS-centered reads
def bamtocountStart(bed,Input,windowSize,Up,Down):
    bamf=pysam.AlignmentFile(Input,"rb")
    region=bamTodict(bamf)
    num=bamf.count()
    allgene=GenomebedCount(bed,region,num,windowSize,Up,Down)
    return allgene

#{'Loc****':[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],'Loc':...}
#calculate gene TES-centered reads
def bamtocountEnd(bed,Input,windowSize,Up,Down):
    bamf=pysam.AlignmentFile(Input,"rb")
    region=bamTodict(bamf)
    num=bamf.count()
    allgene=GenomebedCountEnd(bed,region,num,windowSize,Up,Down)
    return allgene

#allvalue [[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8]]
#calculate coordinate start-centered reads
def bamTocountStart(bed,Input,windowSize,Up,Down):
    bamf=pysam.AlignmentFile(Input,"rb")
    region=bamTodict(bamf)
    num=bamf.count()
    allvalue=bedCount(bed,region,num,windowSize,Up,Down)
    return allvalue

#allvalue [[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8]]
#calculate coordinate end-centered reads
def bamTocountEnd(bed,Input,windowSize,Up,Down):
    bamf=pysam.AlignmentFile(Input,"rb")
    region=bamTodict(bamf)
    num=bamf.count()
    allvalue=bedCountEnd(bed,region,num,windowSize,Up,Down)
    return allvalue

#bodyvalue [[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8]]
#calculate up start end down reads
def bamtobodycount(bed,Input,windowSize,Up,Down,BlockNUM):
    bamf=pysam.AlignmentFile(Input,"rb")
    region=bamTodict(bamf)
    num=bamf.count()
    bodyvalue=BodyCount(bed,region,num,windowSize,Up,Down,BlockNUM)
    return bodyvalue

#Genomebodyvalue {'Loc****':[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.8],'Loc':...}
#calculate genebody reads
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

if "-c" in Option:

    if "--reference-point" in Option:
        ylength=int(up/wdsize)
        left=int(0)
        medium=ylength
        right=int((up+down)/wdsize)
        if "-r" in Option and str(referencePoint)=="start":
            if "-g" not in Option:
                out_profile=open(outf[0],"w")
                out_matrix=open(outf[1],"w")
                out_figure=open(figurename,"w")
                figurelist=[]
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
                 #       for value in tup2:
                  #          print(value,file=out_profile)			
                        for tmp1 in Values:
                            for tmp2 in tmp1:
                                print(tmp2,end="\t",file=out_matrix)
                            print(file=out_matrix)
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
                          #  for value in tup2:
                           #     print(value,file=out_profile)
                            for tmp1 in Values:			
                                for tmp2 in tmp1:
                                    print(tmp2,end="\t",file=out_matrix)
                                print(file=out_matrix)
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
        elif "-r" in Option and str(referencePoint)=="end":
            if "-g" not in Option:
                out_profile=open(outf[0],"w")
                out_matrix=open(outf[1],"w")
                out_figure=open(figurename,"w")
                figurelist=[]
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
                 #       for value in tup2:
                  #          print(value,file=out_profile)
                        for tmp1 in Values:
                            for tmp2 in tmp1:
                                print(tmp2,end="\t",file=out_matrix)
                            print(file=out_matrix)
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
                   #         for value in tup2:
                    #            print(value,file=out_profile)
                            for tmp1 in Values:
                                for tmp2 in tmp1:
                                    print(tmp2,end="\t",file=out_matrix)
                                print(file=out_matrix)
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

    elif "--scale-regions" in Option:
        left=int(0)
        mediumA=int(up/wdsize)
        mediumB=mediumA+blocknum
        right=int(down/wdsize)+mediumA+blocknum
        if "-g" not in Option:
            out_profile=open(outf[0],"w")
            out_matrix=open(outf[1],"w")
            out_figure=open(figurename,"w")
            figurelist=[]
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
             #       for value in tup2:
              #          print(value,file=out_profile)
                    for tmp1 in Values:
                        for tmp2 in tmp1:
                            print(tmp2,end="\t",file=out_matrix)
                        print(file=out_matrix)
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
                 #       for value in tup2:
                  #          print(value,file=out_profile)
                        for tmp1 in Values:
                            for tmp2 in tmp1:
                                print(tmp2,end="\t",file=out_matrix)
                        print(file=out_matrix)
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

    if "--reference-point" in Option:
        ylength=int(up/wdsize)
        left=int(0)
        medium=ylength
        right=int((up+down)/wdsize)
        if "-r" in Option and str(referencePoint)=="start":
            if "-g" not in Option:
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
                        #values=bamTocountStart(bedlist[Bed],bamlist[Bam],wdsize,up,down)
                        values=bedCount(bedlist[Bed],region,num,wdsize,up,down)
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
                biaotou = []
                for i in range(0,right+2):
                    biaotou.append("a")
                print(*biaotou,sep="\t",end="\n",file=out_matrix)
                for Bam in range(len(bamlist)):
                    for Bed in range(len(bedlist)):
                        for Gene in range(len(genelist)):
                            values=bamtocountStart(bedlist[Bed],bamlist[Bam],wdsize,up,down)
                            finalmat=Listgene(genelist[Gene],values)
                            tup1=np.array(finalmat)
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
                plt.xticks([left,medium,right],["-"+L+"kb","TSS",R+"kb"])
                plt.savefig(out_figure)     
                out_figure.close()
                out_profile.close()
                out_matrix.close()
        elif "-r" in Option and str(referencePoint)=="end":
            if "-g" not in Option:
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
            else:
                out_profile=open(outf[0],"w")
                out_matrix=open(outf[1],"w")
                out_figure=open(figurename,"w")
                figurelist=[]
                biaotou = []
                for i in range(0,tmpnumber+2):
                    biaotou.append("a")
                print(*biaotou,sep="\t",end="\n",file=out_matrix)
                for Bam in range(len(bamlist)):
                    for Bed in range(len(bedlist)):
                        for Gene in range(len(genelist)):
                            values=bamtocountEnd(bedlist[Bed],bamlist[Bam],wdsize,up,down)
                            finalmat=Listgene(genelist[Gene],values)
                            tup1=np.array(finalmat)
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
        else:
            pass
#print("Usage: python tmp1.py <<--scale-regions or --reference-point>> -b <bamfile> -G <bedfile> -g <genelist> -r <referencePoint {start or end}> -u <upstream_distance> -d <downstream_distance> -w <window_size> -bn <Equal_block_number> -o <outputfile>")

    elif "--scale-regions" in Option:
        left=int(0)
        mediumA=int(up/wdsize)
        mediumB=mediumA+blocknum
        right=int(up/wdsize)+blocknum+int(down/wdsize)
        if "-g" not in Option:
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
                    #values=bamtobodycount(bedlist[Bed],bamlist[Bam],wdsize,up,down,blocknum)
                    values=BodyCount(bedlist[Bed],region,num,wdsize,up,down,blocknum)
                    tup1=np.array(values)
                    tup2=tup1.mean(axis=0)
                    figurelist.append(tup2)
                    #tmpnumber = len(values[0])
                    if len(bamlist) > 1:
                        for tmp1 in values:
                            print(legendname[Bam],"a",*tmp1,sep="\t",end="\n",file=out_matrix)
                    if len(bedlist) > 1:
                        for tmp1 in values:
                            print(legendname[Bed],"a",*tmp1,sep="\t",end="\n",file=out_matrix)
                    if len(bamlist) == len(bedlist) == 1:
                        for tmp1 in values:
                            print(legendname[0],"a",*tmp1,sep="\t",end="\n",file=out_matrix)
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
            plt.xticks([left,mediumA,mediumB,right],["-"+L+"kb","TSS","TTS",R+"kb"])
            plt.savefig(out_figure)
            out_figure.close()
            out_profile.close()
            out_matrix.close()
        else:
            out_profile=open(outf[0],"w")
            out_matrix=open(outf[1],"w")
            out_figure=open(figurename,"w")
            figurelist=[]
            biaotou = []
            for i in range(0,right+2):
                biaotou.append("a")
            print(*biaotou,sep="\t",end="\n",file=out_matrix)
            for Bam in range(len(bamlist)):
                for Bed in range(len(bedlist)):
                    for Gene in range(len(genelist)):
                        values=bamGenometobodycount(bedlist[Bed],bamlist[Bam],wdsize,up,down,blocknum)
                        finalmat=Listgene(genelist[Gene],values)
                        tup1=np.array(finalmat)
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
            plt.xticks([left,mediumA,mediumB,right],["-"+L+"kb","TSS","TTS",R+"kb"])
            plt.savefig(out_figure)     
            out_figure.close()
            out_profile.close()
            out_matrix.close()
    else:
        pass
