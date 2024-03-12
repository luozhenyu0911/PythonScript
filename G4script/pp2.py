from __future__ import print_function
from __future__ import division
#import GenomeMean
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

def bamTodict(bf):
    region={}
    for read in bf.fetch():
        a=[]
        a.append(read.reference_name)
        a.append(read.reference_start)
        a.append(read.reference_end)
        chrom = str(a[0])
        length = int(a[2])-int(a[1])
        remain = length%2
        if int(remain) == 0:
            centerpoint = int((int(a[1]) + int(a[2]))/2)
            if chrom in region:
                pass
            else:
                region[chrom] = {}
            if str(centerpoint) in region[chrom]:
                region[chrom][str(centerpoint)] += 1
            else:
                region[chrom][str(centerpoint)] = 1
        else:
            centerpoint2=int((int(a[1])+int(a[2])+1)/2)
            centerpoint1=int((int(a[1])+int(a[2])-1)/2)
            if chrom in region:
                pass
            else:
                region[chrom] = {}
            if str(centerpoint2) in region[chrom]:
                region[chrom][str(centerpoint2)] += 0.5
            else:
                region[chrom][str(centerpoint2)] = 0.5
            if str(centerpoint1) in region[chrom]:
                region[chrom][str(centerpoint1)] += 0.5
            else:
                region[chrom][str(centerpoint1)] = 0.5
    return region
    bf.close()

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
	    left=int(row2[2])-up
	    right=int(row2[2])+down
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
            left=int(row2[2])-up
            right=int(row2[2])+down
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

def BodyGenomeCount(bedf,readdict,count,wdsize,up,down,BlockNum):
    bed6col=open(bedf,"r")
    sumReads=0
    genomeRead={}
    for bed in bed6col:
        row2=bed.rstrip().split()
        sublist=[]
       # blockSize=((int(row2[2])-int(row2[1]))/BlockNum)+1
      #  reside=BlockNum-((int(row2[2])-int(row2[1]))%BlockNum)
        reside=BlockNum-int(((int(row2[2])-int(row2[1]))%BlockNum))
        blockSize=int(((int(row2[2])-int(row2[1]))+reside)/BlockNum)
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
            for x in range(int(row2[1]),int(row2[2])+reside,blockSize):
                for y in range(x,(x+blockSize)):
                    if str(y) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(y)])
                    else:
                        pass
              #  sumReads=int((sumReads/blockSize)*10)
                sublist.append((int(sumReads)*1000000)/(count*blockSize))
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
            for x in range(int(row2[2]),int(row2[1])-reside,-blockSize):
                for y in range(x,(x-blockSize),-1):
                    if str(y) in readdict[row2[0]]:
                        sumReads+=int(readdict[row2[0]][str(y)])
                    else:
                        pass
           #     sumReads=int((sumReads/blockSize)*10)
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

def BodyCount(bedf,readdict,count,wdsize,up,down,BlockNum):
    bed6col=open(bedf,"r")
    Lis=[]
    sumReads=0
    for bed in bed6col:
	row2=bed.rstrip().split()
	sublist=[]
       # blockSize=((int(row2[2])-int(row2[1]))/BlockNum)+1
       # reside=BlockNum-((int(row2[2])-int(row2[1]))%BlockNum)
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

if "-c" in Option:

    if "--reference-point" in Option:
        Ylength=int(up+down)
        ylength=int(Ylength/wdsize)
        left=int(0)
        medium=int(ylength/2)
        right=int(ylength)
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
        Ylength=int(up+down)
        ylength=int(Ylength/wdsize)
        left=int(0)
        medium=int(ylength/2)
        right=int(ylength)
        if "-r" in Option and str(referencePoint)=="start":
            if "-g" not in Option:
                out_profile=open(outf[0],"w")
                out_matrix=open(outf[1],"w")
                out_figure=open(figurename,"w")
                figurelist=[]
                for Bam in range(len(bamlist)):
                    for Bed in range(len(bedlist)):
                        values=bamTocountStart(bedlist[Bed],bamlist[Bam],wdsize,up,down)
                        tup1=np.array(values)	
                        tup2=tup1.mean(axis=0)
                        figurelist.append(tup2)
                        tmpnumber = len(values[0])
                        biaotou = []
                        for i in range(0,tmpnumber+2):
                            biaotou.append("a")
                        print(*biaotou,sep="\t",end="\n",file=out_matrix)
                        for tmp1 in values:
                            print("a","a",*tmp1,sep="\t",end="\n",file=out_matrix)
                  #      for value in tup2:
                   #         print(value,file=out_profile)	
                   #     for tmp1 in values:
                    #        for tmp2 in tmp1:
                     #           print(tmp2,end="\t",file=out_matrix)
                      #      print(file=out_matrix)
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
                            finalmat=Listgene(genelist[Gene],values)
                            tup1=np.array(finalmat)
                            tup2=tup1.mean(axis=0)
                            figurelist.append(tup2)
                            tmpnumber = len(values[0])
                            biaotou = []
                            for i in range(0,tmpnumber+2):
                                biaotou.append("a")
                            print(*biaotou,sep="\t",end="\n",file=out_matrix)
                            for tmp1 in values:
                                print("a","a",*tmp1,sep="\t",end="\n",file=out_matrix)
                  #          for value in tup2:
                   #             print(value,file=out_profile)
                   #         for tmp1 in finalmat:			
                    #            for tmp2 in tmp1:
                     #               print(tmp2,end="\t",file=out_matrix)
                      #          print(file=out_matrix)
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
                        tup1=np.array(values)
                        tup2=tup1.mean(axis=0)
                        figurelist.append(tup2)
                        tmpnumber = len(values[0])
                        biaotou = []
                        for i in range(0,tmpnumber+2):
                            biaotou.append("a")
                        print(*biaotou,sep="\t",end="\n",file=out_matrix)
                        for tmp1 in values:
                            print("a","a",*tmp1,sep="\t",end="\n",file=out_matrix)
                  #      for value in tup2:
                   #         print(value,file=out_profile)
                 #       for tmp1 in values:
                  #          for tmp2 in tmp1:
                   #             print(tmp2,end="\t",file=out_matrix)
                    #        print(file=out_matrix)
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
                for Bam in range(len(bamlist)):
                    for Bed in range(len(bedlist)):
                        for Gene in range(len(genelist)):
                            values=bamtocountEnd(bedlist[Bed],bamlist[Bam],wdsize,up,down)
                            finalmat=Listgene(genelist[Gene],values)
                            tup1=np.array(finalmat)
                            tup2=tup1.mean(axis=0)
                            figurelist.append(tup2)
                            tmpnumber = len(values[0])
                            biaotou = []
                            for i in range(0,tmpnumber+2):
                                biaotou.append("a")
                            print(*biaotou,sep="\t",end="\n",file=out_matrix)
                            for tmp1 in values:
                                print("a","a",*tmp1,sep="\t",end="\n",file=out_matrix)
                  #          for value in tup2:
                   #             print(value,file=out_profile)
                     #       for tmp1 in finalmat:
                      #          for tmp2 in tmp1:
                       #             print(tmp2,end="\t",file=out_matrix)
                        #        print(file=out_matrix)
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
        right=int(down/wdsize)+mediumA+blocknum
        if "-g" not in Option:
            out_profile=open(outf[0],"w")
            out_matrix=open(outf[1],"w")
            out_figure=open(figurename,"w")
            figurelist=[]
            for Bam in range(len(bamlist)):
                for Bed in range(len(bedlist)):
                    values=bamtobodycount(bedlist[Bed],bamlist[Bam],wdsize,up,down,blocknum)
                    tup1=np.array(values)
                    tup2=tup1.mean(axis=0)
                    figurelist.append(tup2)
                    tmpnumber = len(values[0])
                    biaotou = []
                    for i in range(0,tmpnumber+2):
                        biaotou.append("a")
                    print(*biaotou,sep="\t",end="\n",file=out_matrix)
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
            out_profile=open(outf[0],"w")
            out_matrix=open(outf[1],"w")
            out_figure=open(figurename,"w")
            figurelist=[]
            for Bam in range(len(bamlist)):
                for Bed in range(len(bedlist)):
                    for Gene in range(len(genelist)):
                        values=bamGenometobodycount(bedlist[Bed],bamlist[Bam],wdsize,up,down,blocknum)
                        finalmat=Listgene(genelist[Gene],values)
                        tup1=np.array(finalmat)
                        tup2=tup1.mean(axis=0)
                        figurelist.append(tup2)
                        tmpnumber = len(values[0])
                        biaotou = []
                        for i in range(0,tmpnumber+2):
                            biaotou.append("a")
                        print(*biaotou,sep="\t",end="\n",file=out_matrix)
                        for tmp1 in values:
                            print("a","a",*tmp1,sep="\t",end="\n",file=out_matrix)

                       # for value in tup2:
                        #    print(value,file=out_profile)
                   #     for tmp1 in finalmat:
                    #        for tmp2 in tmp1:
                     #           print(tmp2,end="\t",file=out_matrix)
                      #      print(file=out_matrix)            
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
