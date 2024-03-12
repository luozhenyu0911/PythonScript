from __future__ import print_function
from __future__ import division
from itertools import combinations
import sys
import getopt
#import gzip
import numpy as np
import scipy
import scipy.stats
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def message():
    print("\n"+"Usage: python reads_density.py --reference-point or --scale-regions [options]"+"\n")
    print("Description:"+"\n"+"\t"+"\t"+"This tool creates scores and profile plot over sets of genomic regions. Typically, these regions are genes, but any other regions defined in BED will work."+"\n")
    print("Options:"+"\n")
    print("\t"+"-h, "+"--help"+"\t"+"show this help message and exit.")
    print("\t"+"-b"+"    <CGmap>"+"\t"+"data files, should be gzip format.")
    print("\t"+"-B"+"    <list file>"+"\t"+"multiple input data files list in a text file, one CGmap file in a line, should be gzip.")
    print("\t"+"-c"+"    <int>"+"\t"+"minimun coverage for C sites.")
    print("\t"+"-S"+"    <string>"+"\t"+"strands orientation could be forward reverse or all.")
    print("\t"+"-r"+"    <bedfiles>"+"\t"+"data files, should be original bed files with 6 columns.")
    print("\t"+"-g"+"    <gene list files>"+"\t"+"data file must end with .list, should be a gene ID each row.")
    print("\t"+"-G"+"    <list file>"+"\t"+"multiple input data files list in a text file, each data file must end with .list, should be a gene ID each row.")
    print("\t"+"-m"+"    <int>"+"\t"+"The number of blocks will genebody be divided.")
    print("\t"+"-u"+"    <upstream distance>"+"\t"+"Distance upstream of the reference-point selected.")
    print("\t"+"-d"+"    <downstream distance>"+"\t"+"Distance downstream of the reference-point selected.")
    print("\t"+"-w"+"    <window size>"+"\t"+"Length, in bases, binSize for averaging the score over the regions length.")
    print("\t"+"-f"+"    <fpkm file>"+"\t"+"Text file contains two columns, GeneIDs and expression values.")
    print("\t"+"-l"+"    <float num>"+"\t"+"Threshold for not express group, genes with a FPKM lower than this value were considered as not expressed ")
    print("\t"+"-n"+"    <int num>"+"\t"+"number of groups divided according to the expression levels of all expressed genes")
    print("\t"+"-o"+"    <outputprefix or ouputfilename>"+"\t"+"File name to save the average data file. If mutiple files were submitted, the output file will use this string as prefix.")
    print("\t"+"-F"+"    <Figure name>"+"\t"+"Name of the plot name. (plot file format is PDF)"+"\n")
    #print("\t"+"--grouprev"+"\t"+"If mutiple bamfiles and gene lists were subject, use this option to plot different lists of different genes for each bamfiles"+"\n")
    print("An example usage is:"+"\n")
    print("\t"+"python reads_density.py -b **.bam -r **.bed -g **.list -m TSS -u 1000 -d 1000 -w 50 -o **.txt -F **.pdf"+"\n")
    print("\t"+"python reads_density.py -b **.bam -r **.bed -g **.list -u 1000 -d 1000 -w 50 -l 0.5 -f **.fpkm -n 5 -o **.txt -F **.pdf"+"\n")
    return

opts,args = getopt.getopt(sys.argv[1:],'-h-b:-r:-G:-g:-u:-d:-o:-w:-B:-f:-n:-l:-F:-c:-C:-m:-S:',['help','heatmap='])
for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        message()
    if opt_name in ('-b'):
        bamIn = opt_value
    if opt_name in ('-r'):
        bedIn = opt_value
    if opt_name in ('-G'):
        mgenelist = opt_value
    if opt_name in ('-g'):
        genelist = opt_value
    if opt_name in ('-u'):
        upstream = opt_value
    if opt_name in ('-d'):
        downstream = opt_value
    if opt_name in ('-o'):
        outputfile = opt_value
    if opt_name in ('-w'):
        wdsize = opt_value
    if opt_name in ('-B'):
        bamIns = opt_value
    if opt_name in ('-f'):
        fpkmgroup = opt_value
    if opt_name in ('-n'):
        groupnum = opt_value
    if opt_name in ('-l'):
        minFpkm = opt_value
    if opt_name in ('-F'):
        figuredir = opt_value 
    if opt_name in ('-c'):
        ccov = opt_value
    if opt_name in ('-m'):
        Bnum = int(opt_value)
    if opt_name in ('--heatmap'):
        heatfile = opt_value
    if opt_name in ('-S'):
        strands = str(opt_value)
    #if opt_name not in ('-G,-w,-o'):
        #print("Usage: python tmp1.py -b <bamfile> -G <bedfile> -g <genefile> -u <upstream_distance> -d <downstream_distance> -o <outputfile>")

figurefilepath=figuredir.split()
if "pdf" in figurefilepath[-1]:
    plt.switch_backend('PDF')
else:
    plt.switch_backend('agg')

allopt=[]
for option in opts:
    allopt.append(option[0])

step_f=int(wdsize)

def bamInfo(Input,Strand):
    region={}
    for read in Input:
        col=read.rstrip().split()
        if "forward" in Strand:
            if int(col[7]) >= int(ccov) and "+" in col[2]: 
                if str(col[0]) in region:
                    region[str(col[0])][str(col[1])]=(int(col[6]),int(col[7]))
                else:
                    region[str(col[0])]={}
                    region[str(col[0])][str(col[1])]=(int(col[6]),int(col[7]))
        elif "reverse" in Strand:
            if int(col[7]) >= int(ccov) and "-" in col[2]:
                if str(col[0]) in region:
                    region[str(col[0])][str(col[1])]=(int(col[6]),int(col[7]))
                else:
                    region[str(col[0])]={}
                    region[str(col[0])][str(col[1])]=(int(col[6]),int(col[7]))
        elif "all" in Strand:
            if int(col[7]) >= int(ccov):
                if str(col[0]) in region:
                    region[str(col[0])][str(col[1])]=(int(col[6]),int(col[7]))
                else:
                    region[str(col[0])]={}
                    region[str(col[0])][str(col[1])]=(int(col[6]),int(col[7]))
    return region


def Readscount(bedlist,readdict,windows,BlockNum):
    metC=0
    sumC=0
    allgene={}
    revwins=0-windows
   # reside=BlockNum-int(((int(row2[3])-int(row2[2]))%BlockNum))
   # blockSize=int(((int(row2[3])-int(row2[2]))+reside)/BlockNum)
    for row2 in bedlist:
        reside=BlockNum-int(((int(row2[3])-int(row2[2]))%BlockNum))
        blockSize=int(((int(row2[3])-int(row2[2]))+reside)/BlockNum)
        values=[]
        if str(row2[6]) is "+":
            if row2[0] in readdict:
                for k in range(int(row2[1]),int(row2[2]),windows):
                    for j in range(k,(k+windows)):
                        if str(j) in readdict[row2[0]]:
                            sumC+=readdict[row2[0]][str(j)][1]
                            metC+=readdict[row2[0]][str(j)][0]
                        else:
                            pass
                    if sumC == 0:
                        values.append(np.nan)
                    else:
                        methyratio=metC/sumC
                        values.append(methyratio)
                    metC=0
                    sumC=0
            else:
                #for j in range(int(row2[1]),int(row2[2]),windows):
                 #   values.append("0")
                pass
            if row2[0] in readdict:
              #  reside=BlockNum-int(((int(row2[3])-int(row2[2]))%BlockNum))
              #  blockSize=int(((int(row2[3])-int(row2[2]))+reside)/BlockNum)
                for x in range(int(row2[2]),int(row2[3])+reside,blockSize):
                    for y in range(x,(x+blockSize)):
                        if str(y) in readdict[row2[0]]:
                            sumC+=readdict[row2[0]][str(y)][1]
                            metC+=readdict[row2[0]][str(y)][0]
                        else:
                            pass
                    if sumC == 0:
                        values.append(np.nan)
                    else:
                        methyratio=metC/sumC
                        values.append(methyratio)
                    metC=0
                    sumC=0
            else:
#                for x in range(int(row2[2]),int(row2[3])+reside,blockSize):
 #                   values.append("0")
                pass
            if row2[0] in readdict:
                for k in range(int(row2[3]),int(row2[4]),windows):
                    for j in range(k,(k+windows)):
                        if str(j) in readdict[row2[0]]:
                            sumC+=readdict[row2[0]][str(j)][1]
                            metC+=readdict[row2[0]][str(j)][0]
                        else:
                            pass
                    if sumC == 0:
                        values.append(np.nan)
                    else:
                        methyratio=metC/sumC
                        values.append(methyratio)
                    metC=0
                    sumC=0
            else:
              #  for j in range(int(row2[3]),int(row2[4]),windows):
               #     values.append("0")
                pass
        else:
            if row2[0] in readdict:
                for k in range(int(row2[4]),int(row2[3]),revwins):
                    for j in range(k,(k+revwins),-1):
                        if str(j) in readdict[row2[0]]:
                            sumC+=readdict[row2[0]][str(j)][1]
                            metC+=readdict[row2[0]][str(j)][0]
                        else:
                            pass
                    if sumC == 0:
                        values.append(np.nan)
                    else:
                        methyratio=metC/sumC
                        values.append(methyratio)
                    metC=0
                    sumC=0
            else:
            #    for j in range(int(row2[4]),int(row2[3]),revwins):
             #       values.append("0")
                pass
            if row2[0] in readdict:
               # reside=BlockNum-int(((int(row2[3])-int(row2[2]))%BlockNum))
               # blockSize=int(((int(row2[3])-int(row2[2]))+reside)/BlockNum)
                for x in range(int(row2[3]),int(row2[2])-reside,-blockSize):
                    for y in range(x,(x-blockSize),-1):
                        if str(y) in readdict[row2[0]]:
                            sumC+=readdict[row2[0]][str(y)][1]
                            metC+=readdict[row2[0]][str(y)][0]
                        else:
                            pass
                    if sumC == 0:
                        values.append(np.nan)
                    else:
                        methyratio=metC/sumC
                        values.append(methyratio)
                    metC=0
                    sumC=0
            else:
              #  for x in range(int(row2[3]),int(row2[2])-reside,-blockSize):
               #     values.append("0")
                pass
            if row2[0] in readdict:
                for k in range(int(row2[2]),int(row2[1]),revwins):
                    for j in range(k,(k+revwins),-1):
                        if str(j) in readdict[row2[0]]:
                            sumC+=readdict[row2[0]][str(j)][1]
                            metC+=readdict[row2[0]][str(j)][0]
                        else:
                            pass
                    if sumC == 0:
                        values.append(np.nan)
                    else:
                        methyratio=metC/sumC
                        values.append(methyratio)
                    metC=0
                    sumC=0
            else:
              #  for j in range(int(row2[2]),int(row2[1]),revwins):
               #     values.append("0")
                pass
        values=map(float,values)
        allgene[str(row2[5])]=values
    return allgene
    bedfile.close()
 
def Listgene(geneID,genedict):
    Input=open(geneID,"r")
    Listmat=[]
    for line in Input:
        Id=line.rstrip()
        Listmat.append(genedict[str(Id)])
    return Listmat

def listg(g,gdict):
    Listmat=[]
    for ID in g:
        Listmat.append(gdict[str(ID)])
    return Listmat


def bamtocount(bed,Input,windowsize,strds,blockNum):
    bamfile=open(Input,"rb")
    region=bamInfo(bamfile,strds)
    geneall=Readscount(bed,region,windowsize,blockNum)
    return geneall

def splist(l,n):
    lenth = len(l)
    sz = lenth // n
    c= lenth % n
    lst = []
    i = 0
    while i < n:
        if i < c:
            bs=sz+1
            lst.append(l[i*bs : i*bs+bs])
        else:
            lst.append(l[i*sz+c : i*sz + c + sz])
        i+=1
    return lst

def plotsmooth(values):
    n=len(values)
    num=range(0,n)
    num=np.array(num)
    xnew=np.linspace(num.min(),num.max(),5*n)
    func=interp1d(num,values,kind="cubic")
    ynew=func(xnew)
    plt.plot(xnew,ynew)
    
bedfile=open(bedIn,"r")
bedpos=[]
m1lab="Start"
m2lab="End"
slab="-"+str(int(upstream)/1000)+"kb"
elab=str(int(downstream)/1000)+"kb"
labs=[slab,m1lab,m2lab,elab]
labp=[0,int(upstream)/step_f,int(upstream)/step_f+Bnum,int(upstream)/step_f+Bnum+int(downstream)/step_f]
for genepos in bedfile:
    row=genepos.rstrip().split()
    [chrom,start,end,geneID,score,direct]=row
    startnew=int(start)-int(upstream)
    endnew=int(end)+int(downstream)
    bedpos.append([chrom,startnew,start,end,endnew,geneID,direct])

if "-b" in allopt:
    out=open(outputfile,"w")
    figure=[]
    if "-f" in allopt:
        genevalues=bamtocount(bedpos,bamIn,step_f,strands,Bnum)
        SortId=[]
        FpkmValue=[]
        legendname=[]
        groupsize=100//int(groupnum)
        for gz in range(0,int(groupnum)-1):
            lab=str(gz*groupsize)+"-"+str((gz+1)*groupsize)
            legendname.append(lab)
        legendname.append(str((int(groupnum)-1)*groupsize)+"-"+"100")
        fpkm=open(fpkmgroup,"r")
        for line3 in fpkm:
            rowF=line3.rstrip().split()
            SortId.append(str(rowF[0]))
            FpkmValue.append(float(rowF[1]))
        if "-l" in allopt:
            legendname.append("not express")
            hg=[]
            for i in FpkmValue:
                if i > float(minFpkm):
                    hg.append(i)
            Nexpress=len(hg)
            eg=SortId[1:Nexpress]
            lg=SortId[Nexpress:]
            hgsp=splist(eg,int(groupnum))
            N=len(hgsp)
            for l in hgsp:
                allmat=listg(l,genevalues)
                mat=np.array(allmat)
                average=np.nanmean(mat,axis=0)
                print(*average,sep="\t",end="\n",file=out)
                figureline=list(average)
                figure.append(figureline)
            allmat=listg(lg,genevalues)
            mat=np.array(allmat)
            average=np.nanmean(mat,axis=0)
            print(*average,sep="\t",end="\n",file=out)
            figureline=list(average)
            figure.append(figureline)        
        else:   
            hgsp=splist(SortId,int(groupnum))
            for l in hgsp:
                allmat=listg(l,genevalues)
                mat=np.array(allmat)
                average=np.nanmean(mat,axis=0)
                figureline=list(average)
                figure.append(figureline)
                print(*average,sep="\t",end="\n",file=out)
    else:
        legendname=[]
        genevalues=bamtocount(bedpos,bamIn,step_f,strands,Bnum)
        if '-G' in allopt:
            mlist=open(mgenelist,"r")
            print("region","model","value",sep="\t",end="\n",file=out)
            if '--heatmap' in allopt:
                heatout=open(heatfile,"w")
               # print("region","model","value",sep="\t",end="\n",file=heatout)
                signlist={}
                bedlist=[]
                for line in mlist:
                    tmpsign=[]
                    bedlist.append(line.rstrip())
                    allmat=Listgene(line.rstrip(),genevalues)
                    mat=np.array(allmat)
                    number = np.size(mat,1)
                    for j in range(0,number):
			tmpsign.append(list(mat[:,j]))
		    signlist[line.rstrip()]=tmpsign

                    average=np.nanmean(mat,axis=0)
                    directname=line.rstrip().split("/")
                    prefix=directname[-1].replace(".list","")
                    legendname.append(prefix)
                    figureline=list(average)
                    figure.append(figureline)
                  #  for line in allmat:
                     #   colnum=len(line)
                    #    for i in range(0,colnum-1):
                   #         print(line[i],end="\t",file=heatout)
                  #      print(line[colnum-1],end="\n",file=heatout)
                combialist=range(len(bedlist))
                combiafinal=list(combinations(combialist,2))
                for c in combiafinal:
                    n=0
                    for l in range(0,len(signlist[bedlist[0]])):
                        print(bedlist[c[0]]+"vs"+bedlist[c[1]],n*step_f,str(scipy.stats.ranksums(signlist[bedlist[c[0]]][l],signlist[bedlist[c[1]]][l])).split(",")[1].split("=")[1].replace(")",""),sep="\t",end="\n",file=heatout)
		    n=n+1

            N=0
            for lis in figure:
                n = 0
                for l in lis:
                    print(n*step_f,legendname[N],l,sep="\t",end="\n",file=out)
                    n += 1
                N += 1
        else:
            if '-g' in allopt:
                allmat=Listgene(genelist,genevalues)
                directname=genelist.rstrip().split("/")
                prefix=directname[-1].replace(".list","")
                legendname.append(prefix)
            else:
                allmat=list(genevalues.values())
                directname=bamIn.rstrip().split("/")
                prefix=directname[-1].replace(".bam","")
                legendname.append(prefix)                
            mat=np.array(allmat)
            average=np.nanmean(mat,axis=0)
            figureline=list(average)
            figure.append(figureline)
            print(*average,sep="\t",end="\n",file=out)
    if '-F' in allopt:
        for line in figure:
            plotsmooth(line)
        plt.legend(legendname)
        plt.ylabel("Methylation level")
        plt.xticks(labp,labs)
        plt.savefig(figuredir)
else:
    pass

    #print(*allavr,sep="\t",end="\n",file=out)
