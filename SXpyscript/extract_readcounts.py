from __future__ import print_function
from __future__ import division
import sys
import getopt
import pysam
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def message():
    print("\n"+"Usage: python reads_density.py --reference-point or --scale-regions [options]"+"\n")
    print("Description:"+"\n"+"\t"+"\t"+"This tool creates scores and profile plot over sets of genomic regions. Typically, these regions are genes, but any other regions defined in BED will work."+"\n")
    print("Options:"+"\n")
    print("\t"+"-h, "+"--help"+"\t"+"show this help message and exit.")
    print("\t"+"-c"+"    <Input bam files>"+"\t"+"Control data files, should be sorted bam format.")
    print("\t"+"-C"+"    <list file>"+"\t"+"multiple control data files list in a text file, one bam file in a line, should be sorted bam format.")
    print("\t"+"-b"+"    <bamfiles>"+"\t"+"data files, should be sorted bam format.")
    print("\t"+"-B"+"    <list file>"+"\t"+"multiple input data files list in a text file, one bam file in a line, should be sorted bam format.")
    print("\t"+"-r"+"    <bedfiles>"+"\t"+"data files, should be original bed files with 6 columns.")
    print("\t"+"-g"+"    <gene list files>"+"\t"+"data file must end with .list, should be a gene ID each row.")
    print("\t"+"-G"+"    <list file>"+"\t"+"multiple input data files list in a text file, each data file must end with .list, should be a gene ID each row.")
    print("\t"+"-m"+"    <int>"+"\t"+"Number of blocks genebody will be divided.")
    print("\t"+"-u"+"    <upstream distance>"+"\t"+"Distance upstream of the reference-point selected.")
    print("\t"+"-d"+"    <downstream distance>"+"\t"+"Distance downstream of the reference-point selected.")
    print("\t"+"-w"+"    <window size>"+"\t"+"Length, in bases, binSize for averaging the score over the regions length.")
    print("\t"+"-n"+"    <int num>"+"\t"+"number of groups divided according to the expression levels of all expressed genes")
    print("\t"+"-o"+"    <outputprefix or ouputfilename>"+"\t"+"File name to save the average data file. If mutiple files were submitted, the output file will use this string as prefix.")
    print("\t"+"-F"+"    <Figure name>"+"\t"+"Name of the plot name. (plot file format is PDF)"+"\n")
    print("\t"+"--grouprev"+"\t"+"If mutiple bamfiles and gene lists were subject, use this option to plot different lists of different genes for each bamfiles"+"\n")
    print("An example usage is:"+"\n")
    print("\t"+"python reads_density.py -b **.bam -r **.bed -g **.list -m TSS -u 1000 -d 1000 -w 50 -o **.txt -F **.pdf"+"\n")
    print("\t"+"python reads_density.py -b **.bam -r **.bed -g **.list -u 1000 -d 1000 -w 50 -N 100 -l 0.5 -f **.fpkm -n 5 -o **.txt -F **.pdf"+"\n")
    return

opts,args = getopt.getopt(sys.argv[1:],'-h-b:-r:-G:-g:-u:-d:-o:-w:-B:-f:-n:-l:-F:-c:-C:-m:',['help','grouprev','heatmap='])
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
    if opt_name in ('-B'):
        bamIns = opt_value
    if opt_name in ('-n'):
        groupnum = opt_value
    if opt_name in ('-c'):
        inIn = opt_value
    if opt_name in ('-C'):
        inIns = opt_value
    if opt_name in ('-m'):
        Bnum = int(opt_value)

allopt=[]
for option in opts:
    allopt.append(option[0])

def bamInfo(Input):
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

def Readscount(bedlist,readdict,rn):
    allgene={}
    for row2 in bedlist:
        sumReads=0
        if row2[0] in readdict:
            for k in range(int(row2[1]),int(row2[2]) + 1):
                if str(k) in readdict[row2[0]]:
                    sumReads+=int(readdict[row2[0]][str(k)])
                else:
                    pass
                norReads=(sumReads/(rn*(int(row2[2])-int(row2[1]))))*1000000000
                values=norReads
        else:
            values=0
        values=float(values)
        allgene[str(row2[3])]=values
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

def bamtocount(bed,Input):
    bamfile=pysam.AlignmentFile(Input,"rb")
    region=bamInfo(bamfile)
    N=bamfile.count()
    geneall=Readscount(bed,region,N)
    return geneall
    
bedfile=open(bedIn,"r")
bedpos=[]
for genepos in bedfile:
    row=genepos.rstrip().split()
    [chrom,start,end,geneID]=row
    bedpos.append([chrom,start,end,geneID])
 
if "-b" in allopt:
    out=open(outputfile,"w")
    if "-c" in allopt:
        inputvalues=bamtocount(bedpos,inIn)
        controlmat=list(inputvalues.values())
        genevalues=bamtocount(bedpos,bamIn)
        allmat=list(genevalues.values())
        tmat=np.array(allmat)
        cmat=np.array(controlmat)
        mat=tmat-cmat
        mat[mat<0]=0
        mat=np.array(list(mat))
        keys=genevalues.keys()
        for i,j,k,l in zip(keys,list(cmat),list(tmat),list(mat)):
            print(i,j,k,l,sep="\t",end="\n",file=out)
        out.close()
else:
    pass
