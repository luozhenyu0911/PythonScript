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
    print("\t"+"-f"+"    <fpkm file>"+"\t"+"Text file contains two columns, GeneIDs and expression values.")
    print("\t"+"-l"+"    <float num>"+"\t"+"Threshold for not express group, genes with a FPKM lower than this value were considered as not expressed ")
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
        inIn = opt_value
    if opt_name in ('-C'):
        inIns = opt_value
    if opt_name in ('-m'):
        Bnum = int(opt_value)
    if opt_name in ('--heatmap'):
        heatfile = opt_value
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

def bamInfo(bf):
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

def Readscount(bedlist,readdict,rn,windows,BlockNum):
    sumReads=0
    allgene={}
    revwins=0-windows
    for row2 in bedlist:
        values=[]
        if str(row2[6]) is "+":
            if row2[0] in readdict:
                for k in range(int(row2[1]),int(row2[2]),windows):
                    for j in range(k,(k+windows)):
                        if str(j) in readdict[row2[0]]:
                            sumReads+=int(readdict[row2[0]][str(j)])
                        else:
                            pass
                    norReads=sumReads/(rn*windows)*10000000
                    values.append(norReads)
                    sumReads=0
            else:
                for j in range(int(row2[1]),int(row2[2]),windows):
                    values.append("0")
            if row2[0] in readdict:
                reside=BlockNum-int(((int(row2[3])-int(row2[2]))%BlockNum))
                blockSize=int(((int(row2[3])-int(row2[2]))+reside)/BlockNum)
                for x in range(int(row2[2]),int(row2[3])+reside,blockSize):
                    for y in range(x,(x+blockSize)):
                        if str(y) in readdict[row2[0]]:
                            sumReads+=int(readdict[row2[0]][str(y)])
                        else:
                            pass
                    values.append((int(sumReads)*10000000)/(rn*blockSize))
                    sumReads=0
            else:
                for x in range(int(row2[2]),int(row2[3])+reside,blockSize):
                    values.append("0")
            if row2[0] in readdict:
                for k in range(int(row2[3]),int(row2[4]),windows):
                    for j in range(k,(k+windows)):
                        if str(j) in readdict[row2[0]]:
                            sumReads+=int(readdict[row2[0]][str(j)])
                        else:
                            pass
                    norReads=sumReads/(rn*windows)*10000000
                    values.append(norReads)
                    sumReads=0
            else:
                for j in range(int(row2[3]),int(row2[4]),windows):
                    values.append("0")
        else:
            if row2[0] in readdict:
                for k in range(int(row2[4]),int(row2[3]),revwins):
                    for j in range(k,(k+revwins),-1):
                        if str(j) in readdict[row2[0]]:
                            sumReads+=int(readdict[row2[0]][str(j)])
                        else:
                            pass
                    norReads=sumReads/(rn*windows)*10000000
                    values.append(norReads)
                    sumReads=0
            else:
                for j in range(int(row2[4]),int(row2[3]),revwins):
                    values.append("0")
            if row2[0] in readdict:
                reside=BlockNum-int(((int(row2[3])-int(row2[2]))%BlockNum))
                blockSize=int(((int(row2[3])-int(row2[2]))+reside)/BlockNum)
                reblockSize=0-blockSize
                for x in range(int(row2[3]),int(row2[2])-reside,reblockSize):
                    for y in range(x,(x-blockSize),-1):
                        if str(y) in readdict[row2[0]]:
                            sumReads+=int(readdict[row2[0]][str(y)])
                        else:
                            pass
                    values.append((int(sumReads)*10000000)/(rn*blockSize))
                    sumReads=0
            else:
                for x in range(int(row2[3]),int(row2[2])-reside,reblockSize):
                    values.append("0")
            if row2[0] in readdict:
                for k in range(int(row2[2]),int(row2[1]),revwins):
                    for j in range(k,(k+revwins),-1):
                        if str(j) in readdict[row2[0]]:                                                  sumReads+=int(readdict[row2[0]][str(j)])
                        else:
                            pass
                    norReads=sumReads/(rn*windows)*10000000
                    values.append(norReads)
                    sumReads=0
            else:
                for j in range(int(row2[2]),int(row2[1]),revwins):
                    values.append("0")
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


def bamtocount(bed,Input,windowsize,BN):
    bamfile=pysam.AlignmentFile(Input,"rb")
    region=bamInfo(bamfile)
    N=bamfile.count()
    geneall=Readscount(bed,region,N,windowsize,BN)
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
m1lab="TSS"
m2lab="TTS"
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
    if "-c" in allopt:
        inputvalues=bamtocount(bedpos,inIn,step_f,Bnum)
    if "-f" in allopt:
        genevalues=bamtocount(bedpos,bamIn,step_f,Bnum)
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
                if "-c" in allopt:
                    allmat=listg(l,genevalues)
                    controlmat=listg(l,inputvalues)
                    mat=np.array(allmat)
                    cmat=np.array(controlmat)
                    mat=mat-cmat
                    mat[mat<0]=0
                    allmat=list(mat)
                else:
                    allmat=listg(l,genevalues)
                    mat=np.array(allmat)
                average=mat.mean(axis=0)
                print(*average,sep="\t",end="\n",file=out)
                figureline=list(average)
                figure.append(figureline)
            allmat=listg(lg,genevalues)
            mat=np.array(allmat)
            average=mat.mean(axis=0)
            print(*average,sep="\t",end="\n",file=out)
            figureline=list(average)
            figure.append(figureline)        
        else:   
            hgsp=splist(SortId,int(groupnum))
            for l in hgsp:
                if "-c" in allopt:
                    allmat=listg(l,genevalues)
                    controlmat=listg(l,inputvalues)
                    mat=np.array(allmat)
                    cmat=np.array(controlmat)
                    mat=mat-cmat
                    mat[mat<0]=0
                    allmat=list(mat)
                else:
                    allmat=listg(l,genevalues)
                    mat=np.array(allmat)
                average=mat.mean(axis=0)
                figureline=list(average)
                figure.append(figureline)
                print(*average,sep="\t",end="\n",file=out)
    else:
        legendname=[]
        genevalues=bamtocount(bedpos,bamIn,step_f,Bnum)
        if '-G' in allopt:
            mlist=open(mgenelist,"r")
            for line in mlist:
                if "-c" in allopt:
                    allmat=Listgene(line.rstrip(),genevalues)
                    controlmat=Listgene(line.rstrip(),inputvalues)
                    mat=np.array(allmat)
                    cmat=np.array(controlmat)
                    mat=mat-cmat
                    mat[mat<0]=0
                    allmat=list(mat)
                else:
                    allmat=Listgene(line.rstrip(),genevalues)
                    mat=np.array(allmat)
                average=mat.mean(axis=0)
                print(*average,sep="\t",end="\n",file=out)
                directname=line.rstrip().split("/")
                prefix=directname[-1].replace(".list","")
                legendname.append(prefix)
                figureline=list(average)
                figure.append(figureline)
        else:
            if '-g' in allopt:
                if "-c" in allopt:
                    allmat=Listgene(genelist,genevalues)
                    controlmat=Listgene(genelist,inputvalues)
                    mat=np.array(allmat)
                    cmat=np.array(controlmat)
                    mat=mat-cmat
                    mat[mat<0]=0
                    allmat=list(mat)
                else:
                    allmat=Listgene(genelist,genevalues)
                directname=genelist.rstrip().split("/")
                prefix=directname[-1].replace(".list","")
                legendname.append(prefix)
            else:
                if "-c" in allopt:
                    allmat=list(genevalues.values())
                    controlmat=list(inputvalues.values())
                    mat=np.array(allmat)
                    cmat=np.array(controlmat)
                    mat=mat-cmat
                    mat[mat<0]=0
                    allmat=list(mat)
                else:
                    allmat=list(genevalues.values())
                directname=bamIn.rstrip().split("/")
                prefix=directname[-1].replace(".bam","")
                legendname.append(prefix)                
            mat=np.array(allmat)
            if '--heatmap' in allopt:
                heatout=open(heatfile,"w")
                for line in allmat:
                    colnum=len(line)
                    for i in range(0,colnum-1):
                        print(line[i],end="\t",file=heatout)
                    print(line[colnum-1],end="\n",file=heatout)
            average=mat.mean(axis=0)
            figureline=list(average)
            figure.append(figureline)
            print(*average,sep="\t",end="\n",file=out)
    if '-F' in allopt:
        for line in figure:
            plotsmooth(line)
        plt.legend(legendname)
        plt.ylabel("Normalized reads count")
        plt.xticks(labp,labs)
        plt.savefig(figuredir)
else:
    pass

if "--grouprev" in allopt:
    if '-B' in allopt:
        bamdict=[]
        bamlist=open(bamIns,"r")
        legendname=[]
        outlist=[]
        ptitle=[]
        for line2 in bamlist:
            bamdir=line2.rstrip()
            genevalues=bamtocount(bedpos,bamdir,step_f,Bnum)
            directname=bamdir.split("/")
            prefix=directname[-1].replace(".bam","")
            outputname=outputfile+"_"+prefix+".txt"
            outlist.append(outputname)
            bamdict.append(genevalues)
            ptitle.append(prefix)
        if "-C" in allopt:
            inputdict=[]
            inputlist=open(inIns,"r")
            for line2 in inputlist:
                inputdir=line2.rstrip()
                inputvalues=bamtocount(bedpos,inputdir,step_f,Bnum)
                inputdict.append(inputvalues)
        bamnum=len(bamdict)
        mutifigure=[]
        for i in range(0,bamnum):
            figure=[]
            output=open(outlist[i],"w")
            if '-G' in allopt:
                mlist=open(mgenelist,"r")
                for line in mlist:
                    indir=line.rstrip()
                    directory=line.rstrip().split("/")
                    prefix=directory[-1].replace(".list","")
                    legendname.append(prefix)
                    if "-C" in allopt:
                        allmat=Listgene(indir,bamdict[i])
                        controlmat=Listgene(indir,inputdict[i])
                        mat=np.array(allmat)
                        cmat=np.array(controlmat)
                        mat=mat-cmat
                        mat[mat<0]=0
                        allmat=list(mat)
                    else:
                        allmat=Listgene(indir,bamdict[i])
                        mat=np.array(allmat)
                    average=mat.mean(axis=0)
                    figureline=list(average)
                    figure.append(figureline)
                    print(*average,end="\n",sep="\t",file=output)
                output.close()
                mutifigure.append(figure)
            else:
                pass
        if '-F' in allopt:
            plotnum=len(mutifigure)
            plt.figure(figsize=(15,(plotnum//3+1)*5))
            for i in range(0,plotnum):
                plt.subplot((plotnum//3+1),3,(i+1))
                for line in mutifigure[i]:
                    plotsmooth(line)
                plt.legend(legendname)
                plt.title(ptitle[i])
                plt.ylabel("Normalized reads count")
                plt.xticks(labp,labs)
            plt.savefig(figuredir)
    else:
        pass

else:
    if '-B' in allopt:
        bamdict=[]
        bamlist=open(bamIns,"r")
        legendname=[]
        for line2 in bamlist:
            bamdir=line2.rstrip()
            genevalues=bamtocount(bedpos,bamdir,step_f,Bnum)
            directname=bamdir.split("/")
            prefix=directname[-1].replace(".bam","")
            legendname.append(prefix)
            bamdict.append(genevalues)
        if "-C" in allopt:
            inputdict=[]
            inputlist=open(inIns,"r")
            for line2 in inputlist:
                inputdir=line2.rstrip()
                inputvalues=bamtocount(bedpos,inputdir,step_f,Bnum)
                inputdict.append(inputvalues)
        if '-G' in allopt:
            mutifigure=[]
            ptitle=[]
            mlist=open(mgenelist,"r")
            for line in mlist:
                figure=[]
                indir=line.rstrip()
                directory=line.rstrip().split("/") 
                prefix=directory[-1].replace(".list","")
                outputname=outputfile+"_"+prefix+".txt"
                output=open(outputname,"w")
                bamnum=len(bamdict)
                ptitle.append(prefix)
                for i in range(0,bamnum):
                    if "-C" in allopt:
                        allmat=Listgene(indir,bamdict[i])
                        controlmat=Listgene(indir,inputdict[i])
                        mat=np.array(allmat)
                        cmat=np.array(controlmat)
                        mat=mat-cmat
                        mat[mat<0]=0   
                        allmat=list(mat)
                    else:
                        allmat=Listgene(indir,bamdict[i])
                        mat=np.array(allmat)
                    average=mat.mean(axis=0)
                    figureline=list(average)
                    figure.append(figureline)
                    print(*average,end="\n",sep="\t",file=output)
                output.close()
                mutifigure.append(figure)
            if '-F' in allopt:
                plotnum=len(mutifigure)
                plt.figure(figsize=(15,((plotnum+1)//3+1)*5))
                for i in range(0,plotnum):
                    plt.subplot(((plotnum+1)//3+1),3,(i+1))
                    for line in mutifigure[i]:
                        plotsmooth(line)
                    plt.legend(legendname)
                    plt.title(ptitle[i])
                    plt.ylabel("Normalized reads count")
                    plt.xticks(labp,labs)
                plt.savefig(figuredir)
        else:
            figure=[]
            out=open(outputfile,"w")
            bamnum=len(bamdict)
            for i in range(0,bamnum):
                if '-g' in allopt:
                    if "-C" in allopt:
                        allmat=Listgene(genelist,bamdict[i])
                        controlmat=Listgene(genelist,inputdict[i])
                        mat=np.array(allmat)
                        cmat=np.array(controlmat)
                        mat=mat-cmat
                        mat[mat<0]=0
                        allmat=list(mat)
                    else:
                        allmat=Listgene(genelist,bamdict[i])
                        mat=np.array(allmat)
                else:
                    if "-C" in allopt:
                        allmat=list(bamdict[i].values())
                        controlmat=list(inputdict[i].values())
                        mat=np.array(allmat)
                        cmat=np.array(controlmat)
                        mat=mat-cmat
                        mat[mat<0]=0
                        allmat=list(mat)
                    else:
                        allmat=list(bamdict[i].values())
                        mat=np.array(allmat)
                mat=np.array(allmat)
                average=mat.mean(axis=0)
                figureline=list(average)
                figure.append(figureline)
                print(*average,end="\n",sep="\t",file=out)
            if '-F' in allopt:
                for line in figure:
                    plotsmooth(line)
                plt.legend(legendname)
                plt.ylabel("Normalized reads count")
                plt.xticks(labp,labs)
                plt.savefig(figuredir)



    #print(*allavr,sep="\t",end="\n",file=out)
