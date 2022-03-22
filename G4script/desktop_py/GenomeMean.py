from __future__ import print_function
import sys,pysam
import pybedtools as BedTool
from numpy import array, mean, var

def GWbamTodict(bf):
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

def GWGenomeValue(genomefile,wdsize,readdict):
    VALUE = []
    BedTool.BedTool().window_maker(g=genomefile,w=wdsize).moveto("genome_wd.bed")
    f = open("genome_wd.bed","r")
    genomeList = []
    sumReads = 0
    for line in f:
        col = line.rstrip().split()
        start=int(col[1])
        end=int(col[2])
        length=end-start+1
        if col[0] in readdict:
            for i in range(start+1,end+1):
                if str(i) in readdict[col[0]]:
                    sumReads+=int(readdict[col[0]][str(i)])
                else:
                    pass
        else:
            pass
        genomeList.append(sumReads)
    f.close()
    data = array(genomeList)
    M_value = mean(data)
    V_value = std(data)
    VALUE.extend([M_value,V_value])
    return VALUE

def GWRunCalculate(bf,Gf,ws):
    bamf=pysam.AlignmentFile(bf,"rb")
    region = GWbamTodict(bamf)     
    value = GWGenomeValue(Gf,ws,region)
    return value

