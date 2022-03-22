from __future__ import print_function
from __future__ import division
from pybedtools import BedTool
import pybedtools as pybt
import sys
import pysam
import numpy as np
import subprocess
import scipy
import gzip
import scipy.stats



def bamInfo(Input):
    region={}
    for read in Input:
        col=read.rstrip().split()
        if int(col[7]) >= 5:
            if str(col[0]) in region:
                region[str(col[0])][str(col[1])]=(int(col[6]),int(col[7]))
            else:
                region[str(col[0])]={}
                region[str(col[0])][str(col[1])]=(int(col[6]),int(col[7]))
    return region


def Readscount(bedlist,readdict):
    metC=0
    sumC=0
    values=[]
    for row in bedlist:
        row2=row.rstrip().split("\t")
        if row2[0] in readdict:
            for k in range(int(row2[1]),int(row2[2])):
                if str(k) in readdict[row2[0]]:
                    sumC+=readdict[row2[0]][str(k)][1]
                    metC+=readdict[row2[0]][str(k)][0]
                else:
                    pass
            if sumC == 0:
                methyratio=0
                values.append(methyratio)
            else:
                methyratio=metC/sumC
                values.append(methyratio)
        else:
            pass
        metC=0
        sumC=0
    return values

bamfile=gzip.open(sys.argv[1],"rb")
bed1=open(sys.argv[2],"r")

region=bamInfo(bamfile)
obvalue=Readscount(bed1,region)

tmpplist=[]
exvalue = []
name = sys.argv[1].split(".")[0]
InBed=BedTool(sys.argv[2])
for i in range(0,10):
    InBed.shuffle(g=sys.argv[3],noOverlapping=True).moveto(name+"tmp.bed")
    randomb=open(name+"tmp.bed","r")
    extmpv=Readscount(randomb,region)
    extm=np.array(extmpv)
    tmpmean=np.mean(extm)
    exvalue.append(tmpmean)
    tmpplist.append(extmpv)
    randomb.close()
    subprocess.call(["rm",name+"tmp.bed"])

ob=np.array(obvalue)
ex=np.array(exvalue)
mean1=np.mean(ob)
mean2=np.mean(ex)

print(mean1/mean2)
aa=np.array(tmpplist)
average1=list(np.mean(aa,axis=0))
print(scipy.stats.mannwhitneyu(obvalue,average1))
print(scipy.stats.ranksums(obvalue,average1))
bamfile.close()
bed1.close()
