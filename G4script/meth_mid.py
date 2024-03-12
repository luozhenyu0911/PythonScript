from __future__ import print_function
from __future__ import division
from pyfasta import Fasta
from pybedtools import BedTool
from itertools import combinations
import pybedtools as pybt
import re
import sys
import subprocess
import numpy as np
import scipy
import scipy.stats
#import matplotlib.pyplot as plt

meth = open(sys.argv[1],"r")
bedlist=sys.argv[2].split(",")
wdsize = int(sys.argv[3])
up = int(sys.argv[4])
down=int(sys.argv[5])
of1=open(sys.argv[6],"w")
of2=open(sys.argv[7],"w")

methl = {}
for line in meth:
    l = line.rstrip().split("\t")
    if int(l[7]) > 5:
        if str(l[0]) in methl:
            methl[str(l[0])][str(l[1])]=(int(l[6]),int(l[7]))
        else:
            methl[str(l[0])]={}
            methl[str(l[0])][str(l[1])]=(int(l[6]),int(l[7]))
    else:
        pass

figurelist=[]
signlist={}
for bedd in bedlist:
    alllist=[]
    tmpsign=[]
    bed = open(bedd,"r")
    for Line in bed:
        l = Line.rstrip().split("\t")
        chorm = str(l[0])
        start = int(l[1])
        end = int(l[2])
        value=[]
        for i in range(start-up,start+down,wdsize):
            sumC=0
            metC=0
            for j in range(i,(i+wdsize)):
                if str(j) in methl[l[0]]:
                    sumC+=methl[l[0]][str(j)][1]
                    metC+=methl[l[0]][str(j)][0]
            if sumC == 0:
                value.append(np.nan)
            else:
                methyratio=metC/sumC
                value.append(methyratio)
        alllist.append(value)
    dd = np.array(alllist)
    number = np.size(dd,1)
    for j in range(0,number):
        tmpsign.append(list(dd[:,j]))
    signlist[bedd]=tmpsign
    gg = np.nanmean(dd,axis=0)
    figurelist.append(gg)
    bed.close()
#print(len(signlist[bedlist[0]][0]),len(signlist[bedlist[1]][0]),len(signlist[bedlist[2]][0]))
print("region","model","value",sep="\t",end="\n",file=of1)
N = 0
for lis in figurelist:
    n = 0
    for l in lis:
        print(n*wdsize,bedlist[N],l,sep="\t",end="\n",file=of1)
        n += 1
    N += 1

print("region","model","value",sep="\t",end="\n",file=of2)
combialist=range(len(bedlist))
combiafinal=list(combinations(combialist,2))
for c in combiafinal:
    n=0
    for l in range(0,len(signlist[bedlist[0]])):
        print(bedlist[c[0]]+"vs"+bedlist[c[1]],n*wdsize,str(scipy.stats.ranksums(signlist[bedlist[c[0]]][l],signlist[bedlist[c[1]]][l])).split(",")[1].split("=")[1].replace(")",""),sep="\t",end="\n",file=of2)
        n=n+1


of1.close()
of2.close()
meth.close()
