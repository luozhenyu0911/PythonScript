from __future__ import print_function
import sys
import re

a=open(sys.argv[1],"r")
out=open(sys.argv[2],"w")

for line in a:
    if sys.argv[3] in line:
        row=line.rstrip().split("\t")
        row2=re.split(r"[:-]",row[2])
        motifS=int(row2[3])+int(row[3])-1
        motifE=int(row2[3])+int(row[4])-1
        print(row2[2],motifS,motifE,row2[0],row[9],row[5],sep="\t",end="\n",file=out)
