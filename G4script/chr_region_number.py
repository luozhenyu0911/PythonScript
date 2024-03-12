from __future__ import print_function
from __future__ import division
import sys

site = open(sys.argv[1],"r")
bed = open(sys.argv[2],"r")
of = open(sys.argv[3],"w")

allsite = {}
for i in site:
    l = i.rstrip().split("\t")
    mid = int(round((int(l[2])-int(l[1]))/2))+int(l[1])
    if l[0] in allsite:
        pass
    else:
        allsite[l[0]] = []
    allsite[l[0]].append(mid)

#print("type","length",sep="\t",file=of)
for line in bed:
    d = line.rstrip().split("\t")
    number = 0
    for i in allsite[d[0]]:
        if int(d[1]) <= int(i) <= int(d[2]):
            number += 1
        else:
            pass
    print(sys.argv[3].split(".")[0],d[1],number,sep="\t",end="\n",file=of)


site.close()
bed.close()
of.close()
