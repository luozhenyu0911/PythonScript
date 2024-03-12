from __future__ import print_function
from __future__ import division
import sys
import pysam 

bf=open(sys.argv[1],"r")
bed = open(sys.argv[2],"r")
wdsize = int(sys.argv[3])

region={}

for line in bf:
    l = line.rstrip().split()
    if l[0] in region:
    else:
        region[l[0]] = {}




bf.close()
#allread = (int(bf.count()))/1000000
print("type","length",sep="\t")

for Line in bed:
    l = Line.rstrip().split("\t")
    chorm = str(l[0])
    length = int(l[1])
    dev = length//wdsize
    start=0
    end=length-(length%wdsize)
    wdsize2 = length%wdsize
    num=0
    for i in range(start,end,wdsize):
        for j in range(i,(i+wdsize)):
            if str(j) in region[l[0]]:
                num += int(region[l[0]][str(j)])
        if num > 0:
            bizhi = num/wdsize
        else:
            bizhi = 0
        print(i,bizhi,sep="\t")
        num=0
    num=0
    for j in range(end,length):
        if str(j) in region[l[0]]:
            num += int(region[l[0]][str(j)])
    if num > 0:
        bizhi = num/wdsize2
    else:
        bizhi=0
    print(dev,bizhi,sep="\t")
         
bed.close()
