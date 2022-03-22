from __future__ import print_function
from __future__ import division
import sys
import pysam 

bf=pysam.AlignmentFile(sys.argv[1],"rb")
bed = open(sys.argv[2],"r")
wdsize = int(sys.argv[3])

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

bf.close()
#allread = (int(bf.count()))/1000000
print("group","type","length",sep="\t")

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
        print("G4",i,bizhi,sep="\t")
        num=0
    num=0
    for j in range(end,length):
        if str(j) in region[l[0]]:
            num += int(region[l[0]][str(j)])
    if num > 0:
        bizhi = num/wdsize2
    else:
        bizhi=0
    print("G4",length,bizhi,sep="\t")
         
bed.close()
