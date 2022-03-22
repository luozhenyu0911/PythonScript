from __future__ import print_function 
from __future__ import division
import sys
import pysam
bf=pysam.AlignmentFile(sys.argv[1],"rb")
bed=open(sys.argv[2],"r")

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
	
sumReads=0
count = (bf.count())/10000000
#for l in bedlist:
 #   bedf = l.rstrip()
  #  legend = bedf[:-4]
   # bed = open(bedf,"r")
for line2 in bed:
    col=line2.rstrip().split()
    start=int(col[1])
    end=int(col[2])
    length=end-start+1
    if col[0] in region:
        for i in range(start+1,end+1):
            if str(i) in region[col[0]]:
                sumReads+=int(region[col[0]][str(i)])
            else:
                pass
    bizhi = float(sumReads/(length*count))
    print(sys.argv[2],bizhi,end="\n",sep="\t") 
    sumReads=0
bed.close()
#bedlist.close()
