'''
python this.py methgzlist bedlist outfile
'''
from __future__ import print_function
from __future__ import division
import sys
import gzip
# import subprocess
# import scipy
# import scipy.stats

methf=open(sys.argv[1],'r')
bedf=open(sys.argv[2],'r')
outf=open(sys.argv[3],'w')
bamlist=[]
for i in methf:
    bamlist.append(i.strip())
bedlist=[]
for i in bedf:
    bedlist.append(i.strip())

def CGmapTodict(Input):
    region={}
    for line in Input:
        col=line.rstrip().split()
        chr=str(col[0])
        # allcontext=str(col[1])
        site=str(col[2])
        # context=str(col[3])
        allC=int(col[7])
        methyC=int(col[6])
        if allC >= 5:
            if chr in region:
                region[chr][site]=(methyC,allC)
            else:
                region[chr]={}
                region[chr][site]=(methyC,allC)
    return region

def Readscount(bed,readdict):
    # values=[]
    allnum=0
    validnum=0
    for row in bed:
        metC=0
        sumC=0
        allnum+=1
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
                # values.append(methyratio)
            else:
                methyratio=metC/sumC
                # values.append(methyratio)
            if methyratio > 0.1:
                validnum+=1
            else:
                pass
        else:
            pass
    return allnum,validnum

for Bam in range(len(bamlist)):
    bamf=gzip.open(bamlist[Bam],'rb')
    MethyDict = CGmapTodict(bamf)
    for Bed in range(len(bedlist)):
        with open(bedlist[Bed],'r') as bedfile:
            allnum,validnum=Readscount(bedfile,MethyDict)
            print(bedlist[Bed].strip().split('.')[0],allnum,validnum,sep="\t",file=outf)
                

methf.close()
bedf.close()
outf.close()
