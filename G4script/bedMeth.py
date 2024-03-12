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
    CG={}
    CHG={}
    CHH={}
    C={}
    for line in Input:
        col=line.rstrip().split()
        chr=str(col[0])
        allcontext=str(col[1])
        site=str(col[2])
        context=str(col[3])
        allC=int(col[7])
        methyC=int(col[6])
        if allC >= 2:
            if chr in C:
                C[chr][site]=(methyC,allC)
            else:
                C[chr]={}
                C[chr][site]=(methyC,allC)
            if context == 'CG':
                if chr in CG:
                    CG[chr][site]=(methyC,allC)
                else:
                    CG[chr]={}
                    CG[chr][site]=(methyC,allC)
            elif context == 'CHG':
                if chr in CHG:
                    CHG[chr][site]=(methyC,allC)
                else:
                    CHG[chr]={}
                    CHG[chr][site]=(methyC,allC)
            elif context == 'CHH':
                if chr in CHH:
                    CHH[chr][site]=(methyC,allC)
                else:
                    CHH[chr]={}
                    CHH[chr][site]=(methyC,allC)
            else:
                pass
    region={'C':C,'CG':CG,'CHG':CHG,'CHH':CHH}
    return region
def Readscount(bedlist,readdict):
    metC=0
    sumC=0
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
    else:
        methyratio=metC/sumC
    return methyratio 
context_l=['C','CG','CHG','CHH']
for Bam in range(len(bamlist)):
    bamf=gzip.open(bamlist[Bam],'rb')
    MethyDict = CGmapTodict(bamf)
    for Bed in range(len(bedlist)):
        for context in context_l:
            with open(bedlist[Bed],'r') as bedfile:
                values=Readscount(bedfile,MethyDict[context])
                print(bamlist[Bam].strip().split('.')[0],bedlist[Bed].strip().split('.')[0],context,values,sep="\t",file=outf)

methf.close()
bedf.close()
outf.close()
