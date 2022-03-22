from __future__ import print_function 
from __future__ import division
from pybedtools import BedTool
import sys
import subprocess
import numpy

#"---observed value---"
fa=open(sys.argv[1].split(".")[0]+".fa","w")
subprocess.call(["bedtools","getfasta","-fi","/gss1/home/zqq20190923/index/Nitab-v4.5_genome_Chr_Edwards2017.fasta","-bed",sys.argv[1]],stdout=fa)
subprocess.call(["sh",sys.argv[2]+".sh",str(sys.argv[1].split(".")[0]+".fa"),str(sys.argv[1].split(".")[0]+sys.argv[2]+".bed")])
tmplist=[]
tmpf=open(sys.argv[1].split(".")[0]+sys.argv[2]+".bed","r")
for line in tmpf:
    l = line.rstrip().split()
    tmplist.append(l[0])
ob = len(tmplist)
fa.close()
tmpf.close()
#subprocess.call(["rm",sys.argv[1].split(".")[0]+sys.argv[2]+"tmp"])

#"---expected value---"
excepted = []
#name = sys.argv[1].split(".")[0]
InBed=BedTool(sys.argv[1])
for i in range(0,3):
    InBed.shuffle(g=sys.argv[3],noOverlapping=True).moveto(sys.argv[1].split(".")[0]+sys.argv[2]+"tmp.bed")
    tmpfa=open(sys.argv[1].split(".")[0]+sys.argv[2]+"tmp.fa","w")
    subprocess.call(["bedtools","getfasta","-fi","/gss1/home/zqq20190923/index/Nitab-v4.5_genome_Chr_Edwards2017.fasta","-bed",str(sys.argv[1].split(".")[0]+sys.argv[2]+"tmp.bed")],stdout=tmpfa)
    subprocess.call(["sh",str(sys.argv[2]+".sh"),str(sys.argv[1].split(".")[0]+sys.argv[2]+"tmp.fa"),str(sys.argv[1].split(".")[0]+sys.argv[2]+"tmp")])
    tmpplist=[]
    tmpff=open(sys.argv[1].split(".")[0]+sys.argv[2]+"tmp","r")
    for line in tmpff:
        l = line.rstrip().split()
        tmpplist.append(l[0])
    exv = len(tmpplist)
    excepted.append(exv)
    tmpfa.close()
    tmpff.close()
    subprocess.call(["rm",str(sys.argv[1].split(".")[0]+sys.argv[2]+"tmp.fa")])
    subprocess.call(["rm",str(sys.argv[1].split(".")[0]+sys.argv[2]+"tmp.bed")])
    subprocess.call(["rm",str(sys.argv[1].split(".")[0]+sys.argv[2]+"tmp")])
e=numpy.array(excepted)
exmean=e.sum()/len(excepted)
print(ob,exmean,ob/exmean,sep="\t",end="\n")    
    
    
    
    
    
    
    
