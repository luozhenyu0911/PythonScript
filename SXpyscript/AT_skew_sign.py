# -*- coding: UTF-8 -*-
from __future__ import print_function
from __future__ import division
from pyfasta import Fasta
import re
import subprocess
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
plt.switch_backend('PDF')

parser = argparse.ArgumentParser(description="This tool used to caculate AT skew for genome",formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument('--bedfiles','-b',type = str,help="-b <str>    Input bed files,(Support for multiple files, Separated by commas)",required = True)
parser.add_argument('--fastafile','-F',type = str,help="-F <str>    Input genome fasta file",required = True)
parser.add_argument('--wdsize','-w',type = int,help="-w <int>    Window size(can be 10bp)",required = True)
parser.add_argument('--blocknum','-n',type = int,help="-n <str>    how many block for region",required = True)
parser.add_argument('--upstream','-u',type = int,help="-u <int>    region upstream",required = True)
parser.add_argument('--downstream','-d',type = int,help="-d <int>    region downstream",required = True)
parser.add_argument('--strand','-s',type = str,help="-s <str>    strand should be (Sign),this is only complement sequence in another strand(not like TSS direction)\n"
                    "\tNote: If <Sign> is specified, which strand be used depends on the + or - in 6th of bed file, + means forward strand, - means reverse strand.",required = True)
parser.add_argument('--legend','-l',type = str,help="-l <str>    plot legends(Support for multiple files, Separated by commas)",required = True)
parser.add_argument('--outfile','-o',type = str,help="-o <str>    Output AT skew file",required = True)
parser.add_argument('--figurefile','-f',type = str,help="-f <str>    Output figure file",required = True)
args = parser.parse_args()

bedlist=(args.bedfiles).split(",")
legend=(args.legend).split(",")

def plotsmooth(Values):
    n=len(Values)
    num=range(0,n)
    num=np.array(num)
    xnew=np.linspace(num.min(),num.max(),5*n)
    func=interp1d(num,Values,kind="cubic")
    ynew=func(xnew)
    plt.plot(xnew,ynew)

def GCValue(bedfile,wdsize,fastafile,BlockNum,up,down,strand):
    fa = Fasta(fastafile)
    f = open(bedfile,"r")
    left=int(0)
    mediumA=int(up/wdsize)
    mediumB=mediumA+BlockNum
    right=int(down/wdsize)+mediumA+BlockNum
    alllist = []
    for Line in f:
        L = []
        M = []
        R = []
        line=Line.rstrip().split()
        reside=BlockNum-int(((int(line[2])-int(line[1]))%BlockNum))
        blockSize=int(((int(line[2])-int(line[1]))+reside)/BlockNum)
        sublist=[]
        if str(strand) == "Sign":
            sign=line[5]
            if sign == "+":
                for l in range(int(int(line[1])-up),int(line[1]),wdsize):
                    left = str(l)
                    right = str(l+wdsize)
                    a = str("{'chr':"+"'"+line[0]+"'"+","+"'start':"+left+","+"'stop':"+right+"}")
                    L.append(a)
                for l in range(int(line[1]),int(line[2])+reside,blockSize):
                    left = str(l)
                    right = str(l+blockSize)
                    a = str("{'chr':"+"'"+line[0]+"'"+","+"'start':"+left+","+"'stop':"+right+"}")
                    M.append(a)
                for l in range(int(line[2]),int(int(line[2])+down),wdsize):     
                    left = str(l)
                    right = str(l+wdsize)
                    a = str("{'chr':"+"'"+line[0]+"'"+","+"'start':"+left+","+"'stop':"+right+"}")
                    R.append(a)
                for i in L:
                    convertI = eval(i)
                    seq = fa.sequence(eval(i), one_based=False)
                    Anum = len(re.findall("A",seq,re.IGNORECASE))
                    Tnum = len(re.findall("T",seq,re.IGNORECASE))
                    if Anum == 0 and Tnum == 0:
                        value = 0
                    else:
                        value = float((Anum-Tnum)/(Anum+Tnum))
                    sublist.append(value)
                for i in M:
                    convertI = eval(i)
                    seq = fa.sequence(eval(i), one_based=False)
                    Anum = len(re.findall("A",seq,re.IGNORECASE))
                    Tnum = len(re.findall("T",seq,re.IGNORECASE))
                    if Anum == 0 and Tnum == 0:
                        value = 0
                    else:
                        value = float((Anum-Tnum)/(Anum+Tnum))
                    sublist.append(value)
                for i in R:
                    convertI = eval(i)
                    seq = fa.sequence(eval(i), one_based=False)
                    Anum = len(re.findall("A",seq,re.IGNORECASE))
                    Tnum = len(re.findall("T",seq,re.IGNORECASE))
                    if Anum == 0 and Tnum == 0:
                        value = 0
                    else:
                        value = float((Anum-Tnum)/(Anum+Tnum))
                    sublist.append(value)     
            elif sign == '-':
                #原本是这样的，左、中、右3部分没有按负链从右到左来，现在改为负链一致从右到左(5' -> 3')
                #不对每个bin的序列反向是因为提取序列已经是反向互补了
                for l in range(int(int(line[2])+up),int(line[2]),-wdsize):
                    left = str(l-wdsize)
                    right = str(l)
                    a = str("{'chr':"+"'"+line[0]+"'"+","+"'start':"+left+","+"'stop':"+right+","+"'strand':"+"'-'"+"}")
                    L.append(a)
                for l in range(int(line[2]),int(line[1])-reside,-blockSize):
                    left = str(l-blockSize)
                    right = str(l)
                    a = str("{'chr':"+"'"+line[0]+"'"+","+"'start':"+left+","+"'stop':"+right+","+"'strand':"+"'-'"+"}")
                    M.append(a)
                for l in range(int(line[1]),int(int(line[1])-down),-wdsize):     
                    left = str(l-wdsize)
                    right = str(l)
                    a = str("{'chr':"+"'"+line[0]+"'"+","+"'start':"+left+","+"'stop':"+right+","+"'strand':"+"'-'"+"}")
                    R.append(a)
                for i in L:
                    convertI = eval(i)
                    seq = fa.sequence(eval(i), one_based=False)
                    Anum = len(re.findall("A",seq,re.IGNORECASE))
                    Tnum = len(re.findall("T",seq,re.IGNORECASE))
                    if Anum == 0 and Tnum == 0:
                        value = 0
                    else:
                        value = float((Anum-Tnum)/(Anum+Tnum))
                    sublist.append(value)
                for i in M:
                    convertI = eval(i)
                    seq = fa.sequence(eval(i), one_based=False)
                    Anum = len(re.findall("A",seq,re.IGNORECASE))
                    Tnum = len(re.findall("T",seq,re.IGNORECASE))
                    if Anum == 0 and Tnum == 0:
                        value = 0
                    else:
                        value = float((Anum-Tnum)/(Anum+Tnum))
                    sublist.append(value)
                for i in R:
                    convertI = eval(i)
                    seq = fa.sequence(eval(i), one_based=False)
                    Anum = len(re.findall("A",seq,re.IGNORECASE))
                    Tnum = len(re.findall("T",seq,re.IGNORECASE))
                    if Anum == 0 and Tnum == 0:
                        value = 0
                    else:
                        value = float((Anum-Tnum)/(Anum+Tnum))
                    sublist.append(value)
            alllist.append(sublist)      
        else:
            pass
    f.close()
    return alllist

of = open(args.outfile,"w")
outfigure=open(args.figurefile,"w")
figurelist=[]
left=int(0)
mediumA=int(args.upstream/args.wdsize)
mediumB=mediumA+args.blocknum
right=int(args.downstream/args.wdsize)+mediumA+args.blocknum

for bed in bedlist:
    Alllist = GCValue(bed,args.wdsize,args.fastafile,args.blocknum,args.upstream,args.downstream,args.strand)
    tup=np.array(Alllist)    
    tup2=tup.mean(axis=0)
    tup1=list(tup2)
    figurelist.append(tup1)
print("region","model","value",sep="\t",end="\n",file=of)
N = 0
for lis in figurelist:
    n = 0
    for l in lis:
        print(n*int(args.wdsize),legend[N],l,sep="\t",end="\n",file=of)
        n += 1
    N += 1
for line in figurelist:
    plotsmooth(line)
plt.legend(legend)
plt.ylabel("AT skew")
plt.xticks([left,mediumA,mediumB,right],["-"+str(args.upstream)+"bp","start","end",str(args.downstream)+"bp"])
plt.savefig(outfigure)     
outfigure.close()
of.close()

#print(convertI["chr"],convertI["start"],convertI["stop"],value,sep="\t",end="\n",file=GCconfile)
#GCValue(args.bedfile,args.wdsize,args.fastafile,args.blocknum,args.upstream,args.downstream,args.outfile,args.figurefile)





    
