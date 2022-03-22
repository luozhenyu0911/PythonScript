# -*- coding: UTF-8 -*-
from __future__ import print_function
from __future__ import division
#from numba import jit
import sys, getopt
import pysam 
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
plt.switch_backend('PDF')
##@jit

def message():
    print("\n"+"Usage: python **.py --reference-point or --scale-regions [options]"+"\n")
    print("Description:"+"\n"+"\t"+"\t"+"This tool output bamfile reads midpoint."+"\n")
    print("Options:"+"\n")
    print("\t"+"-v, "+"--version"+"\t"+"show program's version number and exit")
    print("\t"+"-h, "+"--help"+"\t"+"show this help message and exit.")
    print("\t"+"-b"+"    <bamfiles>")
    print("\t"+"-w"+"    <window size>")
    print("\t"+"-o"+"    <outputfiles distance file>")
    print("\t"+"-C"+"    <chromsome name>"+"\t"+"Chromosome that you want to count."+"\n")
    print("An example usage is:"+"\n")
    print("\t"+"python **.py -b **.bam -C chromfile -w 1000 -o **bamfile.distance"+"\n")
    return

opts,args = getopt.getopt(sys.argv[1:],'-h-v-b:-C:-w:-o:-F:',['help','version'])
for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        message()
    if opt_name in ('-v','--version'):
        print("**.py 1.0 2020-06-01")
    if opt_name in ('-b'):
        bamfile = opt_value
    if opt_name in ('-C'):
        Chromfile = opt_value
    if opt_name in ('-w'):
        window = opt_value
        Wdsize=int(window)
    if opt_name in ('-F'):
        figurename = opt_value
    if opt_name in ('-o'):
        outputfile = opt_value
    if opt_name not in ('-h,--help,-v,--version,-b,-w,-o,-C'):
        message()

Option={}
for option in opts:
    Option[option[0]]=option[1]

'''
ReadInfo {query_name:['chr1',ref_start,ref_end],...}

midpoint {'1':{query_name:midpoint,
               query_name:midpoint},
          '2':{query_name:midpoint,
               query_name:midpoint},
               ...,
          '10':{query_name:midpoint,
               query_name:midpoint}}
或者
midpoint {readname:[],
          chr:[],
          midpoint:[]}
'''
#一条染色体上怎么做：给定一个window,把染色体均匀分成若干个相等的window大小，排列组合计算每一个window内所有reads中点之间的距离
#保存到字典的reads如何排序？
#字典排不了序，只能转换为pandas 构建一个行名为read name,列名为reads中点的数据框
#read name保存到一个列表，midpoint保存到一个列表，直接pd.DataFrame({'midpoint':midlist},index=readname)
#排序之后，按照分好的窗口，遍历每一个窗口，对每一个窗口下的reads排列组合统计中点距离，结果累加到pandas数据框中，行名为reads中点之间距离，列名为统计到这个距离的reads对
# def bamTodict(Input):
#     ReadInfo={}
#     for read in Input.fetch():
#         ID=read.query_name
#         if ID in ReadInfo:
#             if int(read.reference_start) < int(ReadInfo[ID][1]):
#                 ReadInfo[ID][1]=int(read.reference_start)
#             else:
#                 ReadInfo[ID][2]=int(read.reference_end)
#         else:
#             a=[]
#             a.append(str(read.reference_name))
#             a.append(int(read.reference_start))
#             a.append(int(read.reference_end))
#             ReadInfo[ID]=a
#     region={}
#     Readname=[]
#     Chr=[]
#     Midpoint=[]
#     for reads,info in ReadInfo.items():
#         readname=reads
#         chrom=info[0]
#         start=info[1]
#         end=info[2]
#         Readname.append(reads)
#         Chr.append(chrom)
#         Midpoint.append(int((start+end)/2))
#     region['Chr']=Chr
#     region['midpoint']=Midpoint
#     Df=pd.DataFrame(region,index=Readname)
#     Dfsort=Df.sort_values(by=['Chr','midpoint'],ascending=[True,True])
#     return Dfsort

# def MidpointDistanceFrequency(readdict,wdsize,chromfile):
#     DistanceDict={}
#     Chr=[]
#     with open(chromfile,'r') as f:
#         for line in f:
#             line=line.strip()
#         Chr.append(line)
#     for chr in Chr:
#         chrmidpoint=readdict.loc[readdict['Chr']==chr]
#         start=min(chrmidpoint['midpoint'])
#         end=max(chrmidpoint['midpoint'])
#         remainder=abs(end-start)%wdsize
#         endnew=(wdsize-remainder)+end
#         for i in range(start,endnew,wdsize):
#             tmpDf=chrmidpoint.loc[(chrmidpoint['midpoint']>=i) & 
#                                   (chrmidpoint['midpoint']<i+wdsize),
#                                  'midpoint']
#             readname=tmpDf.index
#             for combns in itertools.combinations(readname,2):
#                 readA=combns[0]
#                 readB=combns[1]
#                 distance=abs(tmpDf[readA]-tmpDf[readB])
#                 if str(distance) in DistanceDict:
#                     DistanceDict[str(distance)]+=1
#                 else:
#                     DistanceDict[str(distance)]=1
#     output=pd.Series(DistanceDict)
#     return output

# def MidpointDistanceFrequency(readdict,wdsize,chromfile):
#     DistanceDict={}
#     Chr=[]
#     with open(chromfile,'r') as f:
#         for line in f:
#             line=line.strip()
#         Chr.append(line)
#     for chr in Chr:
#         chrmidpoint=readdict.loc[readdict['Chr']==chr,'midpoint']
#         readname=chrmidpoint.index
#         for combns in itertools.combinations(readname,2):
#             readA=combns[0]
#             readB=combns[1]
#             distance=abs(chrmidpoint[readA]-chrmidpoint[readB])
#             if distance <= wdsize:
#                 if str(distance) in DistanceDict:
#                     DistanceDict[str(distance)]+=1
#                 else:
#                     DistanceDict[str(distance)]=1
#             else:
#                 pass
#     output=pd.Series(DistanceDict)
#     return output

def bamTodict(Input):
    # first: record every read start end and fragment
    ReadInfo={}
    for read in Input.fetch():
        ID=read.query_name
        if ID in ReadInfo:
            if int(read.reference_start) < int(ReadInfo[ID][1]):
                ReadInfo[ID][1]=int(read.reference_start)
            else:
                ReadInfo[ID][2]=int(read.reference_end)
        else:
            a=[]
            a.append(str(read.reference_name))
            a.append(int(read.reference_start))
            a.append(int(read.reference_end))
            ReadInfo[ID]=a
    # second: count midpoint per bp
    region={}
    for value in ReadInfo.values():
        chrom=value[0]
        start=value[1]
        end=value[2]
        centerPoint=int((start+end)/2)
        if chrom in region:
            pass
        else:
            region[chrom]={}
        if str(centerPoint) in region[chrom]:
            region[chrom][str(centerPoint)]+=1
        else:
            region[chrom][str(centerPoint)]=1
    return region
    
def MidpointDistanceFrequency(readdict,wdsize,chromfile):
    DistanceDict={}
    Chr=[]
    with open(chromfile,'r') as f:
        for line in f:
            line=line.strip()
        Chr.append(line)
    for chr in Chr:
        CoverageDict=readdict[chr]
        keys=CoverageDict.keys()#错在单独遍历每个位点
        for combns in itertools.combinations(keys,2):
            positionA=combns[0]
            positionB=combns[1]
            distance=abs(int(positionA)-int(positionB))
            if distance <= wdsize:
                DistanceCount=CoverageDict[positionA]*CoverageDict[positionB]
                if str(distance) in DistanceDict:
                    DistanceDict[str(distance)]+=DistanceCount
                else:
                    DistanceDict[str(distance)]=0
                    DistanceDict[str(distance)]+=DistanceCount
            else:
                pass
    output=pd.Series(DistanceDict)
    return output

if __name__ == '__main__':
    bam=pysam.AlignmentFile(bamfile,'rb')
    MidDataFrame=bamTodict(bam)
    distance=MidpointDistanceFrequency(MidDataFrame,Wdsize,Chromfile)
    distance.to_csv(outputfile,sep='\t')

# def plotsmooth(Values):
#     n=len(Values)
#     num=range(0,n)
#     num=np.array(num)
#     xnew=np.linspace(num.min(),num.max(),5*n)
#     func=interp1d(num,Values,kind="cubic")
#     ynew=func(xnew)
#     plt.plot(xnew,ynew)


