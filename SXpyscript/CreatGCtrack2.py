# -*- coding: UTF-8 -*-
from __future__ import print_function
from __future__ import division
from scipy.signal import savgol_filter
from pyfasta import Fasta
import pyBigWig
import numpy as np
import argparse

example_text = '''example:
 2020.10.16
 python *py -b *.fa -o *.bw
'''
parser = argparse.ArgumentParser(description="This tool creat genome G/C content bigwig from a fasta.", 
                                 formatter_class= argparse.RawTextHelpFormatter,
                                 prog='base_maker',
                                 epilog=example_text)
parser.add_argument('--Fasta','-F',type= str,help="<genome fasta file>",required= True)
parser.add_argument('--binSize','-b',type= int,help="<bin size>",required= True)
#parser.add_argument('--smoothLength',type= int,help="<The smooth length defines a window, larger than the binSize, to average the number of reads.>",required= True)
parser.add_argument('--output','-o',type=str,help="<output bigwig file>",required= True)
args = parser.parse_args()

fa=args.Fasta
bs=args.binSize
bw=args.output

def computeGenomeGC(fasta,binsize,outbigwig):
    bigwig=pyBigWig.open(outbigwig,'w') #创建一个bw
    fa=Fasta(fasta,key_fn=lambda key: key.split()[0])
    #向bw中添加header，必须要一次性添加才不报错
    header=[(key,len(fa[key])) for key in fa.keys()]
    bigwig.addHeader(header)
    for key in fa.keys():
        chroms=[]
        starts=[]
        ends=[]
        values=[]
        seq=fa[key]
        chrlen=len(fa[key])
        for start in range(0,chrlen,binsize):
            if start+binsize<chrlen:
                chroms.append(key)
                starts.append(start)
                ends.append(start+binsize)
                bin=seq[start:start+binsize].upper()
                G=bin.count('G')
                C=bin.count('C')
                values.append((G+C)/len(bin))
            else:
                chroms.append(key)
                starts.append(start)
                ends.append(chrlen)
                bin=seq[start:chrlen].upper()
                G=bin.count('G')
                C=bin.count('C')
                values.append((G+C)/len(bin))
        #Savitzky-Golay卷积平滑算法，smooth reads
        values_fit = savgol_filter(values, 7, 5)
        bigwig.addEntries(chroms, starts, ends=ends, values=values_fit)
    bigwig.close()

if __name__ == '__main__':
    computeGenomeGC(fa,bs,bw)
    
"""
scipy.signal.savgol_filter(x, window_length, polyorder)

x为要滤波的信号
window_length即窗口长度
取值为奇数且不能超过len(x)。它越大，则平滑效果越明显；越小，则更贴近原始曲线。
polyorder为多项式拟合的阶数。
它越小，则平滑效果越明显；越大，则更贴近原始曲线。
"""