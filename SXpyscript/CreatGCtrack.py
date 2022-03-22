# -*- coding: UTF-8 -*-
from __future__ import print_function
from __future__ import division
from pyfasta import Fasta
import pyBigWig
import numpy as np
import argparse

example_text = '''example:
 2020.08.09
 python *py -b *.fa -o *.bw
'''
parser = argparse.ArgumentParser(description="This tool creat genome G/C content bigwig from a fasta.", 
                                 formatter_class= argparse.RawTextHelpFormatter,
                                 prog='base_maker',
                                 epilog=example_text)
parser.add_argument('--Fasta','-F',type= str,help="<genome fasta file>",required= True)
parser.add_argument('--binSize','-b',type= int,help="<bin size>",required= True)
parser.add_argument('--output','-o',type=str,help="<output bigwig file>",required= True)
args = parser.parse_args()

fa=args.Fasta
bs=args.binSize
bw=args.output

def computeGenomeGC(fasta,binsize,outbigwig):
    chroms=[]
    starts=[]
    ends=[]
    values=[]
    bigwig=pyBigWig.open(outbigwig,'w') #创建一个bw
    fa=Fasta(fasta,key_fn=lambda key: key.split()[0])
    #向bw中添加header，必须要一次性添加才不报错
    header=[(key,len(fa[key])) for key in fa.keys()]
    bigwig.addHeader(header)
    for key in fa.keys():
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
    bigwig.addEntries(chroms, starts, ends=ends, values=values)
    bigwig.close()

if __name__ == '__main__':
    computeGenomeGC(fa,bs,bw)