# -*- coding: UTF-8 -*-
from __future__ import print_function
from __future__ import division
from pyfasta import Fasta
import numpy as np
import argparse

example_text = '''example:
 2020.07.22
 python *py -b *.fa -o *.dinucleotide
'''
parser = argparse.ArgumentParser(description="This tool compute each dinucleotide frequency of fasta.", 
                                 formatter_class= argparse.RawTextHelpFormatter,
                                 prog='base_maker',
                                 epilog=example_text)
parser.add_argument('--Fasta','-F',type= str,help="<genome fasta files>",required= True)
parser.add_argument('--outfile','-o',type=str,help="<output dinucleotide frequency>",required= True)
args = parser.parse_args()

fasta=args.Fasta
outf=args.outfile

Dinuc=['AA','AT','AC','AG',
       'TA','TT','TC','TG',
       'CA','CT','CC','CG',
       'GA','GT','GG','GC']

Pattern=[['AA'],['AT'],['AC'],['AG'],
         ['TA'],['TT'],['TC'],['TG'],
         ['CA'],['CT'],['CC'],['CG'],
         ['GA'],['GT'],['GG'],['GC']]

def GenomeFrequency(Fa,Pattern):
    gsize=0
    for i in Fa.keys():
        gsize+=len(str(Fa[i]))
    Count=0
    for pattern in Pattern:
        countlist=[str(Fa[chr]).count(pattern) for chr in Fa.keys()]
        countarray=np.array(countlist)
        Sum=countarray.sum()
        Count+=Sum
    frequency=Count/gsize
    return(frequency)

fa=Fasta(fasta,key_fn=lambda key: key.split()[0])
out=open(outf,'w')
for dinuc,pattern in zip(Dinuc,Pattern):
    genomefrequency=GenomeFrequency(fa,pattern)
    print(dinuc,genomefrequency,sep='\t',file=out)
out.close()