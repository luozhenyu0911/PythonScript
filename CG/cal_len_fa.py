#coding=utf-8
from __future__ import print_function 
import sys
import argparse
example_text = '''example:
    python this.py -f xx.fa -m max_fa
'''
parser = argparse.ArgumentParser(description="To calculate the length of fasta that file have multiple lines.",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--list_file','-l',type=str,help="xx.fa file",required= True,metavar='')
parser.add_argument('--output', '-o',type= str,help="output file",required= True,metavar='')
args = parser.parse_args()

list_file=args.list_file
outputf=args.output

def get_all_fa(fa_file):
    seq = {}
    name1=''
    for line in fa_file:
        if line.startswith('>'):
            name1 = line.strip()
            seq[name1] = ''
        else:
            seq[name1] += line.strip()
    return seq

output2 = open(outputf, 'w')
with open(list_file,'r') as list_f:
    for line in list_f:
        line = line.strip()
        fname = line.split(".")[0]+".fa"
        with open(line,'r') as fa_file:
            seq_dict = get_all_fa(fa_file)
            with open(fname,'w') as output1:
                for name, seq in seq_dict.items():
                    len_fa = len(seq.strip())
                    print(name+"_"+str(len_fa), seq, sep='\n',file = output1)
                    print(name+"_"+str(len_fa), len_fa, sep='\t',file = output2)

output2.close()
