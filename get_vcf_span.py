#coding=utf-8
from __future__ import print_function 
import sys
import argparse
example_text = '''example:
    python this.py -v xx.vcf -o output
'''
parser = argparse.ArgumentParser(description="To get the span between two mutation sites in xx.vcf.",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--vcf_file','-v',type=str,help="vcf_file",required= True,metavar='')
parser.add_argument('--output', '-o',type= str,help="output file",required= True,metavar='')
args = parser.parse_args()

vcf_file=args.vcf_file
output=args.output
def get_distance(vcf_file):
    chr1 = "chr"
    pos1 = 0
    distance_vcf={}
    for line in vcf_file:
        if line.startswith("#"):
            pass
        else:
            lines = line.strip().split('\t')
            chr2 = lines[0]
            pos2 = int(lines[1])
            if chr2 == chr1:
                distance_vcf[chr2].append(pos2-pos1)
                pos1 = pos2
            else:
                distance_vcf[chr2]=[]
                chr1 = chr2
    return distance_vcf

outf = open(output,'w')
print("Chr","span",sep='\t', file=outf)
with open(vcf_file,'r') as vcff:
    distance_vcf = get_distance(vcff)
    for chr, pos in distance_vcf.items():
        for dis in pos:
            print(chr, dis, sep='\t', file=outf)
outf.close()
