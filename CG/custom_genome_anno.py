#coding=utf-8
from __future__ import print_function
from __future__ import division
import argparse
import re
example_text = '''example:
    python this.py -f xx.fa -m max_fa
'''
parser = argparse.ArgumentParser(description="To get custom genome annotation file gff.",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--input_gff','-i',type=str,help="the raw gff file ",required= True,metavar='')
parser.add_argument('--output_gff', '-o',type= str,help="the custom gff file",required= True,metavar='')
args = parser.parse_args()

input_gff=args.input_gff
output_gff=args.output_gff

outputf = open(output_gff, 'w')

with open(input_gff, 'r') as gff:
    pattern = re.compile('ID=gene-(.*?);')
    for line in gff:
        if line.startswith("#"):
            print(line.strip(), file = outputf)
        elif "region" in line.strip().split('\t')[2]:
            pass
        else:
            lines = line.strip().split('\t')
            if lines[8].startswith("ID=gene-") and "gene" in lines[2]:
                start = lines[3]
                end = lines[4]
                chr = re.findall(pattern, line)[0]
                lines[3] = int(lines[3]) - int(start) + 1
                lines[4] = int(lines[4]) - int(start) + 1
                lines[0] = chr
                print(*lines, sep='\t', file = outputf)
            elif int(lines[4]) > int(end) or int(lines[3]) < int(start):
                pass
            else:
                lines[3] = int(lines[3]) - int(start) + 1
                lines[4] = int(lines[4]) - int(start) + 1
                lines[0] = chr
                print(*lines, sep='\t', file = outputf)
                
outputf.close()
