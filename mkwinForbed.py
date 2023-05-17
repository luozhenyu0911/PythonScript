#coding=utf-8
from __future__ import print_function
import argparse
example_text = '''example:
    python this.py -i input.txt  -o output.txt
'''
parser = argparse.ArgumentParser(description="  ",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--bed_file','-b',type=str,help=" ",required= True,metavar='')
parser.add_argument('--output', '-o',type= str,help=" ",required= True,metavar='')

args = parser.parse_args()

bed_file=args.bed_file
output=args.output

outputf = open(output, 'w')
with open(bed_file, 'r') as bedf:
    for line in bedf:
        lines = line.strip().split('\t')
        chr = lines[0]
        start = int(lines[1])
        end = int(lines[2])
        bc = lines[3]
        span = (end - start)//10
        for i in range(10):
            print(chr, start+span*i, start+span*(i +1), bc, "bin"+str(i+1), sep = '\t', file = outputf)
            
outputf.close()
