#coding=utf-8
from __future__ import print_function
import argparse
import pysam
example_text = '''example:
    python this.py -i input.txt  -o output.txt
'''
parser = argparse.ArgumentParser(description="  ",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--input','-i',type=str,help=" ",required= True,metavar='')
parser.add_argument('--bam','-b',type=str,help=" ",required= True,metavar='')
parser.add_argument('--output', '-o',type= str,help=" ",required= True,metavar='')

args = parser.parse_args()

input=args.input
output=args.output
bam=args.bam

def list_all(file1):
    all_list = []
    for line in file1:
        line = line.strip()
        if line not in all_list:
            all_list.append(line)
    return all_list

bamfile = pysam.AlignmentFile(bam, "rb")
outputf = pysam.AlignmentFile(output, "wb", template=bamfile)

with open(input, 'r') as inputf:
    all_list = list_all(inputf)
    
for read in bamfile.fetch():
    bx = read.get_tag('BX')
    if bx in all_list:
        outputf.write(read)

bamfile.close()
outputf.close()
