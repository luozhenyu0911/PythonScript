#coding=utf-8
from __future__ import print_function 
from __future__ import division
import sys
import argparse
example_text = '''example:
    python this.py 
'''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

group=parser.add_mutually_exclusive_group(required= True)
group.add_argument('--part_count','-p',type=int,help="How many parts to divide the file into",metavar='')
group.add_argument('--line_count','-l',type=int,help='Split the file by how many lines',metavar='')

parser.add_argument('--inputfile','-i',type=str,help="input file",required= True,metavar='')
parser.add_argument('--prefix', '-o',type= str,help="prefix+part.txt",required= True,metavar='')
args = parser.parse_args()

part_count=args.part_count
line_count=args.line_count
inputfile=args.inputfile
prefix=args.prefix

with open(inputfile, 'r') as inputf:
    infile = inputf.readlines()
    lines_count = len(infile)
    
    if part_count:
        start = 0
        size = lines_count//part_count
        end = size
        for i in range(part_count):
            outputf = prefix + str(i+1) +".txt"
            with open(outputf, 'w') as outf:
                outf.write(''.join(infile[start:end]))
            start += size
            end += size
        with open(prefix + str(part_count+1) +".txt", 'w') as outf:
    #         print(start,end)
            end = start+lines_count%part_count
    #         print(lines_count%part_count,start,end)
            outf.write(''.join(infile[start:end]))
    #     os.system("tail -n " + str(nrow%size) + " " + inpath + ">" prefix + str(part_count+1) +".id")
    elif line_count:
        start = 0
        part_count = lines_count//line_count
        size = line_count
        end = size
        for i in range(part_count):
            outputf = prefix + str(i+1) +".txt"
            with open(outputf, 'w') as outf:
                outf.write(''.join(infile[start:end]))
            start += size
            end += size
        with open(prefix + str(part_count+1) +".txt", 'w') as outf:
            end = start+lines_count%line_count
            outf.write(''.join(infile[start:end]))