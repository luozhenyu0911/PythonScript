#coding=utf-8
from __future__ import print_function 
from __future__ import division
import sys
import argparse
example_text = '''example:
    python this.py ..
'''
parser = argparse.ArgumentParser(description=".",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--inputsum','-s',type=str,help="Expression_matrix",required= True,metavar='')#  注意点：'--expression_matrix'，引号中间不能有空格
parser.add_argument('--inputpart','-p',type=str,help="Expression_matrix",required= True,metavar='')#  注意点：'--expression_matrix'，引号中间不能有空格
parser.add_argument('--output1', '-o',type= str,help="A six-column genome bed file is needed ",required= True,metavar='')
args = parser.parse_args()

#记得习惯性空一格
input1=args.inputsum
input2=args.inputpart
output1=args.output1

region={}
inputf1=open(input1,'r')
inputf2=open(input2,'r')
outputf=open(output1,'w')
for line in inputf1:
    if line.startswith('#'):
        pass
    else:
        lines=line.strip().split('\t')
        if lines[0] in region:
            region[lines[0]].append(lines[1])
        else:
            region[lines[0]]=[]
        
for line in inputf2:
    if line.startswith('#'):
        pass
    else:
        lines=line.strip().split('\t')
        if lines[1] in region[lines[0]]:
            print(line.strip(),file=outputf)
inputf1.close()
inputf2.close()
outputf.close()



