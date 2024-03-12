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

parser.add_argument('--input','-i',type=str,help=" ",required= True,metavar='')
parser.add_argument('--output', '-o',type= str,help=" ",required= True,metavar='')

args = parser.parse_args()

input=args.input
output=args.output

with open(output, 'w') as outputf:
    with open(input, 'r') as inputf:
        for line in inputf:
            if not line.startswith('Barcode'):
                lines = line.strip().split('\t')
                res=lines[3].strip('[')
                res=res.strip(']')
                res=res.split(',')
                for id in list(res):
                    print(id.strip().strip("'"), file = outputf)
