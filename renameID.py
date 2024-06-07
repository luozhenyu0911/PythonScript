#coding=utf-8
from __future__ import print_function
import argparse
example_text = '''example:
    python this.py 
'''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--sample_id','-s',type=str,help="id",required= True,metavar='')
parser.add_argument('--output', '-o',type= str,help="",required= True,metavar='')
parser.add_argument('--input', '-i',type= str,help="",required= True,metavar='')
parser.add_argument('--depth', '-d',type= int,help="",required= True,metavar='')
args = parser.parse_args()

sample_id=args.sample_id
output=args.output
input=args.input
depth=args.depth


with open(output, 'w') as outf:
    with open(input, 'r') as configf:
        for line in configf:
            if line.strip(" ").startswith("id"):
                print('    id: "'+sample_id+'"',file = outf)
            elif line.strip(" ").startswith("depth"):
                print('    depth: '+str(depth),file = outf)
            else:
                outf.write(line)
