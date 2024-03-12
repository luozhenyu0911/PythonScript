#coding=utf-8

import argparse
example_text = '''example:
    python this.py 
'''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--annotation','-a',type=str,help="",required= True,metavar='')
parser.add_argument('--idlist', '-i',type= str,help="",required= True,metavar='')
parser.add_argument('--prefix', '-p',type= str,help="",required= True,metavar='')
args = parser.parse_args()

annotation=args.annotation
idlist=args.idlist
prefix=args.prefix

from collections import defaultdict

all_list = defaultdict(list)

with open(idlist, 'r') as all_id:
    for line in all_id:
        lines = line.strip().split("\t")
        all_list[lines[1]].append(lines[0])

with open(prefix, 'w') as outf:
    with open(annotation, 'r') as anno_f:
        for line in anno_f:
            lines = line.lstrip(">").split(" ")
            if lines[0] in all_list:
                for id in all_list[lines[0]]:
                    print(id,lines[0], line.strip(),sep='\t', file = outf)
                # print(*all_list[lines[0]]:,lines[0], line.strip(),sep='\t', file = outf)


































