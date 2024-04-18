#coding=utf-8

import argparse
example_text = '''example:
    python this.py 
'''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--res_blast','-b',type=str,help="",required= True,metavar='')
parser.add_argument('--seq_annotation', '-a',type= str,help="",required= True,metavar='')
parser.add_argument('--outputf', '-o',type= str,help="",required= True,metavar='')
args = parser.parse_args()

res_blast=args.res_blast
seq_annotation=args.seq_annotation
outputf=args.outputf

from collections import defaultdict

all_list = defaultdict(list)
all_seq = []

with open(res_blast, 'r') as all_id:
    for line in all_id:
        if not line.startswith("#"):
            lines = line.strip().split("\t")
            if lines[0] not in all_seq:
                all_list[lines[1]].append(line.strip())
                all_seq.append(lines[0])

with open(outputf, 'w') as outf:
    print("# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score",file = outf)
    with open(seq_annotation, 'r') as anno_f:
        for line in anno_f:
            lines = line.lstrip(">").split(" ")
            if lines[0] in all_list:
                for id in all_list[lines[0]]:
                    print(line.strip(),id,sep='\t', file = outf)
                # print(*all_list[lines[0]]:,lines[0], line.strip(),sep='\t', file = outf)
