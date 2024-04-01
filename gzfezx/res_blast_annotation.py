#coding=utf-8
from __future__ import print_function 
from collections import defaultdict
from collections import Counter
import re
import pickle
import sys
import argparse
example_text = '''example:
    python this.py 
'''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--res_blast','-r',type=str,help=" ",metavar='')
parser.add_argument('--virus_anno','-v',type=str,help=' ',metavar='')
parser.add_argument('--evalue','-e',type=float,help=' ',metavar='')
parser.add_argument('--reads_count','-c',type=int,help=' ',metavar='')

parser.add_argument('--prefix', '-p',type= str,help="prefix+part.txt",required= True,metavar='')
args = parser.parse_args()

res_blast=args.res_blast
virus_anno=args.virus_anno
evalue=args.evalue
prefix=args.prefix
reads_count=args.reads_count


# select the blast result and make dictionary

def name_dict(inputf, evalue = 1e-5):
    adict = defaultdict(list)
    for line in inputf:
        if line.strip().startswith('#'):
            pass
        else:
            lines = line.strip().split('\t')
            if float(lines[-2]) < evalue:
                if lines[0] not in adict[lines[1]]:
                    adict[lines[1]].append(lines[0])
    return adict

def annotatoinDict(inputf):
    import pickle
    geneid_dict = {}
    with open(inputf, 'r') as inputf:
        for line in inputf:
            if line.strip().startswith(">"):
                geneid = re.match(r'>(\w+.\d) (.+)', line.strip()).group(1)
                annotation = re.match(r'>(\w+.\d) (.+)', line.strip()).group(2)
                geneid_dict[geneid]=annotation
    with open(prefix+'.geneid_anno.pickle','wb') as f1:
        pickle.dump(geneid_dict, f1)
    return geneid_dict
    

sel_gene = open(prefix+'.geneid_anno.txt', 'w')
stats = open(prefix+'.geneid_anno.stats.txt', 'w')
print("Ref id", "seq.count","description",sep='\t', file =sel_gene)
print('filtered(y/n)',"total virus", "total seq",sep='\t', file =stats)

seq_count = 0
seq_count_filter = 0
with open(res_blast,'r') as inputf:
    if evalue:
        adict = name_dict(inputf,evalue) # dict {virus Ref id: [seq1,seq2...]}
    else:
        adict = name_dict(inputf)
    virus_count = len(adict)
    geneid_dict = annotatoinDict(virus_anno) # dict {virus Ref id : annotation}
    if reads_count:
        with open(prefix+".G"+str(reads_count)+'.geneid_anno.txt', 'a') as outf:、
            print("Ref id", "seq.count","description",sep='\t', file =outf)
            i = 0
            for k in adict:
                print(k, len(adict[k]),geneid_dict[k],sep='\t', file =sel_gene)
                seq_count += len(adict[k])
                if len(adict[k]) > reads_count:
                    print(k, len(adict[k]),geneid_dict[k],sep='\t', file =outf)
                    seq_count_filter += len(adict[k])
                    i+=1
        print('unfiltered',virus_count, seq_count,sep='\t', file =stats)
        print('filtered',i, seq_count_filter,sep='\t', file =stats)
    else:
        for k in adict:
            print(k, len(adict[k]),geneid_dict[k],sep='\t', file =sel_gene)
            seq_count += len(adict[k])
            print('unfiltered',virus_count, seq_count,sep='\t', file =stats)
    
