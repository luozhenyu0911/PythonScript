#coding=utf-8
from __future__ import print_function 
from __future__ import division
import sys
import linecache
import argparse
example_text = '''example:
    python this.py -b bc_list.txt -p work_path -o output.txt
'''
parser = argparse.ArgumentParser(description="  ",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--bc_list','-b',type=str,help=" ",required= True,metavar='')#  注意点：'--expression_matrix'，引号中间不能有空格
parser.add_argument('--output', '-o',type= str,help=" ",required= True,metavar='')
parser.add_argument('--work_path', '-p',type= str,help=" ",required= True,metavar='')

args = parser.parse_args()

#记得习惯性空一格
bc_list=args.bc_list
work_path=args.work_path
output=args.output

def get_tatal_len(filename):
    #filenames=filename.strip().split('/')[-3].split('_')[1]
    filenames=filename.strip().split('/')[-3]
    text = linecache.getline(filename, 18)
    texts = int(text.strip().split(' ')[-1])
    tmp2=linecache.getline(filename, 28)
    Misassembled=float(tmp2.strip().split(' ')[-1])
    tmp = linecache.getline(filename, 38)
    mismatches = float(tmp.strip().split(' ')[-1])
    return filenames,texts,Misassembled,mismatches

bc_list=open(bc_list,'r')
output=open(output,'w')
print("bc_id","fragment_len","denovo_len",'ratio',"Misassembled contigs length",'mismatches per 100 kbp',"Nreads",sep='\t',file=output)

for line in bc_list:
#     print(line)
    lines=line.strip().split('\t')
    bc_id=lines[0]
    actual_len= int(lines[4])-int(lines[3])
    #filename= work_path+"/frag_"+bc_id+"/quast_contig/report.txt"
    filename= work_path+"/"+bc_id+"/quast_contig/report.txt"
    tmp_id,denovo_len,Misassembled,mismatches=get_tatal_len(filename)
    if denovo_len >= actual_len:
        # print(bc_id,actual_len,denovo_len,"1",*lines,sep='\t')
        print(bc_id,actual_len,denovo_len,"1",Misassembled,mismatches,args.bc_list,sep='\t',file=output)
    else:
        rate=denovo_len/actual_len
        # print(bc_id,actual_len,denovo_len,rate,*lines,sep='\t')
        print(bc_id,actual_len,denovo_len,rate,Misassembled,mismatches, args.bc_list,sep='\t',file=output)


bc_list.close()
output.close()

