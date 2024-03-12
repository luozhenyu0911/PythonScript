#coding=utf-8
from __future__ import print_function 
from __future__ import division
import sys
import pandas as pd
import argparse
example_text = '''example:
    python this.py -e D4merge_TPM.txt -g ~/index/Sitalica_312_v2_gene6col.bed -p D4_TPM
'''
parser = argparse.ArgumentParser(description="The script is to divide the expression matrix into high, middle, and low files.",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                          #格式化显示方式(描述程序用途的字符串)
                                epilog=example_text)

parser.add_argument('--expression_matrix','-e',type=str,help="Expression_matrix",required= True,metavar='')
parser.add_argument('--genome_bed_file', '-g',type= str,help="A six-column genome bed file is needed ",required= True,metavar='')
parser.add_argument('--prefix', '-p',type= str,help="the prefix of output file",required= True,metavar='')
args = parser.parse_args()

expr_=args.expression_matrix
genome_bed=args.genome_bed_file
prefix=args.prefix

f1=pd.read_csv(expr_,sep='\t',names=["Gene ID","TPM"])
# f1.sort_values('TPM',ascending=False,inplace=True)
# len_f1=f1.shape[0]
# f1.index = range(len(f1))  
#df.reset_index(drop=True)  重建索引，并删除原来的index（不然原来的index会变成新的一列） 
# remainder = int(len_f1%3)
# num = int(len_f1//3)
# f1h = pd.DataFrame(f1.loc[:num-1,'Gene ID']) 
# f1m = pd.DataFrame(f1.loc[num:2*num-1,'Gene ID']) 
# f1l = pd.DataFrame(f1.loc[2*num:,'Gene ID']) 
f1h = pd.DataFrame(f1[f1['TPM']>10]['Gene ID']) 
f1m = pd.DataFrame(f1[(f1['TPM']>=1) & (f1['TPM']<=10) ]['Gene ID']) 
f1l = pd.DataFrame(f1[f1['TPM']<1]['Gene ID']) 
f2=pd.read_csv(genome_bed,sep='\t',names=['chr', 'start', 'end', 'Gene ID','score','strand'])

def get_bed(genefile):
    name=pd.merge(genefile,f2,left_on='Gene ID', right_on='Gene ID', how="inner")
    order = ['chr', 'start', 'end', 'Gene ID','score','strand']
    name = name[order]
    return name
    
f1hbed=get_bed(f1h)
f1mbed=get_bed(f1m)
f1lbed=get_bed(f1l)
f1hbed.to_csv(prefix+'_high.bed',sep='\t',header=None,index=None)
f1mbed.to_csv(prefix+'_mid.bed',sep='\t',header=None,index=None)
f1lbed.to_csv(prefix+'_low.bed',sep='\t',header=None,index=None)
