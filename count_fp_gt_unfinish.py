#coding=utf-8
from __future__ import print_function 
import argparse
import gzip
example_text = '''example:
    python this.py -i fp.vcf/fp.vcf.gz -b bench_mark.vcf/bench_mark.vcf.gz -o fp_in_bench.txt -s statistc.txt
'''
parser = argparse.ArgumentParser(description="The script is to divide the expression matrix into high, middle, and low files.",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--input','-i',type=str,help="fp.vcf",required= True,metavar='')
parser.add_argument('--bench_mark', '-b',type= str,help="bench_mark.vcf",required= True,metavar='')
parser.add_argument('--output', '-o',type= str,help="fp_in_bench.txt",required= True,metavar='')
parser.add_argument('--statistc', '-s',type= str,help="statistc.txt",required= True,metavar='')

args = parser.parse_args()

#记得习惯性空一格
input=args.input
bench_mark=args.bench_mark
output=args.output
statistc=args.statistc


def region_GT(vcffile):
    region={}
    for line in vcffile:
        if line.strip().startswith('#'):
            pass
        else:
            lines = line.strip().split('\t')
            chr=lines[0].strip()
            site=lines[1].strip()
            gt=lines[8].strip()
            gttype=lines[9].strip().split(':')[0]
            chr_site=str(chr)+"_"+str(site)
            # region.setdefault(chr_site,[]).append(site,gt)
            region[chr_site]=gttype
    return region

def statistic_gt(region):
    gt_dic = {}
    for key, value in region.items():
        gt_dic[value] = gt_dic.get(value,0)+1
    return gt_dic

if str(input).endswith('gz'):
    fp_file= gzip.open(input,'rb')
else:
    fp_file= open(input,'r')
stat =open(statistc,'w')

fp_region=region_GT(fp_file)
gt_dic=statistic_gt(fp_region)
print("# summary genome type in fp.vcf",file=stat)

for k,v in gt_dic.items():
    print(k,v,sep='\t',file=stat)

fp2benchmak = open(output,'w')
print("chr_site","FP_GT","Bench_MK_GT",sep='\t',file=fp2benchmak)

if str(bench_file).endswith('gz'):
    bench_file = gzip.open(bench_mark,'rb')
else:
    bench_file = open(bench_mark,'r')


fp_GT_dict={}
bench_in_fp={}
for line in bench_file:
    if line.strip().startswith('#'):
        pass
    else:
        lines = line.strip().split('\t')
        chr=lines[0].strip()
        site=lines[1].strip()
        gt=lines[8].strip()
        gttype=lines[9].strip().split(':')[0]
        chr_site=str(chr)+"_"+str(site)
        if chr_site in fp_region:
            print(chr_site,fp_region[chr_site],gttype,sep='\t', file=fp2benchmak)
            fp_GT_dict.setdefault(fp_region[chr_site],[]).append(gttype)
            bench_in_fp[chr_site]='gttype'
        else:
            pass

for k in fp_region:
    if k not in bench_in_fp:
        print(k,fp_region[k],'nothing',sep='\t', file=fp2benchmak)
        fp_GT_dict.setdefault(fp_region[k],[]).append('nothing')
    else:
        pass

# print(fp_GT_dict)
count_dict={}
for k,v in fp_GT_dict.items():
    count_dict[k]={}
    for gt_type in v:
        if gt_type in count_dict[k]:
            count_dict[k][gt_type]+=1
        else:
            count_dict[k][gt_type]=1
#         print(count_dict)
print()
print("# change of the GT in fp.vcf in benchmark.vcf",file=stat)
for k,v in count_dict.items():
    for i,j in v.items():
        print(k,i,j,sep='\t',file=stat)
