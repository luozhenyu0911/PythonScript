#coding=utf-8
from __future__ import print_function 
import sys
import argparse
example_text = '''example:
    python this.py -v xx.vcf -o output -q 30
'''
parser = argparse.ArgumentParser(description="To get the span between two mutation sites in xx.vcf.",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--vcf_file','-v',type=str,help="vcf_file",required= True,metavar='')
parser.add_argument('--snp_quality','-q',type=str,help="snp_quality",required= False,metavar='')
parser.add_argument('--snp_depth','-d',type=str,help="snp_depth",required= False,metavar='')
parser.add_argument('--output', '-o',type= str,help="output file",required= True,metavar='')
args = parser.parse_args()

vcf_file=args.vcf_file
output=args.output
    
def get_distance_qual(vcf_file, qual):
    chr1 = "chr"
    pos1 = 0
    distance_vcf={}
    for line in vcf_file:
        if line.startswith("#"):
            pass
        elif float(line.strip().split('\t')[5]) < qual:
            pass
        else:
            lines = line.strip().split('\t')
            chr2 = lines[0]
            pos2 = int(lines[1])
            if chr2 == chr1:
                distance_vcf[chr2][chr2+":"+str(pos1)+"-"+str(pos2)]=int(pos2-pos1)
                pos1 = pos2
            else:
                distance_vcf[chr2]={}
                chr1 = chr2
                pos1 = pos2
    return distance_vcf

def get_distance_qual_depth(vcf_file, qual, depth):
    chr1 = "chr"
    pos1 = 0
    distance_vcf={}
    for line in vcf_file:
        if line.startswith("#"):
            pass
        elif float(line.strip().split('\t')[5]) <= qual and int(line.strip().split('\t')[9].split(":")[2]) <= depth:
            pass
        else:
            lines = line.strip().split('\t')
            chr2 = lines[0]
            pos2 = int(lines[1])
            if chr2 == chr1:
                distance_vcf[chr2][chr2+":"+str(pos1)+"-"+str(pos2)]=int(pos2-pos1)
                pos1 = pos2
            else:
                distance_vcf[chr2]={}
                chr1 = chr2
                pos1 = pos2
    return distance_vcf

def get_distance(vcf_file):
    chr1 = "chr"
    pos1 = 0
    distance_vcf={}
    for line in vcf_file:
        if line.startswith("#"):
            pass
        else:
            lines = line.strip().split('\t')
            chr2 = lines[0]
            pos2 = int(lines[1])
            if chr2 == chr1:
                distance_vcf[chr2][chr2+":"+str(pos1)+"-"+str(pos2)]=int(pos2-pos1)
                pos1 = pos2
            else:
                distance_vcf[chr2]={}
                pos1 = pos2
                chr1 = chr2
    return distance_vcf


outf = open(output,'w')
print("Chr",'position',"span",sep='\t', file=outf)
with open(vcf_file,'r') as vcff:
    if args.snp_quality and args.snp_depth :
        qual= int(args.snp_quality)
        depth = int(args.snp_depth)
        distance_vcf = get_distance_qual_depth(vcff, qual, depth)
    elif args.snp_quality:
        qual= int(args.snp_quality)
        distance_vcf = get_distance_qual(vcff, qual)
    else:
        distance_vcf = get_distance(vcff)
    for chr, pos in distance_vcf.items():
        for region, dis in pos.items():
            print(chr, region,dis, sep='\t',file=outf)
outf.close()
