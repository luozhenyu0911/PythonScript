#coding=utf-8
from __future__ import print_function 
from __future__ import division
from collections import defaultdict
import gzip
import argparse
example_text = '''example:
   python seq_venn_cov.py -v 00114031182M22BFF2.V.allmapped.bed -b 00114031182M22BFF2.readnV.human.bed -g viral.1.1.genomic.bed -c viral_genomic_V.all.per-base.bed.gz -p 00114031182M22BFF2
'''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--viral_bed','-v',type=str,help="the viral bed file made from bam file with bam2bed tool",required= True,metavar='')
parser.add_argument('--viral_cov_file', '-c',type= str,help="the viral coverage file made from bam file with mosdepth tool",required= True,metavar='')
parser.add_argument('--genome_bed_file', '-g',type= str,help="the genome bed file made from reference genome",required= True,metavar='')
parser.add_argument('--host_bed_file', '-b',type= str,help="the host bed file made from host bam file with bam2bed tool",required= True,metavar='')
parser.add_argument('--prefix', '-p',type= str,help="the prefix of output file",required= True,metavar='')
args = parser.parse_args()



def gene_len(bedf):
    """
    make a dictionary of gene lengths from a ref bed file
    """
    gene_len_dict = defaultdict(int)
    with open(bedf, 'r') as f:
        for line in f:
            gene, start, end, *_ = line.strip().split('\t')
            gene_len_dict[gene] = int(end) - int(start)
    return gene_len_dict

def act_len(bedgz):
    """
    make a dictionary of actual region lengths from a bed.gz file
    """
    with gzip.open(bedgz, 'rt') as f:
        act_len_dict = defaultdict()
        for line in f:
            if int(line.strip().split('\t')[3]) > 0:
                chrom, start, end, count  = line.strip().split('\t')
                act_len_dict[chrom] = act_len_dict.get(chrom, 0) +  (int(end) - int(start))
    return act_len_dict

def gene_cov(expect_dict, actual_dict):
    """
    calculate the coverage of each gene
    """
    cov_dict = defaultdict(float)
    for gene in actual_dict:
        cov_dict[gene] = actual_dict[gene] / expect_dict[gene]
    return cov_dict

def seq_comm(viral_bedf, host_bedf):
    viral_seq = set()
    host_seq = set()
    with open(viral_bedf, 'r') as f1, open(host_bedf, 'r') as f2:
        for line in f1:
            viral_seq.add(line.strip().split('\t')[3])
        for line in f2:
            host_seq.add(line.strip().split('\t')[3])
        shared_seq = viral_seq.intersection(host_seq)
    return shared_seq

if __name__ == '__main__':
    viral_bedf = args.viral_bed
    viral_cov_file = args.viral_cov_file
    genome_bed_file = args.genome_bed_file
    host_bed_file = args.host_bed_file
    prefix = args.prefix
    # make a dictionary of gene lengths from the ref bed file
    gene_len_dict = gene_len(genome_bed_file)

    # make a dictionary of actual region lengths from the viral bed.gz file
    act_len_dict = act_len(viral_cov_file)

    # make a dictionary of viral coverage from the viral cov file
    cov_dict = gene_cov(gene_len_dict, act_len_dict)

    with open(prefix + '_all_viral_cov.txt', 'w') as outf:
        for virus in cov_dict:
            print(virus, cov_dict[virus], file=outf, sep='\t')

    shared_seq = seq_comm(viral_bedf, host_bed_file)
    with open(prefix + '_shared_virus_cov.txt', 'w') as outf:
        with open(viral_bedf, 'r') as f:
            virus_seq_dict = defaultdict(list)
            for line in f:
                seq_name = line.strip().split('\t')[3]
                if seq_name in shared_seq:
                    virus = line.strip().split('\t')[0]
                    virus_seq_dict[virus].append(seq_name)
            for virus in virus_seq_dict:
                # for seq_name in virus_seq_dict[virus]:
                cov = cov_dict[virus]
                print(virus, cov,*virus_seq_dict[virus], file=outf, sep='\t')
    with open(prefix + '_shared_human_bed.txt', 'w') as outf:

        with open(host_bed_file, 'r') as f:
            for line in f:
                seq_name = line.strip().split('\t')[3]
                if seq_name in shared_seq:
                    print(line.strip(), file=outf)
