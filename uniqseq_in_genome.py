# -*- coding: utf-8 -*-
from __future__ import print_function 
import subprocess
import argparse
# import itertools as it
import pysam
example_text = '''example:
    python3 this.py -l 25 -r hs38.fa -b window_5k.bed -l 25 -s 200 -o output --tools_bowtie
'''
parser = argparse.ArgumentParser(description="To get the unique sequence in whole genome.",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)
parser.add_argument('--reference', '-r',type= str,help="input the reference genome.fa",required= True,metavar='')
parser.add_argument('--unique_length','-l',type= int,help="set the unique sequence length/size bp",default= 25,metavar='')
parser.add_argument('--search_range','-s',type= int,help="the search range bp in 5k region chunk's left/right", default= 200, metavar='')
parser.add_argument('--bed', '-b',type= str,help="the split window bed file",required= True,metavar='')
parser.add_argument('--output_prefix', '-o',type= str,help="the prefix of output file",required= True,metavar='')

group=parser.add_mutually_exclusive_group(required= True)
group.add_argument('--tools_bowtie', action='store_true',help="chose tools to search")
group.add_argument('--tools_regex', action='store_true', help="chose tools to search")
group.add_argument('--tools_blast', action='store_true', help="chose tools to search")

args = parser.parse_args()
reference = args.reference
unique_length = args.unique_length
search_range = args.search_range
bed = args.bed  
output = args.output_prefix 

bowtie = args.tools_bowtie 
regex = args.tools_regex 
blast = args.tools_blast 

def tools_bowtie(seq, reference, unique_length, search_range, output):
    with open(output+"_"+str(unique_length)+"_"+str(search_range)+'_bowtie.txt', 'a') as outputf:
        tools = '/home-02/zhenyuluo/anaconda3/envs/python3/bin/bowtie'
        subprocess.call([tools, "-c", reference, seq], stdout=outputf)
    
def tools_regex(seq, reference, unique_length, search_range, output):
    with open(output+"_"+str(unique_length)+"_"+str(search_range)+'_regex.txt', 'a') as outputf:
        tools = '/home-02/zhenyuluo/script/G4script/fastaRegexFinder_v2.py'
        subprocess.call(["python3", tools, "-r", seq, "-f", reference, "-q"], stdout=outputf)
    
def tools_blast(unique_length, search_range, output):
    # with open(output+str(unique_length)+str(search_range)+'_blast.txt', 'a') as outputf
    outputf = output+"_"+str(unique_length)+"_"+str(search_range)+'_blast.txt'
    tools ="/home-02/zhenyuluo/anaconda3/envs/python3/bin/blastn"
    db="/research/rv-02/home/zhenyuluo/reference/blast/GRCh38_no_alt/GRCh38_no_alt.fa"
    subprocess.call([blastn, "-task", "blastn", "-query", 'for_blastn_len'+str(unique_length)+".txt", "-db", db, "-outfmt", "7", "-evalue", "0.00001", "-num_threads", "100", "-out", outputf])


ref_fa = pysam.Fastafile(reference)
with open(bed, 'r') as bedf:
    for line in bedf:
        lines = line.strip().split('\t')
        chr = lines[0]
        start = int(lines[1])
        end = int(lines[2])
        # search left of the 5k region
        for i in range(search_range+1 -unique_length):
            start_start = start + i
            start_end = start + unique_length+i
            # print(chr, start_start, start_end, sep='\t')
            seq = ref_fa.fetch(chr, start_start, start_end)
            if set(seq) != {'N'}:
                ## chose tools to search
                if bowtie:
                    tools_bowtie(seq, reference, unique_length, search_range, output)
                elif regex:
                    tools_regex(seq, reference, unique_length, search_range, output)
                elif blast:
                    with open('for_blastn_len'+str(unique_length)+".txt", 'a') as for_blastn:
                        print(">"+chr+":"+str(start_start)+"-"+str(start_end), file=for_blastn)
                        print(seq, file=for_blastn)
            
        # search right of the 5k region
        for i in range(search_range+1 -unique_length):
            end_start = end-i - unique_length
            end_end = end-i
            # print(chr, end_start, end_end, sep='\t')
            seq = ref_fa.fetch(chr, end_start, end_end)
            if set(seq) != {'N'}:
                ## chose tools to search
                if bowtie:
                    tools_bowtie(seq, reference, unique_length, search_range, output)
                elif regex:
                    tools_regex(seq, reference, unique_length, search_range, output)
                elif blast:
                    with open('for_blastn_len'+str(unique_length)+".txt", 'a') as for_blastn:
                        print(">"+chr+":"+str(end_start)+"-"+str(end_end), file=for_blastn)
                        print(seq, file=for_blastn)
        with open('done.txt', 'a') as done:
            print(line.strip(), "search done", sep = '\t', file=done)
ref_fa.close()

if blast:
    tools_blast(unique_length, search_range, output)

