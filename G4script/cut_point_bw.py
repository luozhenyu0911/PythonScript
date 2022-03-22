from __future__ import print_function
from __future__ import division
import sys
import argparse
import subprocess
import pysam

parser = argparse.ArgumentParser(description="This tool can convert the bam into the cut-point bw", formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('--genomesize','-g',type=str,help="Genome size file",required= True)
parser.add_argument('--bam', '-b',type= str,help="Bam file for open chromatin assay",required= True)
parser.add_argument('--outprefix', '-o',type= str,help="the prefix of output file",required= True)
args = parser.parse_args()

bamIn=args.bam
genomesize=args.genomesize
outprefix=args.outprefix

opbam=pysam.AlignmentFile(bamIn,"rb")
region={}
for read in opbam.fetch():
    chrom=read.reference_name
    if chrom in region:
        pass
    else:
        region[chrom]={}
    cutpoint=read.reference_end+1 if read.is_reverse else read.reference_start+1
    if str(cutpoint) in region[chrom]:
        region[chrom][str(cutpoint)]+=1
    else:
        region[chrom][str(cutpoint)]=1

with open(outprefix+".wig",'w') as wig:
    for Chr in region:
        wig.write("variableStep chrom="+Chr+"\n")
        for cutpoint in region[Chr]:
            wig.write(cutpoint+"\t"+str(region[Chr][cutpoint])+"\n")
    #convert wig into bigwig
    bw = outprefix+".bw"
subprocess.call(["wigToBigWig",outprefix+".wig",genomesize,bw])
subprocess.call(["rm",outprefix+".wig"])

