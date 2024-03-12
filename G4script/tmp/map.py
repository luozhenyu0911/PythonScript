from __future__ import print_function
import sys,pysam,subprocess
import argparse

parser = argparse.ArgumentParser(description="This tool can be used to mapping reads for RNA-seq, DNase-seq, ChIP-seq. you can get bam, unique bam, bw file...", formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('--seq','-s',type=str,help="(RNA or Chip or DNase) Fill in some one of the three",required= True)
parser.add_argument('--fastq','-f',type=str,help="Input fastq files name <Absolute path> (Support for multiple files, Separated by commas),format should be .gz",required= True)
parser.add_argument('--index','-i',type=str,help="Index filename prefix",required= True)
parser.add_argument('--prefix','-p',type=str,help="Sample name prefix <for sam, bam, bw...>",required= True)
parser.add_argument('--GFF','-F',type=str,help="reference annotation to include in the merging (GTF/GFF3)only for RNAseq mapping",required= True)
parser.add_argument('--genomesize', '-g',type= str,help="Genome size <for convert bw format file(for igv browser)>",required= True)
args = parser.parse_args()

seq=args.seq
fastqfile=args.fastq
fastq = fastqfile.split(",")
index=args.index
prefix=args.prefix
GFF=args.GFF
genomesize=str(args.genomesize)

def RNA_seq_P(index,fastq,prefix,GFF):
    subprocess.call(["hisat2","-p","20","--dta","-x",index,"-1",fastq[0],"-2",fastq[1],"-S",prefix+".sam",">"+prefix+"_map.log","2>&1"])
    subprocess.call(["samtools","sort","-@","20","-O","BAM","-o",prefix+".bam",prefix+".sam"])
    subprocess.call(["rm",prefix+".sam"])
    subprocess.call(["samtools","index","-@","20",prefix+".bam"])
    subprocess.call(["stringtie","-e","-p","20","-G",GFF,"-o",prefix+".gtf","-l",prefix,"-A",prefix+"_genefpkm.txt",prefix+".bam"])
    subprocess.call(["bamCoverage","-p","20","-b",prefix+".bam","-o",prefix+".bw","--centerReads","-bs","20","--smoothLength","60","--normalizeUsingRPKM"])
    return

def DNase_seq_P(index,fastq,prefix):
    subprocess.call(["bowtie2","-p","20","-x",index,"-1",fastq[0],"-2",fastq[1],"-S",prefix+".sam",">"+prefix+"_map.log","2>&1"])
    subprocess.call(["samtools","sort","-@","20","-O","BAM","-o",prefix+".bam",prefix+".sam"])
    subprocess.call(["rm",prefix+".sam"])
    subprocess.call(["samtools","index","-@","20",prefix+".bam"])
    subprocess.call(["python","/gss1/home/zqq20190923/python_script/unique_reads_for_paired.py",prefix+".bam",prefix+"_uniq.bam"])
    subprocess.call(["samtools","index",prefix+"_uniq.bam"])
   # subprocess.call(["python","/gss1/home/zqq20190923/Popera-master/Popera.py","-d",prefix+"_uniq.bam","-n",prefix,"-t","5","--threads=20"])
    return

def ChIP_seq_P(index,fastq,prefix):
    subprocess.call(["bowtie2","-p","20","-x",index,"-1",fastq[0],"-2",fastq[1],"-S",prefix+".sam",">"+prefix+"_map.log","2>&1"])
    subprocess.call(["samtools","sort","-@","20","-O","BAM","-o",prefix+".bam",prefix+".sam"])
    subprocess.call(["rm",prefix+".sam"])
    subprocess.call(["samtools","index","-@","20",prefix+".bam"])
    subprocess.call(["python","/gss1/home/zqq20190923/python_script/unique_reads_for_paired.py",prefix+".bam",prefix+"_uniq.bam"])
    subprocess.call(["samtools","index",prefix+"_uniq.bam"])
    return

if seq is "RNA":
    RNA_seq_P(index,fastq,prefix,GFF)
elif seq is "Chip":
    DNase_seq_P(index,fastq,prefix)
    if '-g' in Option:
        subprocess.call(["bamCoverage","-p","20","-b",prefix+"_uniq.bam","-o",prefix+"_uniq.bw","--centerReads","--binSize","20","--smoothLength","60","--normalizeTo1x",genomesize])
elif seq is "DNase":
    ChIP_seq_P(index,fastq,prefix)
    if '-g' in Option:
        subprocess.call(["bamCoverage","-p","20","-b",prefix+"_uniq.bam","-o",prefix+"_uniq.bw","--centerReads","--binSize","20","--smoothLength","60","--normalizeTo1x",genomesize])

