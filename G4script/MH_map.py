from __future__ import print_function
import sys
import pysam
import subprocess

subprocess.call(["bowtie2","-p","10","-N","1","-x",sys.argv[4],"-1",sys.argv[1],"-2",sys.argv[2],"-S",sys.argv[3]+"raw.sam"])

subprocess.call(["samtools","sort","-@","10","-O","BAM","-o",sys.argv[3]+"raw.bam",sys.argv[3]+"raw.sam"])
subprocess.call(["samtools","view","-@","10","-q","20","-O","BAM","-o",sys.argv[3]+"tmp.bam",sys.argv[3]+"raw.sam"])
subprocess.call(["samtools","sort","-O","BAM","-@","10","-n","-o",sys.argv[3]+"sort_tmp.bam",sys.argv[3]+"tmp.bam"])
subprocess.call(["samtools","fixmate","-@","10","-m",sys.argv[3]+"sort_tmp.bam",sys.argv[3]+"fixmate.bam"])
subprocess.call(["samtools","sort","-@","10","-O","BAM","-o",sys.argv[3]+"fixmatesort.bam",sys.argv[3]+"fixmate.bam"])
subprocess.call(["samtools","markdup","-@","10","-r","-O","BAM",sys.argv[3]+"fixmatesort.bam",sys.argv[3]+"rmdup.bam"])
subprocess.call(["samtools","index","-@","10",sys.argv[3]+"rmdup.bam"])

subprocess.call(["samtools","view","-@","10","-c",sys.argv[3]+"raw.sam"])
subprocess.call(["samtools","view","-@","10","-F","4","-c",sys.argv[3]+"raw.sam"])
subprocess.call(["samtools","view","-@","10","-c",sys.argv[3]+"rmdup.bam"])
subprocess.call(["bamCoverage","-p","10","-b",sys.argv[3]+"rmdup.bam","-o",sys.argv[3]+"rmdup.bw","--centerReads","--binSize","20","--smoothLength","60","--normalizeUsing","RPKM"])
subprocess.call(["rm",sys.argv[3]+"raw.sam"])
subprocess.call(["rm",sys.argv[3]+"tmp.bam"])
subprocess.call(["rm",sys.argv[3]+"sort_tmp.bam"])
subprocess.call(["rm",sys.argv[3]+"fixmate.bam"])
subprocess.call(["rm",sys.argv[3]+"fixmatesort.bam"])
