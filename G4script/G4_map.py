from __future__ import print_function
import sys
import pysam
import subprocess

subprocess.call(["bwa","mem","-t","20","-M",sys.argv[4],sys.argv[1],sys.argv[2],"-o",sys.argv[3]+"raw.sam"])
subprocess.call(["samtools","sort","-@","20","-O","BAM","-o",sys.argv[3]+"raw.bam",sys.argv[3]+"raw.sam"])
subprocess.call(["samtools","view","-@","20","-F","4","-q","10","-O","BAM","-o",sys.argv[3]+"tmp.bam",sys.argv[3]+"raw.sam"])
subprocess.call(["samtools","sort","-O","BAM","-@","20","-o",sys.argv[3]+"sort_tmp.bam",sys.argv[3]+"tmp.bam"])
subprocess.call(["java","-jar","/gss1/biosoft/picard.jar","MarkDuplicates","I="+sys.argv[3]+"sort_tmp.bam","O="+sys.argv[3]+"rmdup.bam","M=marked_dup_metrics.txt","REMOVE_DUPLICATES=true"])
subprocess.call(["samtools","index","-@","20",sys.argv[3]+"rmdup.bam"])
bf = pysam.AlignmentFile(sys.argv[3]+"rmdup.bam","rb")

with pysam.AlignmentFile(sys.argv[3]+"_uniq.bam", "wb", template = bf) as outf:
    for read in bf.fetch():
        mes = read.query_alignment_length
        if mes > 50:
            outf.write(read)
        else:
            pass
        
bf.close()
outf.close()
subprocess.call(["samtools","index","-@","20",sys.argv[3]+"_uniq.bam"])
subprocess.call(["samtools","view","-@","20","-c",sys.argv[3]+"raw.sam"])
subprocess.call(["samtools","view","-@","20","-F","4","-c",sys.argv[3]+"raw.sam"])
subprocess.call(["samtools","view","-@","20","-c",sys.argv[3]+"_uniq.bam"])
subprocess.call(["bamCoverage","-p","20","-b",sys.argv[3]+"_uniq.bam","-o",sys.argv[3]+"_uniq.bw","--centerReads","--binSize","20","--smoothLength","60","--normalizeUsing","RPKM"])
subprocess.call(["rm",sys.argv[3]+"raw.sam"])
subprocess.call(["rm",sys.argv[3]+"tmp.bam"])
subprocess.call(["rm",sys.argv[3]+"sort_tmp.bam"])
subprocess.call(["rm",sys.argv[3]+"rmdup.bam"])
