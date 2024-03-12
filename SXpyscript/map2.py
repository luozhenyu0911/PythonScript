from __future__ import print_function
import sys,pysam,subprocess,getopt

def message():
    print("\n"+"Usage: python map.py [--RNA-seqP] or [--RNA-seqS] or [--DNase-seqS] or [--DNase-seqP] or [--ChIP-seqS] or [--ChIP-seqP] [-h] [-v] [options] ..."+"\n")
    print("Description:"+"\n"+"\t"+"\t"+"This tool can be used to mapping reads for RNA-seq, DNase-seq, ChIP-seq. you can get bam, unique bam, bw file..."+"\n")
    print("\t"+"-v, "+"--version"+"\t"+"show program's version number and exit")
    print("\t"+"-h, "+"--help"+"\t"+"show this help message and exit."+"\n")
    print("Required parameter:"+"\n")
    print("\t"+"-f, "+"--fastq"+"\t"+"<fastq files>"+"\t"+"Input fastq files name <Absolute path> (Support for multiple files, Separated by commas),format should be .gz")
    print("\t"+"-i, "+"--index"+"\t"+"<index prefix>"+"\t"+"Index filename prefix")
    print("\t"+"-p, "+"--prefix"+"\t"+"<sample name prefix>"+"\t"+"Sample name prefix <for sam, bam, bw...>"+"\n")
    print("\t"+"-F, "+"--GFF"+"\t"+"<reference GFF file>"+"\t"+"reference annotation to include in the merging (GTF/GFF3)only for RNAseq mapping")
    print("Optional parameter:"+"\n")
    print("\t"+"-g, "+"--genomesize"+"\t"+"<genome size>"+"\t"+"Genome size <for convert bw format file(for igv browser)>")
    print("\n"+"An example usage is:"+"\n")
    print("\t"+"python map.py --ChIP-seqP -f sample1pair1.fq.gz,sample1pair2.fq.gz -i * -p * -g *")
    return

opts,args = getopt.getopt(sys.argv[1:],'-h-v-f:-p:-i:-g:-F:',['help','version','fastq','prefix','index','genomesize','GFF','RNA-seqP','RNA-seqS','DNase-seqS','DNase-seqP','ChIP-seqS','ChIP-seqP'])
for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        message()
    if opt_name in ('-v','--version'):
        print("**.py 1.0 2018-09-02")
    if opt_name in ('-f'):
        fastqfile = opt_value
        fastq = fastqfile.split(",")
    if opt_name in ('-p'):
        prefix = opt_value
    if opt_name in ('-i'):
        index = opt_value
    if opt_name in ('-F'):
        GFF = opt_value
    if opt_name in ('-g'):
        GS = opt_value
        genomesize = str(GS)
    if opt_name not in ('-h,-v,-f,-p,-i,-g,-F,--help,--version,--prefix,--fastq,--index,--GFF,--genomesize,--RNA-seqP,--RNA-seqS,--DNase-seqS,--DNase-seqP,--ChIP-seqS','--ChIP-seqP'):
        pass
Option={}
for option in opts:
    Option[option[0]]=option[1]

def RNA_seq_P(index,fastq,prefix,GFF):
    subprocess.call(["hisat2","-p","40","--dta","-x",index,"-1",fastq[0],"-2",fastq[1],"-S",prefix+".sam",">"+prefix+"_map.log","2>&1"])
    subprocess.call(["samtools","sort","-@","40","-O","BAM","-o",prefix+".bam",prefix+".sam"])
    subprocess.call(["rm",prefix+".sam"])
    subprocess.call(["samtools","index","-@","40",prefix+".bam"])
    subprocess.call(["stringtie","-e","-p","40","-G",GFF,"-o",prefix+".gtf","-l",prefix,"-A",prefix+"_genefpkm.txt",prefix+".bam"])
    subprocess.call(["bamCoverage","-b",prefix+".bam","-o",prefix+".bw","--centerReads","-bs","20","--smoothLength","60","--normalizeUsingRPKM","-p","40"])
    return

def RNA_seq_S(index,fastq,prefix,GFF):
    subprocess.call(["hisat2","-p","15","--dta","-x",index,"-U",fastq[0],"-S",prefix+".sam",">"+prefix+"_map.log","2>&1"])
    subprocess.call(["samtools","sort","-@","20","-O","BAM","-o",prefix+".bam",prefix+".sam"])
    subprocess.call(["rm",prefix+".sam"])
    subprocess.call(["samtools","index",prefix+".bam"])
    subprocess.call(["stringtie","-e","-p","10","-G",GFF,"-o",prefix+".gtf","-l",prefix,"-A",prefix+"_genefpkm.txt",prefix+".bam"])
    subprocess.call(["bamCoverage","-b",prefix+".bam","-o",prefix+".bw","--centerReads","-bs","20","--smoothLength","60","--normalizeUsingRPKM"])
    return

def DNase_seq_S(index,fastq,prefix):
   # subprocess.call(["gunzip",fastq[0]])
    subprocess.call(["bowtie","-p","15","-n","1","-m","1","-S",index,fastq[0],prefix+".sam"])
    subprocess.call(["samtools","sort","-@","20","-O","BAM","-o",prefix+".bam",prefix+".sam"])
    subprocess.call(["rm",prefix+".sam"])
    subprocess.call(["samtools","index",prefix+".bam"])
    subprocess.call(["samtools","view","-@","15","-h","-F","4","-O","BAM",prefix+".bam","-o",prefix+"_uniq.bam"])
    subprocess.call(["samtools","index",prefix+"_uniq.bam"])
   # subprocess.call(["python","/gss1/home/lxx20180918/software/Popera/Popera.py","-d",prefix+"_uniq.bam","-n",prefix,"-t","5","--threads=20"])
    return

def DNase_seq_P(index,fastq,prefix):
    subprocess.call(["bowtie2","-p","20","-x",index,"-1",fastq[0],"-2",fastq[1],"-S",prefix+".sam",">"+prefix+"_map.log","2>&1"])
    subprocess.call(["samtools","sort","-@","20","-O","BAM","-o",prefix+".bam",prefix+".sam"])
    subprocess.call(["rm",prefix+".sam"])
    subprocess.call(["samtools","index",prefix+".bam"])
    subprocess.call(["python","/gss1/home/lxx20180918/python_script/unique_reads_for_paired.py",prefix+".bam",prefix+"_uniq.bam"])
    subprocess.call(["samtools","index",prefix+"_uniq.bam"])
   # subprocess.call(["python","/gss1/home/lxx20180918/software/Popera/Popera.py","-d",prefix+"_uniq.bam","-n",prefix,"-t","5","--threads=20"])
    return

def ChIP_seq_P(index,fastq,prefix):
    subprocess.call(["bowtie2","-p","20","-x",index,"-1",fastq[0],"-2",fastq[1],"-S",prefix+".sam",">"+prefix+"_map.log","2>&1"])
    subprocess.call(["samtools","sort","-@","20","-O","BAM","-o",prefix+".bam",prefix+".sam"])
    subprocess.call(["rm",prefix+".sam"])
    subprocess.call(["samtools","index","-@","20",prefix+".bam"])
    #subprocess.call(["python","/gss1/home/lxx20180918/python_script/unique_reads_for_paired.py",prefix+".bam",prefix+"_uniq.bam"])
    subprocess.call(["samtools","view","-@","20","-q","20","-o",prefix+"_uniq.bam",prefix+".bam"])
    subprocess.call(["samtools","index","-@","20",prefix+"_uniq.bam"])
    return

def ChIP_seq_S(index,fastq,prefix):
    subprocess.call(["bowtie2","-p","40","-x",index,"-U",fastq[0],"-S",prefix+".sam",">"+prefix+"_map.log","2>&1"])
    subprocess.call(["samtools","sort","-@","40","-O","BAM","-o",prefix+".bam",prefix+".sam"])
    subprocess.call(["rm",prefix+".sam"])
    subprocess.call(["samtools","index","-@","40",prefix+".bam"])
    subprocess.call(["python","/gss1/home/zpy/Script/python_script/unique_reads_for_bowtie2_single.py",prefix+".bam",prefix+"_uniq.bam"])
    subprocess.call(["samtools","index","-@","40",prefix+"_uniq.bam"])
    return

if '--RNA-seqP' in Option:
    RNA_seq_P(index,fastq,prefix,GFF)
elif '--RNA-seqS' in Option:
    RNA_seq_S(index,fastq,prefix,GFF)
elif '--DNase-seqS' in Option:
    DNase_seq_S(index,fastq,prefix)
    if '-g' in Option:
        subprocess.call(["bamCoverage","-b",prefix+"_uniq.bam","-o",prefix+"_uniq.bw","--centerReads","--binSize","20","--smoothLength","60","--normalizeTo1x",genomesize])
elif '--DNase-seqP' in Option:
    DNase_seq_P(index,fastq,prefix)
    if '-g' in Option:
        subprocess.call(["bamCoverage","-b",prefix+"_uniq.bam","-o",prefix+"_uniq.bw","--centerReads","--binSize","20","--smoothLength","60","--normalizeTo1x",genomesize])
elif '--ChIP-seqP' in Option:
    ChIP_seq_P(index,fastq,prefix)
    if '-g' in Option:
        subprocess.call(["bamCoverage","-b",prefix+"_uniq.bam","-o",prefix+"_uniq.bw","--centerReads","--binSize","20","--smoothLength","60","--normalizeTo1x",genomesize])
    else:
        pass
elif '--ChIP-seqS' in Option:
    ChIP_seq_S(index,fastq,prefix)
    if '-g' in Option:
        subprocess.call(["bamCoverage","-b",prefix+"_uniq.bam","-o",prefix+"_uniq.bw","--centerReads","--binSize","20","--smoothLength","60","--normalizeTo1x",genomesize,"-p","40"])
       # subprocess.call(["macs14","callpeak","-t",prefix+"_uniq.bam","-f","BAM","-g",genomesize,"-n",prefix,"--bw","200","--nomodel"])
    else:
        pass
else:
    pass
