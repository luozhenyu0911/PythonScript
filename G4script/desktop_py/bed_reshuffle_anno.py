from __future__ import print_function
from __future__ import division
from pybedtools import BedTool
import pybedtools as pybt
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import subprocess
import argparse
r=robjects.r

parser = argparse.ArgumentParser(description="This tool used to do simulation for peak annotation using random peak region",formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument('--Input','-i',type = str,help="-i <str>    Input bed file",required = True)
parser.add_argument('--txdb','-d',type = str,help="-d <str>    Input txdb file",required = True)
parser.add_argument('--random_num','-n',type = int,help="-n <int>    Random number for test",required = True)
parser.add_argument('--genome_size','-g',type = str,help="-g <str>    Genome size file in text format",required = True)
parser.add_argument('--output','-o',type = str,help="-o <str>    Output file",required = True)
args = parser.parse_args()

ChIPseeker = importr('ChIPseeker')
r['require']('GenomicFeatures')
Db=r['loadDb'](args.txdb)
simustat={}
InBed=BedTool(args.Input)

for i in range(0,args.random_num):
    InBed.shuffle(g=args.genome_size).moveto("tmp.bed")
    x1=ChIPseeker.annotatePeak("tmp.bed", tssRegion=robjects.IntVector([-1000, 0]),level="gene",overlap="all", TxDb=Db,assignGenomicAnnotation = True,genomicAnnotationPriority = robjects.StrVector(["Promoter","5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"]))
    anno=r['table'](r['gsub']("\\(.*\\)","",r['as.data.frame'](x1).rx2('annotation')))
    feature=list(r['names'](anno))
    stat=list(anno)
    for i in range(0,7):
        if str(feature[i]) in simustat:
            simustat[str(feature[i])]+=int(stat[i])
        else:
            simustat[str(feature[i])]=int(stat[i])
    subprocess.call(["rm","tmp.bed"])

outputfile=open(args.output,"w")
for name in feature:
    print(name,simustat[name]/args.random_num,sep="\t",end="\n",file=outputfile)
