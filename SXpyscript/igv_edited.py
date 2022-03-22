from __future__ import print_function
import pyBigWig
import subprocess
import argparse
parser = argparse.ArgumentParser(description="This tool used to plot igv for genome",formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument('--bwfile','-b',type = str,help="-b <str>    Input bigwig files(Support for multiple files, Must be separated by commas)",required = True)
parser.add_argument('--bedfile','-B',type = str,help="-B <str>    Input bed file",required = True)
parser.add_argument('--outfigure','-o',type = str,help="-o <str>    Out figure file prefix",required = True)
args = parser.parse_args()

bwlist=(args.bwfile).split(",")
for f in bwlist:
    bw = pyBigWig.open(f)
    tmp = open(args.outfigure+"tmp.txt","w")
    bed = open(args.bedfile,"r")
    value = []
    for line in bed:
        l = line.rstrip().split("\t")
        value = bw.values(str(l[0]),int(l[1]),int(l[2]))
    print("region","value",sep="\t",end="\n",file=tmp)
    a = 0
    for l in value:
        print(a,l,sep="\t",end="\n",file=tmp)
        a += 1
    tmp.close()
    bed.close()
    subprocess.call(["Rscript","/gss1/home/lxx20180918/python_script/igv.r",args.outfigure+"tmp.txt",args.outfigure+".pdf"])
 #   subprocess.call(["rm",args.outfigure+"tmp.txt"])
