from __future__ import print_function
import pysam
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.description='Extracting insert size of a PAIRED END bam file.'
parser.add_argument("-b","--bam",nargs='+',help="PAIRED END Bam format file. Seperated by spaces.",type=str)
parser.add_argument("-l","--label",nargs='+',help="Labels for each bam file. Seperated by spaces.",type=str)
parser.add_argument("-o","--outFile",help="output fule name.")
args = parser.parse_args()
#print(args)

if args.bam:
    bamlist=args.bam
    if args.outFile:
        outfile=open(args.outFile,'w')
        all=[]
        for bam in bamlist:
            bamfile = pysam.AlignmentFile(bam,'rb')
            sizeList=[]
            for read in bamfile.fetch():
                if 30 <= abs(read.template_length) <= 800:
                    sizeList.append(abs(read.template_length))
            all.append(sizeList)
        if len(args.bam)==len(args.label):
            print("size","number","label",sep='\t',file=outfile)
            for sublist,label in zip(all,args.label):
                for size in set(sublist):
                    print(size,sublist.count(size),label,sep='\t',file=outfile)
            outfile.close()
            
