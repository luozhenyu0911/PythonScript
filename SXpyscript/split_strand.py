from __future__ import print_function
import sys
import pysam

samfile=pysam.AlignmentFile(sys.argv[1],'rb')

outf1=pysam.AlignmentFile(sys.argv[2], "wb", template = samfile)
outf2=pysam.AlignmentFile(sys.argv[3], "wb", template = samfile)
    
for read in samfile.fetch():
    flag=int(read.flag)
    if flag == 163 or flag == 83:
        outf1.write(read)
    elif flag == 99 or flag == 147:
        outf2.write(read)

outf1.close()
outf2.close()

samfile.close()
