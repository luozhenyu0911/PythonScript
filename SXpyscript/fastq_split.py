from __future__ import print_function
import sys
import gzip

fastq = gzip.open(sys.argv[1],'rb')
basename=sys.argv[1].split('.')[0]
read1=basename+'_1.fastq.gz'
read2=basename+'_2.fastq.gz'

"""
gzip.open('','wb')
"""

with gzip.GzipFile(filename=read1,mode='wb',compresslevel=1) as file1,\
     gzip.GzipFile(filename=read2,mode='wb',compresslevel=1) as file2:
    for line in fastq:
        if line.startswith('@'):
            header = line.split(' ')[0].split('.')[-1]
            if header == '1':
                file1.write(line)
                file1.write(next(fastq))
                next(fastq)
                file1.write('+\n')
                file1.write(next(fastq))
            else:
                file2.write(line)
                file2.write(next(fastq))
                next(fastq)
                file2.write('+\n')
                file2.write(next(fastq))

fastq.close()
