import pysam
bamfile = pysam.AlignmentFile('/gss1/home/lxx20180918/maize_MH/Maize-MHS-req1_uniq.bam','rb')
for read in bamfile.fetch('1',30000,40000):
    print read
    print read.reference_start
    print read.reference_name
    print read.query_alignment_start
#if bamfile.count('10',66,67):
#    print "yes"
#else:
#    print "no"
bamfile.close()
