from __future__ import print_function
import pysam
import re
import sys


samfile = pysam.AlignmentFile(sys.argv[1],'rb')

with pysam.AlignmentFile(sys.argv[2], "wb", template = samfile) as outf:
    for read in samfile.fetch():
        tags = str(read.tags)
	asearch = re.search(r"\'AS\'", tags)
        xsearch = re.search(r"\'XS\'", tags)
	psearch = re.search(r"\'CP\'", tags)
	if asearch and psearch:
	    if not xsearch:
	        outf.write(read)


outf.close()
samfile.close()
