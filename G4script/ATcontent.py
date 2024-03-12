"python this.py input.bed sample.fa > output.txt"
from __future__ import print_function
from __future__ import division
from pyfasta import Fasta
import sys
import re
fa = Fasta(sys.argv[2])
bed=open(sys.argv[1])
print("closet MHS","AT content",sep='\t')
for line in bed:
    l = line.rstrip().split()
    name=l[3]
    a = str("{'chr':"+"'"+l[0]+"'"+","+"'start':"+l[1]+","+"'stop':"+l[2]+"}")
    convertI = eval(a)
    seq = fa.sequence(eval(a), one_based=False)
    Anum = len(re.findall("A",seq,re.IGNORECASE))
    Tnum = len(re.findall("T",seq,re.IGNORECASE))
    total = len(seq)
    if total == 0:
        value=0
        print(name,value,sep="\t")
    else:
        value = float((Anum+Tnum)/total)
        print(name,value,sep="\t")
bed.close()
