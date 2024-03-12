from __future__ import print_function
from __future__ import division
from pyfasta import Fasta
import sys
import re




fa = Fasta(sys.argv[2])
#print("type","value",sep="\t")
bed=open(sys.argv[1])
name=sys.argv[1].rstrip().split(".")[0]
for line in bed:
    l = line.rstrip().split()
#    name=l[3]
    a = str("{'chr':"+"'"+l[0]+"'"+","+"'start':"+l[1]+","+"'stop':"+l[2]+"}")
    convertI = eval(a)
    seq = fa.sequence(eval(a), one_based=False)
    Cnum = len(re.findall("C",seq,re.IGNORECASE))
    Gnum = len(re.findall("G",seq,re.IGNORECASE))
    total = len(seq)
    if total == 0:
        value=0
        print(name,value,sep="\t")
    else:
        value = float((Cnum+Gnum)/total)
        print(name,value,sep="\t")
bed.close()
