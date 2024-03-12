from __future__ import print_function
import sys
import pandas as pd
input1=open(sys.argv[1],'r')
refer_bed=open(sys.argv[2],'w')

alist=[]
for i in range(226,510):
    w= str("scaffold_")+str(i)
    alist.append(w)

# for i in range(226,510):
for line in input1:
    line=line.strip()
    line2=line.strip().split('\t')[0]
    # w= str("scaffold_")+str(i)
    if str(line2) not in alist:
        print(line,file=refer_bed)
    else:
        pass
input1.close()
refer_bed.close()
