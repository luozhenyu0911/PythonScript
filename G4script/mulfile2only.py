from __future__ import print_function
from __future__ import division
import pandas as pd
import sys

filelist = open(sys.argv[1],"r")
of = sys.argv[2]
colname=sys.argv[3]

for i,j in enumerate(filelist):
    if i < 1 :
        df15 = pd.read_csv(j.strip(),sep='\t',error_bad_lines=False)
        name = j.strip()
        df15[colname]=name
    else:
        df16 = pd.read_csv(j.strip(),sep='\t',error_bad_lines=False)
        name = j.strip()
        df16[colname]=name
        df15=pd.concat([df15,df16])
df15.to_csv(of,sep='\t',index=None)
filelist.close()
