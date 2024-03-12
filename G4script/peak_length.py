from __future__ import print_function
from __future__ import division
import pandas as pd
import sys

input1=open(sys.argv[1],'r')
output=open(sys.argv[2],'w')
for i,file1 in enumerate(input1):
    if i < 1 :
        name=file1.strip()
        file1=file1.strip()
        df1=pd.read_csv(file1,sep='\t',header=None).iloc[:,1:3]
        df1["sample"]=str(name)
        df1["length"]=df1[2]-df1[1]
        del df1[2]
        del df1[1]
    else:
        name=file1.strip()
        file1=file1.strip()
        df2=pd.read_csv(file1,sep='\t',header=None).iloc[:,1:3]
        df2["sample"]=str(name)
        df2["length"]=df2[2]-df2[1]
        del df2[2]
        del df2[1]
        df1=pd.concat([df1,df2],axis=0)
df1.to_csv(output,sep='\t',index=None)
