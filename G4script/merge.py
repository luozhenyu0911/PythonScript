from __future__ import print_function
import pandas as pd
import sys
down=sys.argv[1]
up=sys.argv[2]
output=sys.argv[3]
down1 = pd.read_csv(sys.argv[1], sep='\t',header=None,names=['Gene'])
up1 = pd.read_csv(sys.argv[2], sep='\t')#,header=None,names=['Gene ID','Gene id']

df1=pd.merge(down1, up1, left_on="Gene", right_on="Gene", how="left")
df1.to_csv(output,sep='\t',index=None)
df1.drop_duplicates(inplace=True)
