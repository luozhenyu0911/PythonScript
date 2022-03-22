from __future__ import print_function
import pandas as pd
import sys
filelist=open(sys.argv[1],"r")
output=open(sys.argv[2],"w")

l1=[]
for f_ in filelist:
    l1.append(f_.strip())
    
for i,j in enumerate(l1):
    if i < 1 :
        f1 = pd.read_csv(j.strip(), sep='\t',header=None,names=['TF','geneid'])
        name = j.strip().split('.')[0]
        f2=pd.DataFrame(f1["TF"].value_counts())
        f2.reset_index()
        f2[name]=f2['TF']/f2['TF'].sum()
        f3=f2.reset_index().drop('TF',axis=1)
    else:
        fa = pd.read_csv(j.strip(), sep='\t',header=None,names=['TF','geneid'])
        name = j.strip().split('.')[0]
        fb=pd.DataFrame(fa["TF"].value_counts())
        fb.reset_index()
        fb[name]=fb['TF']/fb['TF'].sum()
        fc=fb.reset_index().drop('TF',axis=1)
        df_ratings_users = pd.merge(f3, fc, left_on="index", right_on="index", how="outer")
        f3=df_ratings_users
f3 = f3.fillna(0)
f3.to_csv(output,sep='\t',index=None)
filelist.close()
output.close()
