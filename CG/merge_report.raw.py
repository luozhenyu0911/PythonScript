import pandas as pd
import numpy as np
import glob
import argparse

py_files = glob.glob('*/*/summary_report*') 
py_files=sorted(py_files,reverse=False)
print(py_files)

for index1, sum_file in enumerate(py_files):
    if index1 <1:
        with open(sum_file,'r') as sum1:
        #adict[name]={}
            adict={}
            name=sum_file.strip().split("/")[0]
            for line in sum1:
                #print(line)
                lines = line.strip().split(':')
                adict[lines[0]]=lines[1].strip()
                adict["Sample ID"]=name
            # transfer dict to Dataframe
            df1=pd.DataFrame.from_dict(adict, orient='index')
            df1=df1.set_axis(df1.iloc[0],axis=1,inplace=False)
            df1.drop(['Sample ID'], inplace=True)
    else:
        with open(sum_file,'r') as sum1:
            adict={}
            name=sum_file.strip().split("/")[0]
            for line in sum1:
                #print(line)
                lines = line.strip().split(':')
                adict[lines[0]]=lines[1].strip()
                adict["Sample ID"]=name
            df2=pd.DataFrame.from_dict(adict, orient='index')
            df2=df2.set_axis(df2.iloc[0],axis=1,inplace=False)
            df2.drop(['Sample ID'], inplace=True)
            df1=pd.merge(df1, df2,left_index=True,right_index=True,how='outer')

df1.to_csv('./sum_merge.xls',sep='\t')
