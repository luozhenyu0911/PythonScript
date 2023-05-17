import pandas as pd
import numpy as np
import glob
import argparse

example_text = '''example:
    python this.py 
'''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

# parser.add_argument('--path','-p',type=str,help="path to stlfr results", default="/home/ycai/hm/dev/stlfr_results/",metavar='')
group=parser.add_mutually_exclusive_group(required= True)
group.add_argument('--ycai_stlfr_results', '-y', action='store_true',help="/home/ycai/hm/dev/stlfr_results/: path to stlfr results")
group.add_argument('--zyl_stlfr_results', '-z', action='store_true', help="/home/zhenyuluo/work_path/01_stLFR/:  path to stlfr results")

parser.add_argument('--batch', '-b',type= str,help="the batch name (e.g. V350097068) of stlfr",required= True,metavar='')
args = parser.parse_args()

ycai_stlfr_results=args.ycai_stlfr_results
zyl_stlfr_results=args.zyl_stlfr_results
batch=args.batch


def get_sum(py_files, batch):
    for index1, sum_file in enumerate(py_files):
        if index1 <1:
            with open(sum_file,'r') as sum1:
            #adict[name]={}
                adict={}
                name=sum_file.strip().split("/")[-3]
                for line in sum1:
                    if line.strip().startswith("Couldn"):
                        pass
                    #print(line)
                    else:
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
                name=sum_file.strip().split("/")[-3]
                for line in sum1:
                    if line.strip().startswith("Couldn"):
                        pass
                    #print(line)
                    else:
                        lines = line.strip().split(':')
                        adict[lines[0]]=lines[1].strip()
                        adict["Sample ID"]=name
                df2=pd.DataFrame.from_dict(adict, orient='index')
                df2=df2.set_axis(df2.iloc[0],axis=1,inplace=False)
                df2.drop(['Sample ID'], inplace=True)
                df1=pd.merge(df1, df2,left_index=True,right_index=True,how='outer')
    df1.to_csv(batch+'_sum_merge.xls',sep='\t')

if ycai_stlfr_results:
    py_files = glob.glob('/research/rv-02/home/ycai/dev/data_transfer/'+batch+'/*/stLFR_Analysis/summary_report*') 
    py_files=sorted(py_files,reverse=False)
    print(py_files)
    get_sum(py_files, batch)
    
elif zyl_stlfr_results:
    py_files = glob.glob('/home/zhenyuluo/work_path/01_stLFR/'+batch+'/*/stLFR_Analysis/summary_report*') 
    py_files=sorted(py_files,reverse=False)
    print(py_files)
    get_sum(py_files, batch)
    
else:
    print('please chose your personal path')

print("all done! and see "+batch+'_sum_merge.xls')