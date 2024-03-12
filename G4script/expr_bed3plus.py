#coding=utf-8
from __future__ import print_function 
from __future__ import division
import sys
import pandas as pd
input1=sys.argv[1]
refer_bed=sys.argv[2]
# prefix = sys.argv[3]
sampleName=[]
GeneName=""

with open(input1, 'r') as f:  # 打开文件
    lines = f.readlines()  # 读取所有行
    first_line = lines[0]  # 取第一行
    line=first_line.strip().split('\t')
    sampleName=line[1:]
    GeneName=str(line[0])
    
for name1 in sampleName :
    f1=pd.read_csv(input1,sep='\t',usecols=[GeneName,name1])
    f1.sort_values(name1,ascending=False,inplace=True)
    # len_f1=f1.shape[0]
    # f1.index = range(len(f1))  
    #df.reset_index(drop=True)  重建索引，并删除原来的index（不然原来的index会变成新的一列） 
    # remainder = int(len_f1%3)
    # num = int(len_f1//3)
    # f1h = pd.DataFrame(f1.loc[:num-1,'Gene ID']) 
    # f1m = pd.DataFrame(f1.loc[num:2*num-1,'Gene ID']) 
    # f1l = pd.DataFrame(f1.loc[2*num:,'Gene ID']) 
    f1h = pd.DataFrame(f1[f1[name1]>10][GeneName]) 
    f1m = pd.DataFrame(f1[(f1[name1]>=1) & (f1[name1]<=10) ][GeneName]) 
    f1l = pd.DataFrame(f1[f1[name1]<1][GeneName]) 
    f2=pd.read_csv(refer_bed,sep='\t',names=['chr', 'start', 'end', GeneName,'score','strand'])

    def get_bed(genefile):
        name=pd.merge(genefile,f2,left_on=GeneName, right_on=GeneName, how="inner")
        order = ['chr', 'start', 'end', GeneName,'score','strand']
        name = name[order]
        return name
        
    f1hbed=get_bed(f1h)
    f1mbed=get_bed(f1m)
    f1lbed=get_bed(f1l)
    f1hbed.to_csv(name1+'_high.bed',sep='\t',header=None,index=None)
    f1mbed.to_csv(name1+'_mid.bed',sep='\t',header=None,index=None)
    f1lbed.to_csv(name1+'_low.bed',sep='\t',header=None,index=None)
