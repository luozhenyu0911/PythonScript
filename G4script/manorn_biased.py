'''
python this.py AvsBall_MAvalues.xls Abias1 Bbias2 ABCommon plot.txt M_value Pvalue
'''
from __future__ import print_function
from __future__ import division
import pandas as pd
import sys

input_1=sys.argv[1]
bias1=open(sys.argv[2],'w')
bias2=open(sys.argv[3],'w')
bias3=open(sys.argv[4],'w')
plottxt=open(sys.argv[5],'w')
num=float(sys.argv[6])
Pvalue=float(sys.argv[7])
f1=pd.read_csv(input_1,sep='\t')
d1=f1[(f1['M_value']>=num) & (f1['P_value']<=Pvalue)].loc[:,['chr','start','end','M_value','P_value']]
d1['strand']="+"
d2=f1[(f1['M_value']<=-num) & (f1['P_value']<=Pvalue)].loc[:,['chr','start','end','M_value','P_value']]
d2['strand']="+"
d3=f1[(f1['M_value'].abs()<=num) | (f1['P_value']>=Pvalue)].loc[:,['chr','start','end','M_value','P_value']]
d3['strand']="+"
f2=f1[(f1['M_value'].abs()>=num) & (f1['P_value']<=Pvalue)].loc[:,['M_value','A_value','P_value']]


d1.to_csv(bias1,sep='\t',header=None,index=None)
d2.to_csv(bias2,sep='\t',header=None,index=None)
d3.to_csv(bias3,sep='\t',header=None,index=None)
f2.to_csv(plottxt,sep='\t',index=None)
