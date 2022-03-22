from __future__ import print_function 
from __future__ import division
import sys
import scipy
import scipy.stats
from scipy import stats
input = open(sys.argv[1],"r")
out = open(sys.argv[2],"w")

def sum_list(items):  
    sum_numbers = 0  
    for x in items:  
        sum_numbers += float(x)  
    return sum_numbers
    
for line in input:
    a=[]
    b=[]
    c=[]
    d=[]
    e=[]
    f=[]
    a1=[]
    a2=[]
    a3=[]
    b1=[]
    b2=[]
    b3=[]
    c1=[]
    c2=[]
    c3=[]
    d1=[]
    d2=[]
    d3=[]
    e1=[]
    e2=[]
    e3=[]
    f1=[]
    f2=[]
    f3=[]
    all_list=[]
    a1.append(float(line.strip().split('\t')[1]))
    a2.append(float(line.strip().split('\t')[2]))
    a3.append(float(line.strip().split('\t')[3]))
    b1.append(float(line.strip().split('\t')[4]))
    b2.append(float(line.strip().split('\t')[5]))
    b3.append(float(line.strip().split('\t')[6]))
    c1.append(float(line.strip().split('\t')[7]))
    c2.append(float(line.strip().split('\t')[8]))
    c3.append(float(line.strip().split('\t')[9]))
    d1.append(float(line.strip().split('\t')[10]))
    d2.append(float(line.strip().split('\t')[11]))
    d3.append(float(line.strip().split('\t')[12]))
    e1.append(float(line.strip().split('\t')[13]))
    e2.append(float(line.strip().split('\t')[14]))
    e3.append(float(line.strip().split('\t')[15]))
    f1.append(float(line.strip().split('\t')[16]))
    f2.append(float(line.strip().split('\t')[17]))
    f3.append(float(line.strip().split('\t')[18]))
    a=a1+a2+a3
    b=b1+b2+b3
    c=c1+c2+c3
    d=d1+d2+d3
    e=e1+e2+e3
    f=f1+f2+f3
    all_list=a1+a2+a3+b1+b2+b3+c1+c2+c3+d1+d2+d3+e1+e2+e3+f1+f2+f3
    if sum_list(all_list) > 0:
        res=stats.kruskal(a,b,c,d,e,f)
        if res.pvalue < 0.05:
            print(line.strip().split('\t')[0],res.pvalue,sep='\t',file=out)
        else:
            pass
    else:
        pass
input.close()
out.close()  
