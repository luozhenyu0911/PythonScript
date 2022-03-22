from __future__ import print_function
from __future__ import division
import argparse
import pandas as pd
import numpy as np
import gzip
import sys
Ameth = sys.argv[1]
Bmeth = sys.argv[2]
outf = sys.argv[3]

name=['chr','type','site','context','tmp','ratio','all_count','count']
Amethf=pd.read_csv(Ameth, compression='gzip',header=None,sep='\t',names=name)
Amethf['name']=Amethf['chr']+ "_" + Amethf['site'].map(str)

Bmethf=pd.read_csv(Bmeth, compression='gzip',header=None,sep='\t',names=name)
Bmethf['name']=Bmethf['chr']+ "_" +Bmethf['site'].map(str)

all_data=pd.merge(Amethf, Bmethf, left_on="name", right_on="name", how="outer")
all_data.fillna(0,inplace=True)

all_data['ratio_deta']=all_data['ratio_x']-all_data['ratio_y']
all_data['all_count_deta']=all_data['all_count_x']-all_data['all_count_y']
all_data['count_deta']=all_data['count_x']-all_data['count_y']
new_data=all_data[['chr_x','type_x','site_x','context_x','tmp_x','ratio_deta','all_count_deta','count_deta']]
new_data['site_x']=new_data['site_x'].astype("int")
new_data.sort_values(by=['chr_x'],inplace=True)
new_data.to_csv(outf,header=None,index=None,sep='\t',compression="gzip")
