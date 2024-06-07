
import sys
import pandas as pd
import numpy as np

files_str = sys.argv[1]
which_col = sys.argv[2]
files = files_str.split(",")
for i, j in enumerate(files):
    name = j.strip("_pass_svtypr_stats.txt")
    if i == 0:
        df1 = pd.read_csv(j, sep='\t')
        df1["Sample"]=name
        
    else:
        df2 = pd.read_csv(j, sep='\t')
        df2["Sample"]=name
        df1 = pd.concat([df1, df2])