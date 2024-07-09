from __future__ import print_function
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict
import os
import argparse
example_text = '''example:
    python this.py -l filelist -p tmp -n TPM,FPKM
'''
parser = argparse.ArgumentParser(description="The script is to merge the special columns.",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--file_list','-l',type=str,help="the expression file list",required= True,metavar='')
parser.add_argument('--prefix', '-p',type= str,help="the prefix of the output file",required= True,metavar='')
parser.add_argument('--colname', '-n',type= str,help="the column name of the expression file you want to merge, such as FPKM, TPM, split by comma",
                    required= True,metavar='')


def readExpr_1(tsvFileL, typeL=['TPM']):
    '''
    tsvFileL: lists of files waiting for reading
    resultD: a dictionary to save data matrix
            {'TPM':[mat1, mat2,...]
             'FPKM':[mat1, mat2, ...]}
    typeL; list of names for columns to be extracted
    '''
    
    resultD = defaultdict(list)
    for tsvFile in tsvFileL:
        expr = pd.read_csv(tsvFile.strip(), sep='\t', header=0, index_col=0)
        name = os.path.basename(tsvFile.strip().split('.')[0]) #this option is very arbitary
        for _type in typeL: 
            # add _ to type to avoid override Python inner function `type` 
            expr_type = expr.loc[:,[_type]]
            expr_type.columns = [name]
            resultD[_type].append(expr_type)
    return resultD

if __name__ == '__main__':
    plt.switch_backend('PDF')
    args = parser.parse_args()
    fpkmlist = args.file_list
    prefix = args.prefix
    colname = args.colname
    colname_list = [name for name in colname.strip().split(',')]
    with open(fpkmlist, 'r') as file_list:
        big_dict = readExpr_1(file_list, colname_list)
    for _type in colname_list:
        big_df = pd.concat(big_dict[_type], axis=1,sort=True)
        big_df.to_csv(prefix + '_' + _type + '.txt',sep='\t')
        corr = big_df.corr()
        corr.to_csv(prefix + '_' + _type + '_corr.txt',sep='\t')
        picture = prefix + '_' + _type + '_corr_heatmap.pdf'
        sns.clustermap(corr, center=0, cmap="vlag",linewidths=.75,annot=True, figsize=(20, 20))
        plt.xticks(rotation=45)
        plt.savefig(picture)
