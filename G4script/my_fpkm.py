from __future__ import print_function
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import sys
plt.switch_backend('PDF')
fpkmlist=open(sys.argv[1],"r")
output=open(sys.argv[2],"w")
outcor=open(sys.argv[3],"w")
outfig=open(sys.argv[4],"w")
def readExpr_1(tsvFileL, typeL=['FPKM']):
    '''
    tsvFileL: lists of files waiting for reading
    resultD: a dictionary to save data matrix
            {'TPM':[mat1, mat2,...]
             'FPKM':[mat1, mat2, ...]}
    typeL; list of names for columns to be extracted
    '''
    resultD = {}
    for _type in typeL: 
        resultD[_type] = []
    for tsvFile in tsvFileL:
        expr = pd.read_csv(tsvFile.strip(), sep='\t', header=0, index_col=0)
        name = tsvFile.strip() #this option is very arbitary
        for _type in typeL: 
            # add _ to type to avoid override Python inner function `type` 
            expr_type = expr.loc[:,[_type]]
            expr_type.columns = [name]
            resultD[_type].append(expr_type)
    return resultD
exprD = readExpr_1(fpkmlist)
# TPM_mat = exprD['TPM']
FPKM_mat = exprD['FPKM']
test_merge_FPKM = pd.concat(FPKM_mat, axis=1,sort=True)
# dfData = test_merge_FPKM.corr()
test_merge_FPKM.to_csv(output,sep='\t')
output=open(sys.argv[2],"r")
def plotfig(op,of):
    data = pd.read_csv(op,sep='\t')
    dfData = data.corr()
    dfData.to_csv(outcor,sep='\t')
    sns.clustermap(dfData, center=0, cmap="vlag",linewidths=.75,annot=True, figsize=(20, 20))
    plt.xticks(rotation=45)
    plt.savefig(of)
plotfig(output,outfig)
# dfData.to_csv(output,sep='\t')
