from glob import *
import pandas as pd
import numpy as np
dpas = ['_0','_3','_8','_12','_15','_18','_ye','_YS']
gene = dict()
for dpa in dpas:
        dpa2 = dpa.split('_')[1]
        filenames = glob('%s*-gtf'%dpa)
        filenames2 = glob('%s*_'%dpa2)
        for filename in filenames:
            rep = filename.split('-')[1]
            rep = rep.split('_')[0]
            with open(filename+'/genes.fpkm_tracking') as fp:
                for i in fp:
                    inf = i.rstrip().split('\t')
                    if inf[0] != 'tracking_id':
                        inf[0] = inf[0].replace('evm.TU.','')
                        if inf[0] in gene:
                            gene[inf[0]]['%srep%s'%(dpa2,rep)] = float(inf[9])
                        else:
                            gene[inf[0]] = dict()
                            gene[inf[0]]['%srep%s'%(dpa2,rep)] = float(inf[9])
            fp.close()
        
        for filename in filenames2:
            rep = filename.split('-')[1]
            rep = rep.split('_')[0]
            with open(filename+'/genes.fpkm_tracking') as fp:
                for i in fp:
                    inf = i.rstrip().split('\t')
                    if inf[0] != 'tracking_id':
                        inf[0] = inf[0].replace('evm.TU.','')
                        if inf[0] in gene:
                            gene[inf[0]]['%srepadd%s'%(dpa2,rep)] = float(inf[9])
                        else:
                            gene[inf[0]] = dict()
                            gene[inf[0]]['%srepadd%s'%(dpa2,rep)] = float(inf[9])
            fp.close()
        
    print('['+datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")+']', end=' ')
    print('cal pearson')
    data = pd.DataFrame(gene)
    cotton_gene = list(data.columns)
    for i in cotton_gene:
            if sum(data[i]) == 0:
                data = data.drop([i],axis=1)
                    
        df = data.corr()
