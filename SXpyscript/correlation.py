from __future__ import print_function
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import argparse
plt.switch_backend('agg')

parser = argparse.ArgumentParser(description="This tool can compute Correlation from RNA-seq fpkm values, and plot heatmap by seaborn", formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('--fpkmlist', '-i',type= str,help="Multiple fpkm files(compute by stringtie software) write into one file.",required= True)
parser.add_argument('--outmatrix', '-o',type= str,help="Output correlation matrix file.",required= True)
parser.add_argument('--labellist', '-l',type=str,help="{Optional parameter} Labels for heatmap, if this parameter exist, use these as heatmap labels, otherwise use fpkmlists's prefix.")
parser.add_argument('--heatmapfig', '-f',type= str,help="Heatmap figure (Should be PNG format)",required= True)
args = parser.parse_args()

def outdata1(fpl,omat,hfig):
    fpkmlist=open(fpl,"r")
    output=open(omat,"w")
    fpkmall={}
    sampleID=[]
    sampleID.append("GeneID")
    for filedir in fpkmlist:
        fpkmtab=open(filedir.rstrip(),"r")
        prefix=filedir.split('.')
        sampleID.append(prefix[0])
        for line in fpkmtab.readlines():
            row=line.rstrip().split('\t')
            if "FPKM" in row[7]:
                pass
            else:
                if row[0] in fpkmall:
                    fpkmall[row[0]].append(row[7])
                else:
                    fpkmall[row[0]]=[]
                    fpkmall[row[0]].append(row[7])
    print(*sampleID,sep="\t",end="\n",file=output)
    for gene in fpkmall.keys():
        print(gene,end='\t',file=output)         
        print(*fpkmall[gene],sep="\t",end="\n",file=output)
    output.close()     
    fpkmlist.close()
    oput=open(omat,"r")
    ofig=open(hfig,"w")
    data = pd.read_table(oput,sep="\t",delimiter="\t")
    dfData = data.corr()
    sns.set(font_scale=2.2)
         #ax = plt.subplots(figsize=(160, 160))
          #ax.set_xticklabels(dfData,rotation='horizontal')
#    plt.figure(figsize = (300,300))
        #ax.set_xticklabels(dfData,rotation='horizontal')
        #label_y = ax.get_yticklabels()
    sns.clustermap(dfData,cmap="Reds",vmin=.7,vmax=1,linewidths=.75,annot=True,figsize=(20,20))
    plt.savefig(ofig,dpi=300,bbox_inches='tight')
    return

def outdata2(fpl,omat,labl,hfig):
    fpkmlist=open(fpl,"r")
    output=open(omat,"w")
    fpkmall={}
    sampleID=[]
    sampleID.append("GeneID")
    for filedir in fpkmlist:
        fpkmtab=open(filedir.rstrip(),"r")
       # prefix=filedir.split('.')
       # sampleID.append(prefix[0])
        for line in fpkmtab.readlines():
            row=line.rstrip().split('\t')
            if "FPKM" in row[7]:
                pass
            else:
                if row[0] in fpkmall:
                    fpkmall[row[0]].append(row[7])
                else:
                    fpkmall[row[0]]=[]
                    fpkmall[row[0]].append(row[7])
    labels = open(labl,"r")
    for line in labels:
        Line = line.rstrip()
        sampleID.append(Line)
    print(*sampleID,sep="\t",end="\n",file=output)
    for gene in fpkmall.keys():
        print(gene,end='\t',file=output)
        print(*fpkmall[gene],sep="\t",end="\n",file=output)
    output.close()
    fpkmlist.close()
    labels.close()
    oput=open(omat,"r")
    ofig=open(hfig,"w")
    data = pd.read_table(oput)
    dfData = data.corr()
    sns.set(font_scale=3)
         #ax = plt.subplots(figsize=(160, 160))
          #ax.set_xticklabels(dfData,rotation='horizontal')
  #  plt.figure(figsize = (300,300))
        #ax.set_xticklabels(dfData,rotation='horizontal')
    sns.clustermap(dfData,cmap="Reds",vmax=1,linewidths=.75,annot=True,figsize=(20,20),fmt="d")
        #label_y = ax.get_yticklabels()
    plt.savefig(ofig,dpi=300,bbox_inches='tight')
    return

if args.labellist is None:
    outdata1(args.fpkmlist,args.outmatrix,args.heatmapfig)
else:
    outdata2(args.fpkmlist,args.outmatrix,args.labellist,args.heatmapfig)

