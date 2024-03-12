import sys,subprocess,getopt
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('PDF')

def message():
    print("\n"+"Usage: python PeakAnno.py --ggplot2 or --matplotlib [options]"+"\n")
    print("Description:"+"\n"+"\t"+"\t"+"This tool could annotate genomic region of the peak and plot pie figure"+"\n")
    print("Options:")
    print("\t"+"-v, "+"--version"+"\t"+"show program's version number and exit")
    print("\t"+"-h, "+"--help"+"\t"+"show this help message and exit.")
    print("\t"+"--ggplot2"+"  or  "+"--matplotlib"+"\t"+"Optional tools for plot figure")
    print("\t"+"-d, "+"    <TxDb file>"+"\t"+"TxDb via makeTxDbFromGFF function from R packages (ChIPseeker).")
    print("\t"+"-b, "+"    <Bed file>"+"\t"+"Peak bed file.")
    print("\t"+"-a, "+"    <Annotation detail file>"+"\t"+"Genomic feature of the peak details file name.")
    print("\t"+"-p, "+"    <Annotation sum file>"+"\t"+"Genomic feature of the peak after statistics. (Calculate the sum)")
    print("\t"+"-l, "+"    <Figure title>"+"\t"+"Pie figure title name.")
    print("\t"+"-o, "+"    <Figure file>"+"\t"+"Pie figure name (format type is PDF)"+"\n")
    print("An example usage is:"+"\n")
    print("\t"+"python PeakAnno.py --ggplot2 -d ~/index/rice7.db -b wR-loop.bed -a wR-loop_peakanno.txt -p wR-loop_peakannosum.txt -l wR-loop -o wR-loop_peakanno.pdf"+"\n")
    print("\t"+"python PeakAnno.py --matplotlib -d ~/index/rice7.db -b wR-loop.bed -a wR-loop_peakanno.txt -p wR-loop_peakannosum.txt -l wR-loop -o wR-loop_peakanno.pdf"+"\n")
    return

opts,args = getopt.getopt(sys.argv[1:],'-h-v-d:-b:-a:-p:-o:-l:',['help','version','ggplot2','matplotlib'])
for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        message()
    if opt_name in ('-v','--version'):
        print("**.py 1.0 2018-09-06")
    if opt_name in ('-d'):
        Db = opt_value
    if opt_name in ('-b'):
        bed = opt_value
    if opt_name in ('-a'):
        annotationtxt = opt_value
    if opt_name in ('-p'):
        percentfile = opt_value
    if opt_name in ('-o'):
        figure = opt_value
    if opt_name in ('-l'):
        title = opt_value
    if opt_name not in ('-h,-v,-d,-b,-a,-p,-o,-l,--ggplot2,--matplotlib'):
        pass
Option={}
for option in opts:
    Option[option[0]]=option[1]

def plotpiefigure(percentfile,title,figure):
    anno=open(percentfile,"r")
    values=[]
    labels=[]
    for a in anno:
        b=a.rstrip().split("\t")
        values.append(int(b[1]))
        labels.append(b[0])
    anno.close()
    #explode = [0, 0.1, 0, 0]
    plt.axes(aspect=1)
    plt.pie(x=values,labels=None, autopct='%3.1f %%',shadow=False, labeldistance=1.1, startangle = 90,pctdistance = 0.6)
    plt.legend(labels,loc=0,bbox_to_anchor=(1, 0.5))
    plt.title(title,fontsize='large',fontweight='bold')
    plt.subplots_adjust(right=0.7)
    plt.savefig(figure)
    return

if '--matplotlib' in Option:
    subprocess.call(["Rscript","/gss1/home/lxx20180918/python_script/peakanno_for_matplotlib.r",Db,bed,annotationtxt,percentfile])
    plotpiefigure(percentfile,title,figure)
elif '--ggplot2' in Option:
    subprocess.call(["Rscript","/gss1/home/lxx20180918/python_script/peakanno_for_ggplot2.r",Db,bed,annotationtxt,percentfile,title,figure])
else:
    pass
