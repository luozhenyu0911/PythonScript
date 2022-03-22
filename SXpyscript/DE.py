from __future__ import print_function
import sys,pysam,subprocess,getopt

def message():
    print("\n"+"Usage: python DE.py [-h] [-v] [options] ..."+"\n")
    print("Description:"+"\n"+"\t"+"\t"+"This tool can be used to calculated diffexpression genes for two samples, Core program is DESeq2(R environment support) "+"\n")
    print("Options:"+"\n")
    print("\t"+"-v, "+"--version"+"\t"+"show program's version number and exit")
    print("\t"+"-h, "+"--help"+"\t"+"show this help message and exit.")
    print("\t"+"-g, "+"--sample1gtf"+"\t"+"<Sample1 gtf>"+"\t"+"Input files name(Support for multiple biological repetition files, Separated by commas), produced by stringtie program, format should be GTF")
    print("\t"+"-G, "+"--sample2gtf"+"\t"+"<Sample2 gtf>"+"\t"+"Input files name(Support for multiple biological repetition files, Separated by commas), produced by stringtie program, format should be GTF")
    print("\t"+"-l, "+"--sample1label"+"\t"+"<Sample1 label>"+"\t"+"Sample1 label name <for phenodata.csv>")
    print("\t"+"-L, "+"--sample2label"+"\t"+"<Sample2 label>"+"\t"+"Sample2 label name <for phenodata.csv>")
    print("\t"+"-p, "+"--prefixname"+"\t"+"<prepare file preifx name>"+"\t"+"Name prefix for prepareing files(sample.txt, phenodata.csv, gene_count_matrix.csv, transcript_count_matrix.csv)")
    print("\t"+"-o, "+"--diffexpr"+"\t"+"<output file>"+"\t"+"Sample1 vs Sample2 differential expression genes file")
    print("\n"+"An example usage is:"+"\n")
    print("\t"+"python DE.py -g sample1rep1.gtf,sample1rep2.gtf... -G sample2rep1.gtf,sample2rep2.gtf... -l * -L * -p * -o diffexpr.txt")
    return

opts,args = getopt.getopt(sys.argv[1:],'-h-v-g:-G:-l:-L:-p:-o:',['help','version','sample1gtf','sample2gtf','sample1label','sample2label','prefixname','diffexpr'])
for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        message()
    if opt_name in ('-v','--version'):
        print("**.py 1.0 2018-08-31")
    if opt_name in ('-g'):
        gtfsample1 = opt_value
	gtf1 = gtfsample1.split(",")
	g1number = len(gtf1)
        gn1 = str(g1number)
    if opt_name in ('-G'):
        gtfsample2 = opt_value
	gtf2 = gtfsample2.split(",")
	g2number = len(gtf2)
        gn2 = str(g2number)
    if opt_name in ('-l'):
        label1 = opt_value
    if opt_name in ('-L'):
	label2 = opt_value 
    if opt_name in ('-p'):
        prefix = opt_value
    if opt_name in ('-o'):
        diffexprfile = opt_value
    if opt_name not in ('-h,-v,-g,-G,-l,-L,-p,-o,--help,--version,--sample1gtf,--sample2gtf,--sample2label,--sample2label,--prefixname,--diffexpr'):
        message()
Option={}
for option in opts:
    Option[option[0]]=option[1]

def PrepareFile(gs1,gs2,g1n,g2n,g1,g2,l1,l2,df,gnr1,gnr2,pre):
    TmpDict={}
    TmpDict[l1]=gs1
    TmpDict[l2]=gs2
    sampletxt=open(pre+"_sample.txt","w")
    for i in range(0,g1n):
        print(g1[i][:-4],g1[i],sep="\t",end="\n",file=sampletxt)
    for i in range(0,g2n):
        print(g2[i][:-4],g2[i],sep="\t",end="\n",file=sampletxt)
    sampletxt.close()
    phenodata=open(pre+"_phenodata.csv","w")
    print("\"ids\"","\"lines\"",sep=",",end="\n",file=phenodata)
    for key,value in TmpDict.items():
        values=value.split(",")
        for v in values:
            print('"'+v+'"','"'+key+'"',sep=",",end="\n",file=phenodata)
    phenodata.close()
    subprocess.call(["python","/gss1/home/lxx20180918/python_script/prepDE.py","-i",pre+"_sample.txt"])
    subprocess.call(["mv","gene_count_matrix.csv",pre+"_gene_count_matrix.csv"])
    subprocess.call(["mv","transcript_count_matrix.csv",pre+"_transcript_count_matrix.csv"])
    subprocess.call(["Rscript","/gss1/home/lxx20180918/python_script/diff.r",pre+"_gene_count_matrix.csv",pre+"_phenodata.csv",df,gnr1,gnr2])
    return

if Option.has_key('-g') and Option.has_key('-G') and Option.has_key('-l') and Option.has_key('-L') and Option.has_key('-o') and Option.has_key('-p'):
    PrepareFile(gtfsample1,gtfsample2,g1number,g2number,gtf1,gtf2,label1,label2,diffexprfile,gn1,gn2,prefix)
else:
    pass

