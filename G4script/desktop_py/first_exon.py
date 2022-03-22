from __future__ import print_function
import sys 

gff3=open(sys.argv[1],"r")
first_exon=open(sys.argv[2],"w")

first={}
for line in gff3:
    row=line.rstrip().split()
    row2=row[8].split(";")
    row3=row2[0].split(":")
    row4=row3[-1].split("_")
    gene=row4[0]
    start=int(row[3])
    end=int(row[4])
    strand=str(row[6])
    chromsome=str(row[0])
    name = str(row2[1].split("=")[1])

    if "+" in strand:
        if gene in first:
            if start < first[gene][2]:
                first[gene]=[name,chromsome,start,end,"+"]
            else:
                pass
        else:
            first[gene]=[name,chromsome,start,end,"+"]
    if "-" in strand:
        if gene in first:
            if end > first[gene][2]:
                first[gene]=[name,chromsome,start,end,"-"]
            else:
                pass
        else:
            first[gene]=[name,chromsome,start,end,"-"]

for geneid in first.keys():
    print(geneid,*first[geneid],sep="\t",end="\n",file=first_exon)
