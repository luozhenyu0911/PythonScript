# -*- coding: UTF-8 -*-
#计算3个文件交集
from __future__ import print_function
import argparse

parser = argparse.ArgumentParser()
parser.description='Extracting insert size of a PAIRED END bam file.'
parser.add_argument("-f","--file",nargs='+',help="3 files you want to get intersection.",type=str)
parser.add_argument("-l","--label",nargs='+',help="Labels for each file.",type=str)
#parser.add_argument("-o","--outFile",help="output file name.")
args = parser.parse_args()

def writefile(list,label):
    f=open(label,'w')
    for a in list:
        f.write(a+'\n')
    f.close()
        
if args.file:
    filelist=args.file
    if len(filelist)==len(args.label):
        label=args.label
        Combinations=[label[0],label[1],label[2],\
                      "_".join([label[0],label[1]]),\
                      "_".join([label[0],label[2]]),\
                      "_".join([label[1],label[2]]),\
                      "_".join([label[0],label[1],label[2]])]
        labelCombinations=[i+'.ID' for i in Combinations]
        all=[]
        for file in filelist:
            temp=open(file,'r')
            individual=[]
            for line in temp:
                line=line.strip()
                individual.append(line)
            all.append(set(individual))
            temp.close()
        A=all[0]-(all[1]|all[2]) # “|”表示并集
        B=all[1]-(all[0]|all[2])
        C=all[2]-(all[0]|all[1])
        AB=(all[0]&all[1])-all[2] # “&”表示交集
        AC=(all[0]&all[2])-all[1]
        BC=(all[1]&all[2])-all[0]
        ABC=all[0]&all[1]&all[2]
        out=[A,B,C,AB,AC,BC,ABC]
        for i in range(len(out)):
            writefile(out[i],labelCombinations[i])
        
        