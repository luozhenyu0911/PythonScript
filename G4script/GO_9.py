from __future__ import print_function
import sys

f = open(sys.argv[1],"r")
out = open(sys.argv[2],"w")

go = {}
go["0"] = {}
go["1"] = {}
go["2"] = {}
go["3"] = {}
go["4"] = {}
go["5"] = {}
go["6"] = {}
go["7"] = {}
go["8"] = {}
for line in f:
    l = line.rstrip().split("\t")
    if l[0] == "0":
        go["0"][l[1]] = l[2]
    if l[0] == "1":
        go["1"][l[1]] = l[2]
    if l[0] == "2":
        go["2"][l[1]] = l[2]
    if l[0] == "3":
        go["3"][l[1]] = l[2]
    if l[0] == "4":
        go["4"][l[1]] = l[2]
    if l[0] == "5":
        go["5"][l[1]] = l[2]
    if l[0] == "6":
        go["6"][l[1]] = l[2]
    if l[0] == "7":
        go["7"][l[1]] = l[2]
    if l[0] == "8":
        go["8"][l[1]] = l[2]
                
t1 = list(set(go["0"].keys()).union(set(go["1"].keys())))
t2 = list(set(t1).union(set(go["2"].keys())))
t3 = list(set(t2).union(set(go["3"].keys())))
t4 = list(set(t3).union(set(go["4"].keys())))
t5 = list(set(t4).union(set(go["5"].keys())))
t6 = list(set(t5).union(set(go["6"].keys())))
t7 = list(set(t6).union(set(go["7"].keys())))
t8 = list(set(t7).union(set(go["8"].keys())))

print("term","Com-ASO","Com-SAS","Com-SO","PDS-ASO","PDS-SAS","PDS-SO","PEG-ASO","PEG-SAS","PEG-SO",sep="\t",end="\n",file=out)
for i in t4:
    tmp = []
    if i in go["0"]:
        tmp.append(go["0"][i])
    else:
        tmp.append("NA")
    if i in go["1"]:
        tmp.append(go["1"][i])
    else:
        tmp.append("NA")
    if i in go["2"]:
        tmp.append(go["2"][i])
    else:
        tmp.append("NA")
    if i in go["3"]:
        tmp.append(go["3"][i])
    else:
        tmp.append("NA")
    if i in go["4"]:
        tmp.append(go["4"][i])
    else:
        tmp.append("NA")
    if i in go["5"]:
        tmp.append(go["5"][i])
    else:
        tmp.append("NA")
    if i in go["6"]:
        tmp.append(go["6"][i])
    else:
        tmp.append("NA")
    if i in go["7"]:
        tmp.append(go["7"][i])
    else:
        tmp.append("NA")
    if i in go["8"]:
        tmp.append(go["8"][i])
    else:
        tmp.append("NA")
    print(i,tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],tmp[7],tmp[8],sep="\t",end="\n",file=out)


out.close()
f.close()
