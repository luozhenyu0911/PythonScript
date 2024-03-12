import sys

f1=open(sys.argv[1],"r")
f2=open(sys.argv[2],"r")
output=open(sys.argv[3],"w")
list1=[]
for line1 in f1:
    list1.append(line1.strip())
for line2 in f2:
    if line2.strip() in list1:
        output.write(line2)
output.close()
