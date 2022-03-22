from __future__ import print_function 
from __future__ import division
import sys
numbp=int(sys.argv[3])
f1=open(sys.argv[1],'r')
output1=open(sys.argv[2],'w')

#up 2K
for i in f1:
    m=i.strip().split("\t")
    if str(m[5]) == "+":
        if int(m[1]) > numbp:
            m[2]=int(m[1])+numbp
            m[1] = int(m[1])-numbp
            print(*m,sep='\t',file=output1)
        else:
            m[1] = 0
            print(*m,sep='\t',file=output1)
    else:
        m[1]=int(m[2])-numbp
        m[2] = int(m[2])+numbp
        print(*m,sep='\t',file=output1)
f1.close()
output1.close()
