from __future__ import print_function 
from __future__ import division
import sys
mulfile=open(sys.argv[1],"r")
sinfile=open(sys.argv[2],"w")
seq = {}
for line in mulfile:
    if line.startswith('>'):
        name = line.rstrip()
        seq[name] = []
    else:
        seq[name].append(line.strip())
mulfile.close()
for k,v in seq.items():
    print(k.strip(),file=sinfile)
    # for line1 in v:
        # line2=line1.strip()
    print(''.join(v),end='\n',file=sinfile)
sinfile.close()
