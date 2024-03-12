from __future__ import print_function
from __future__ import division
import sys
def fenlei(S1):
	f1=open(S1,'r')
	name=S1.strip().split('.')[0]
	R1=open(name+'.R1','w')
	R2=open(name+'.R2','w')
	R3=open(name+'.R3','w')
	R4=open(name+'.R4','w')
	R5=open(name+'.R5','w')
	R6=open(name+'.R6','w')
	R7=open(name+'.R7','w')
	R8=open(name+'.R8','w')
	R9=open(name+'.R9','w')
	R10=open(name+'.R10','w')
	for line in f1:
		l = line.strip().split('\t')
		value= l[1].strip()
		if 0.3< float(value)<=0.35:
			R1.write(line)
		elif 0.35< float(value)<=0.4:
			R2.write(line)
		elif 0.4< float(value)<=0.45:
			R3.write(line)
		elif 0.45< float(value)<=0.5:
			R4.write(line)
		elif 0.5< float(value)<=0.55:
			R5.write(line)
		elif 0.55< float(value)<=0.6:
			R6.write(line)
		elif 0.6< float(value)<=0.65:
			R7.write(line)
		elif 0.65< float(value)<=0.7:
			R8.write(line)
		elif 0.7< float(value)<=0.75:
			R9.write(line)
		elif 0.75< float(value)<=0.8:
			R10.write(line)
			
	f1.close()
	R1.close()
	R2.close()
	R3.close()
	R4.close()
	R5.close()
	
	R6.close()
	R7.close()
	R8.close()
	R9.close()
	R10.close()
	return
def hangshu(file):
	list1=[]
	f2=open(file,'r')
	for line in f2:
		list1.append(line.strip())
	num=len(list1)
	return num

fenlei(sys.argv[1])
name1=sys.argv[1].strip().split('.')[0]
of=open(name1+'.num','w')
for i in range(1,11):
	id='.R'+str(i)
	print(name1+id,hangshu(name1+id),sep="\t",end="\n",file=of)
of.close()
	
