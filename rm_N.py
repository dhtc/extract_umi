#!/usr/bin/env python
import sys,numpy
keep_line_nums_file=sys.argv[1]
file=sys.argv[2]
with open(keep_line_nums_file,'r')as f:
	keep=numpy.array(f.readlines(),dtype=numpy.int)
	keep2=keep+1
	keep3=keep2+1
	keep4=keep3+1
with open(file,'r')as f:
	reads=f.readlines()
	reads_new=[]
	for n in range(len(keep)):
		reads_new.append(reads[keep[n]])
		reads_new.append(reads[keep2[n]])
		reads_new.append(reads[keep3[n]])
		reads_new.append(reads[keep4[n]])
with open(file[:-3]+'rm_N.fq','w')as f:
	f.writelines(reads_new)