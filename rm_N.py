#!/usr/bin/env python
import sys
keep_line_nums_file=sys.argv[1]
file=sys.argv[2]
with open(keep_line_nums_file,'r')as f:
	keep_line_nums=f.readlines()
	keep2=keep_line_nums+1
	keep3=keep2+1
	keep4=keep3+1
with open(file,'r')as f:
	reads=f.readlines()
	reads_new=[]
	for n in range(len(keep_line_nums)):
		reads_new.append(keep_line_nums[n])
		reads_new.append(keep2)
		reads_new.append(keep3)
		reads_new.append(keep4)
with open(file[:-3]+'rm_N.fq','w')as f:
	f.writelines(reads_new)