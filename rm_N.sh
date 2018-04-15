#!/bin/bash
#This script is created to remove reads which have N in umi.
read -p 'Input the list of files:' -t 30 namelist 
for n in `cat $namelist`
do
	awk 'NR%4==1&&$1!~/N/{print NR-1}' $n>${n}_line_num
	python2 /data5/xqwang/dh/bin/extract_umi/rm_N.py ${n}_line_num $n
done
