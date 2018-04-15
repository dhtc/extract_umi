#!/bin/bash
#This script is created to remove reads which have N in umi.
for n in `cat namelist`
do
	awk 'NR%4==1&&$1~/N/{print NR}' $n>j
	for k in `cat j`
	do
		let m=$k+3
		sed -n "$k,${m}p" $n>>${n%.fq}_trim_N.fq
	done
done
