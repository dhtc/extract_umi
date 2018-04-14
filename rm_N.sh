#!/bin/bash
#This script is created to remove reads which have N in umi.
for n in `cat namelist`
do
	awk 'NR%4==1&&$1~/N/{print NR}' $n>${n}_has_N
	for k in `cat ${n}_has_N`
	do
		let m=$k+3
		sed -n '$k,$mp' $n>${n}_trim_N
	done
done