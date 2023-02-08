#!/bin/bash

i=0;
j=0;
tmp=0;
max=489;
#543;
DIRNAME="out_";
NAMEND=".root";
cd /lustre/hades/user/kjedrzej/iHKM/ic/output/urqmd
while [ $i -le $max ]; do
	tmp=0;
	echo $j;
	while [ $tmp -le 100 ]; do
		if [ $tmp -eq 0 ]; then
			j=$[j+1];
			mkdir $DIRNAME$j/
		fi
		#echo $i$NAMEND;
		cp $i$NAMEND $DIRNAME$j/$tmp$NAMEND
		i=$[i+1];
		tmp=$[tmp+1];
	done
done
