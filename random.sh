#!/bin/bash

for i in `seq 1 1000`; do 
	js=`jobs | wc -l`
	while [ $js -gt 15 ]; do
		js=`jobs | wc -l`
		sleep 1
	done
	echo $i $js
	./janus.py --num 16 --time 100000 --rtime 90000 --seed $i --sigma 0.33 --filebase data/random/$i --output 0 & 
done 

wait 
tail -qn 1 data/random/*out.dat | sort -n > random.dat
