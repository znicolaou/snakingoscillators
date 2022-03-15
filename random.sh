#!/bin/bash
if [ $# -eq 0 ]; then
	echo "need filebase"
	exit
fi
filebase0=$1
args=""
if [ $# -gt 1 ]; then
	args=$2
fi
mkdir -p $filebase0
for i in `seq 1 10000`; do 
	js=`jobs | wc -l`
	while [ $js -gt 15 ]; do
		js=`jobs | wc -l`
		sleep 1
	done
	echo $i $js
	./janus.py --num 10 --time 10000 --rtime 9000 --dt 0.1 --seed $i --sigma 0.18 --sym 1 --filebase $filebase0/$i --output 0 $args  & 
done 

wait 
tail -qn 1 $filebase0/*out.dat | sort -n > $filebase0/random.dat
