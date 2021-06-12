#!/bin/bash

maxjobs=16
ZGN_N=128
ZGN_t=5000
ZGN_t1=4000
ZGN_dt=0.01
ZGN_freq=3.4
ZGN_amp=0.047
ZGN_init=0.5
ZGN_num=128
ZGN_seed0=1

if [ $# -eq 0 ]; then
	echo "need filebase"
	exit
fi
filebase0=$1
args=""
if [ $# -gt 1 ]; then
	args=$2
fi
for s in `seq $ZGN_seed0 $((ZGN_seed0+ZGN_num))`; do
	mkdir -p $filebase0/$s
	./pendula.py --frequency $ZGN_freq --amplitude $ZGN_amp --init $ZGN_init --cycles $ZGN_t --outcycle $ZGN_t1 --dt $ZGN_dt --num $ZGN_N --seed $s --filebase $filebase0/$s/ --verbose 0 &
	js=`jobs | wc -l`
	echo $s $js
	while [ $js -ge $maxjobs ]; do
	       sleep 0.1
		js=`jobs | wc -l`
	done
done
wait

touch ${filebase0}.dat
rm ${filebase0}.dat
for s in `seq $ZGN_seed0 $((ZGN_seed0+ZGN_num))`; do
	cat $filebase0/$s/out.dat | head -n3 | tail -n1 >> ${filebase0}.dat
done
rm -r $filebase0
