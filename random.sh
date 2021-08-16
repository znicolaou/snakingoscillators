for j in `seq 1 10 100`; do for i in `seq $j $((j+9))`; do echo $i; ./janus.py --num 16 --time 100000 --rtime 90000 --dt 1.0 --seed $i --sigma 0.33 --filebase data/random2/$i & done; wait; done &
