#!/bin/bash

for ((e = 10; e <= 20; e++)); do
    n=$[2**$e]
    for ((t=0; t < 3; t++)); do
	# (ln n) / n
	p=$(echo "l($n) / $n" | bc -l)
	echo mpirun -np 2 ./build/genmpig -n $n -p $p -b 65536 -o ./graphs/er/vg.$e.L.$t.ecg 
	time mpirun -np 2 ./build/genmpig -n $n -p $p -b 65536 -o ./graphs/er/vg.$e.L.$t.ecg 
	# (1 + ε)(ln n) / n
	p=$(echo "1.2 * l($n) / $n" | bc -l)
	echo mpirun -np 2 ./build/genmpig -n $n -p $p -b 65536 -o ./graphs/er/vg.$e.L+.$t.ecg 
	time mpirun -np 2 ./build/genmpig -n $n -p $p -b 65536 -o ./graphs/er/vg.$e.L+.$t.ecg 
	# (1 - ε)(ln n) / n
	p=$(echo "0.8 * l($n) / $n" | bc -l)
	echo mpirun -np 2 ./build/genmpig -n $n -p $p -b 65536 -o ./graphs/er/vg.$e.L-.$t.ecg 
	time mpirun -np 2 ./build/genmpig -n $n -p $p -b 65536 -o ./graphs/er/vg.$e.L-.$t.ecg 

	# 1 / n
	p=$(echo "1 / $n" | bc -l)
	echo mpirun -np 2 ./build/genmpig -n $n -p $p -b 65536 -o ./graphs/er/vg.$e.1.$t.ecg 
	time mpirun -np 2 ./build/genmpig -n $n -p $p -b 65536 -o ./graphs/er/vg.$e.1.$t.ecg 

	# (1 + ε) / n
	p=$(echo "1.2 / $n" | bc -l)
	echo mpirun -np 2 ./build/genmpig -n $n -p $p -b 65536 -o ./graphs/er/vg.$e.1+.$t.ecg 
	time mpirun -np 2 ./build/genmpig -n $n -p $p -b 65536 -o ./graphs/er/vg.$e.1+.$t.ecg 

	# (1 - ε) / n
	p=$(echo "0.8 / $n" | bc -l)
	echo mpirun -np 2 ./build/genmpig -n $n -p $p -b 65536 -o ./graphs/er/vg.$e.1-.$t.ecg 
	time mpirun -np 2 ./build/genmpig -n $n -p $p -b 65536 -o ./graphs/er/vg.$e.1-.$t.ecg 
    done;    
done
