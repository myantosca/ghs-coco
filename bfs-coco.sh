#!/bin/bash
#SBATCH -p cosc6326
#SBATCH -t 01:00:00
#SBATCH -N 4
#SBATCH --ntasks-per-node 8

export TS="$(date +%Y-%m-%d.%H.%M.%S)"
export DOUT="./results/bfs-coco/${TS}"
mkdir -p ${DOUT}

for ((k = 32; k > 0; k /= 2)); do
    for fname in $(find ./graphs -name '*.ecg'); do
	bname=$(basename $fname)
	mpirun -np $k ./build/bfs-coco -i $fname 1> ${DOUT}/bfs-coco.$k.$bname.out 2> /dev/null
    done;
done;
