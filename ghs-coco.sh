#!/bin/bash

for ((k = 1; k <= 32; k *= 2)); do
    for ((n = 10; n <= 20; n++)); do
	echo "env N=$n K=$k sbatch ghs-coco-er-kn.sh"
	env N=$n K=$k sbatch ghs-coco-er-kn.sh
	while [[ $(squeue | grep stud08 | wc -l) != '0' ]]; do sleep 1; done;
    done;
done;

for ((k = 1; k <= 32; k *= 2)); do
    echo "env K=$k sbatch ghs-coco-rw-k.sh"
    env K=$k sbatch ghs-coco-rw-k.sh
    while [[ $(squeue | grep stud08 | wc -l) != '0' ]]; do sleep 1; done;
done;
