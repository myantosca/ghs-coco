#!/bin/bash

for ((k = 1; k <= 32; k++)); do
    for ((n = 10; n <= 20; k++)); do
	env N=$n K=$k sbatch ghs-coco-er-kn.sh
	while [[ $(squeue | grep stud08 | wc -l) != '0' ]]; do sleep 1; done;
    done;
done;

for ((k = 1; k <= 32; k++)); do
    env N=$n K=$k sbatch ghs-coco-rw-kn.sh
    while [[ $(squeue | grep stud08 | wc -l) != '0' ]]; do sleep 1; done;
done;
