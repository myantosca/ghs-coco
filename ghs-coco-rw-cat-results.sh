#!/bin/bash

ROOT=$1

for g in "ca-AstroPh" "facebook_combined" "roadNet-TX"; do
    for ((k = 1; k <= 32; k*=2)); do
	find $ROOT -name ghs-coco.$k.$g.* | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/ghs-coco.k$k.$g.csv
    done
done
