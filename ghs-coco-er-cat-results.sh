#!/bin/bash

ROOT=$1

for g in "vg" "mc"; do
    for ((k = 1; k <= 32; k*=2)); do
	find $ROOT -name ghs-coco.$k.$g.*.1-.* | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/ghs-coco.k$k.$g.1-.csv
	find $ROOT -name ghs-coco.$k.$g.*.1.*  | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/ghs-coco.k$k.$g.1.csv
	find $ROOT -name ghs-coco.$k.$g.*.1+.* | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/ghs-coco.k$k.$g.1+.csv
	find $ROOT -name ghs-coco.$k.$g.*.L-.* | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/ghs-coco.k$k.$g.L-.csv
	find $ROOT -name ghs-coco.$k.$g.*.L.*  | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/ghs-coco.k$k.$g.L.csv
	find $ROOT -name ghs-coco.$k.$g.*.L+.* | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/ghs-coco.k$k.$g.L+.csv
    done
done
