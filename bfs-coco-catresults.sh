#!/bin/bash

ROOT=$1

for g in "vg" "mc"; do
    for ((n = 10; n <= 20; n++)); do
	find $ROOT -name bfs-coco.*.$g.$n.1-.* | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/bfs-coco.$g.$n.1-.csv
	find $ROOT -name bfs-coco.*.$g.$n.1.*  | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/bfs-coco.$g.$n.1.csv
	find $ROOT -name bfs-coco.*.$g.$n.1+.* | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/bfs-coco.$g.$n.1+.csv
	find $ROOT -name bfs-coco.*.$g.$n.L-.* | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/bfs-coco.$g.$n.L-.csv
	find $ROOT -name bfs-coco.*.$g.$n.L.*  | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/bfs-coco.$g.$n.L.csv
	find $ROOT -name bfs-coco.*.$g.$n.L+.* | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/bfs-coco.$g.$n.L+.csv
    done
done
