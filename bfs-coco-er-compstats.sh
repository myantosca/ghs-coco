#!/bin/bash

ROOT=$1

for g in "vg" "mc"; do
	find $ROOT -name bfs-coco.*.$g.*.1-.* | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/bfs-coco.compstats.$g.1-.csv
	find $ROOT -name bfs-coco.*.$g.*.1.*  | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/bfs-coco.compstats.$g.1.csv
	find $ROOT -name bfs-coco.*.$g.*.1+.* | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/bfs-coco.compstats.$g.1+.csv
	find $ROOT -name bfs-coco.*.$g.*.L-.* | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/bfs-coco.compstats.$g.L-.csv
	find $ROOT -name bfs-coco.*.$g.*.L.*  | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/bfs-coco.compstats.$g.L.csv
	find $ROOT -name bfs-coco.*.$g.*.L+.* | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/bfs-coco.compstats.$g.L+.csv
done
