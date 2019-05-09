#!/bin/bash

ROOT=$1

for g in "vg" "mc"; do
	find $ROOT -name ghs-coco.*.$g.*.1-.* | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/ghs-coco.compstats.$g.1-.csv
	find $ROOT -name ghs-coco.*.$g.*.1.*  | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/ghs-coco.compstats.$g.1.csv
	find $ROOT -name ghs-coco.*.$g.*.1+.* | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/ghs-coco.compstats.$g.1+.csv
	find $ROOT -name ghs-coco.*.$g.*.L-.* | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/ghs-coco.compstats.$g.L-.csv
	find $ROOT -name ghs-coco.*.$g.*.L.*  | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/ghs-coco.compstats.$g.L.csv
	find $ROOT -name ghs-coco.*.$g.*.L+.* | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/ghs-coco.compstats.$g.L+.csv
done
