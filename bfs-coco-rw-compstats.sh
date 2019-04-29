#!/bin/bash

ROOT=$1

for g in "ca-AstroPh" "facebook_combined" "roadNet-TX"; do
    find $ROOT -name bfs-coco.*.$g.* | xargs head -2q | sed '3~2d' | sort -n -t , -k 1,2 > $ROOT/bfs-coco.compstats.$g.csv
done
