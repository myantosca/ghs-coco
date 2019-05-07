#!/bin/bash
#SBATCH -p cosc6326
#SBATCH -t 01:00:00
#SBATCH -N 4
#SBATCH --ntasks-per-node 8
#SBATCH --nodelist crill-101,crill-102,crill-200,crill-201

export TS="$(date +%Y-%m-%d.%H.%M.%S)"
export DOUT="./er-results/ghs-coco/${K}/${N}/${TS}"
mkdir -p ${DOUT}

for fname in $(find ./graphs/er/${N} -name '*.ecg'); do
    bname=$(basename $fname)
    mpirun -np ${K} ./build/ghs-coco -i $fname 1> ${DOUT}/ghs-coco.${K}.$bname.out 2> /dev/null
done;



