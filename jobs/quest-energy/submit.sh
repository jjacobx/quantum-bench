#!/bin/bash

NQUBS=({32..36})
NODES=(1 2 4 8 16 32 64 128)
FREQS=("low" "medium" "highm1" "high")

for q in ${NQUBS[@]}; do
    # Set minimum number of nodes to fulfill memory requirements
    if [[ $q -gt 33 ]]; then
        NODE_IDX_MIN=$(($q-32))
    else
        NODE_IDX_MIN=0
    fi

    for n in ${NODES[@]:$NODE_IDX_MIN}; do
        for f in ${FREQS[@]}; do
            echo "Submitting job for q=${q}, n=${n}, f=${f}..."

            # Wait if too many jobs are running
            while [ $(squeue -u $USER -h | wc -l) -gt 50 ]; do
                sleep 1
            done

            sbatch --nodes=$n --cpu-freq=$f --export=PROG=qft,QREG_SIZE=$q run-quest-energy.slurm
            sbatch --nodes=$n --cpu-freq=$f --export=PROG=rand,QREG_SIZE=$q,DEPTH=$(($q+1)) run-quest-energy.slurm
        done
    done
done
