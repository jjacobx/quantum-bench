#!/bin/bash

NQUBS=(72 74 76 78 80 82 84 86 88 90)
MDIMS=(256) # 16777216
DFRAC=4
INITS=("|0..0>" "|Wn>")

for q in ${NQUBS[@]}; do
    echo "Submitting job for q=${q}..."

    # Wait if too many jobs are running
    while [ $(squeue -u $USER -h | wc -l) -gt 50 ]; do
        sleep 1
    done

    # Dynamic max dimension
    dm=$((1<<$q/2))
    dm=$(($dm/$DFRAC))
    # sbatch --nodes=1 --export=PROG=bench,QREG_SIZE=$q,MAX_DIM=$dm,DEPTH=$(($q+1)) run-itensor-energy.slurm

    # Static max dimension
    # for m in ${MDIMS[@]}; do
    #    sbatch --nodes=1 --export=PROG=bench,QREG_SIZE=$q,MAX_DIM=$m,DEPTH=$(($q+1)) run-itensor-energy.slurm
    # done

    # QFT
    for i in ${INITS[@]}; do
        sbatch --nodes=1 --export=PROG=qft,QREG_SIZE=$q,INIT=$i,DEPTH=$(($q+1)) run-itensor-energy.slurm
    done

done