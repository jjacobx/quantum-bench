#!/bin/bash

#SBATCH --job-name=itensor-energy
#SBATCH --time=01:00:00
#SBATCH --account=d35-ja
#SBATCH --qos=standard

#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=128

#SBATCH --error=out/run-%j.err
#SBATCH --output=out/run-%j.out
#SBATCH --partition=standard


#module load mkl

#: ${OMP_NUM_THREADS=128}
export OMP_NUM_THREADS=1
export OMP_PLACES=cores

# Print job environment
echo "Job environment: "
echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "SLURM_JOB_NAME=${SLURM_JOB_NAME}"
echo "SLURM_JOB_NUM_NODES=${SLURM_JOB_NUM_NODES}"
echo "SLURM_NTASKS=${SLURM_NTASKS}"
echo "SLURM_CPU_FREQ_REQ=${SLURM_CPU_FREQ_REQ}"
echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
echo


# If not set by sbatch, set the following
: ${PROG=rand}
: ${QREG_SIZE=24}
: ${INIT="|0..0>"}
: ${MAX_DIM=16777216}
: ${CUTOFF="1E-16"}
: ${DEPTH=16}

# Print program environment
echo "Program environment: "
echo "PROG=${PROG}"
echo "QREG_SIZE=${QREG_SIZE}"
echo "INIT=${INIT}"
echo "MAX_DIM=${MAX_DIM}"
echo "CUTOFF=${CUTOFF}"
echo "DEPTH=${DEPTH}"
echo


# Set executable with arguments
DIR=../../itensor-projects/bin
if [ ${PROG} == "bench" ]; then
    EXE="${DIR}/${PROG} --nq ${QREG_SIZE} --init ${INIT} --maxd ${MAX_DIM} --cut ${CUTOFF} --dep ${DEPTH}"
elif [ ${PROG} == "qft" ]; then
    EXE="${DIR}/${PROG} --nq ${QREG_SIZE} --init ${INIT} --maxd ${MAX_DIM} --cut ${CUTOFF}"
else
    echo "Unrecognised program ${PROG}!" 1>&2
    exit 1
fi

# Run executable in parallel
echo "Running ${EXE}: "
srun -Q --hint=nomultithread --distribution=block:block ${EXE}

# Wait for energy results
sleep 10
echo

# Display energy results
sacct -j $SLURM_JOB_ID.0 --format=JobID,Elapsed,ConsumedEnergy
