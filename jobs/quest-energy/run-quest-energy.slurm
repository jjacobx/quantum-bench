#!/bin/bash

#SBATCH --job-name=quest-energy
#SBATCH --time=01:00:00
#SBATCH --account=ta070-kubta070
#SBATCH --qos=standard

#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=128

#SBATCH --error=out/run-%j.err
#SBATCH --output=out/run-%j.out
#SBATCH --partition=standard


: ${OMP_NUM_THREADS=128}
export OMP_NUM_THREADS=${OMP_NUM_THREADS}
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
: ${DEPTH=16}

# Print program environment
echo "Program environment: "
echo "PROG=${PROG}"
echo "QREG_SIZE=${QREG_SIZE}"
echo "DEPTH=${DEPTH}"
echo


# Set executable with arguments
DIR=../../build
if [ ${PROG} == "qft" ]; then
    EXE="${DIR}/${PROG} -q ${QREG_SIZE}"
elif [ ${PROG} == "rand" ]; then
    EXE="${DIR}/${PROG} -q ${QREG_SIZE} -d ${DEPTH}"
else
    echo "Unrecognised program ${PROG}!" 1>&2
    exit 1
fi

# Run executable in parallel
echo "Running ${EXE}: "
srun -Q --hint=nomultithread ${EXE}

# Wait for energy results
sleep 10
echo

# Display energy results
sacct -j $SLURM_JOB_ID.0 --format=JobID,Elapsed,ConsumedEnergy
