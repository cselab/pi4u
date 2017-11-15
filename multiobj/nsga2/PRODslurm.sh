#!/bin/bash -l

#SBATCH --job-name="myJob"
#SBATCH --time=24:00:00
#SBATCH --nodes=32
#SBATCH --ntasks=32
#SBATCH --account=s658
#SBATCH --mail-user=username@ethz.ch
#SBATCH --mail-type=ALL
#SBATCH --constraint=gpu

#======START=====
export CRAY_CUDA_MPS=1
module load daint-gpu
export MPICH_MAX_THREAD_SAFETY=multiple

export TORC_WORKERS=1
export TORC_YIELDTIME=100

IFS=$'\n' read -ra arr -d '' <fish.in
echo "${arr[@]}"

srun -N $SLURM_NNODES -n $SLURM_NTASKS ./nsga2r fish.in

#=====END====
