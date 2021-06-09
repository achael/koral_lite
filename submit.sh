#!/bin/bash
#SBATCH --job-name=ko            # create a short name for your job
#SBATCH --reservation=hackathon2
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # how many instances of your command are run, total, across all nodes
#SBATCH --cpus-per-task=4        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=16G         # memory per cpu-core (4G is default)
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)


#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun -n 1 ./ko_gpu > kogpu.out
