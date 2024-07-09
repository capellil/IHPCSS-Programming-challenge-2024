#!/bin/bash

########
# NOTE #
########
# Lines starting with '#SBATCH' are special comments that SLURM uses to extract
# information about your job. If you insert a space between '#' and 'SBATCH',
# SLURM will think the line is no-longer a special comment and will ignore it.

# The name of the output file
#SBATCH -o ./output_mpi_sneaky.txt

# The partition to use, shared CPU nodes in this case
#SBATCH -p RM-512

# Jobs are capped at 30 seconds (Your code should run for ~10 seconds anyway)
#SBATCH -t 00:00:30

# The number of nodes (at most 2)
#SBATCH -N 2

# The number of MPI processes per node
#SBATCH --ntasks-per-node=1

# The number of OpenMP threads per MPI process
#SBATCH --cpus-per-task=2

# The number of OpenMP threads. If using MPI, it is the number of OpenMP threads
# per MPI process
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Place OpenMP threads on cores
export OMP_PLACES=cores;

# Keep the OpenMP threads where they are
export OMP_PROC_BIND=true;

# Load the modules needed
module load cuda/11.7.1 nvhpc/22.9 openmpi/4.0.5-nvhpc22.9

# Compile everything
make

echo $SLURM_NTASKS

# Execute the program
mpirun -n $SLURM_NTASKS --bind-to none --report-bindings ./bin/main

# After init the new_pagerank takes 0.000000000000.
# After comppute the new_pagerank takes 0.034760951996.
# This iteration takes 0.034765958786.
# 285 iterations achieved in 9.99 seconds