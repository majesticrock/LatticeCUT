#!/bin/bash

## For more information on the following SBATCH commands, please refer to the manual page (man sbatch) in the terminal
## or look up the slurm documentation under https://slurm.schedmed.com/sbatch.html

## Mandatory:
#SBATCH --job-name=sc_cut
#SBATCH --output=/home/althueser/phd/cpp/LatticeCUT/output_sc.txt
#SBATCH --time=8:00:00         ## maximum runtime; hours:minutes:seconds
#SBATCH --partition=med        ## choose queue

#SBATCH --ntasks=1              ## sets number of total mpi processes
#SBATCH --tasks-per-node=1      ## mpi processes spawned per node
#SBATCH --cpus-per-task=12       ## logical cores used per mpi process: should be 1, except if you want to combine OpenMP with$

#SBATCH --mem=63gb              ## sets maximum allowed memory per node
##SBATCH --mem-per-cpu=100mb    ## sets maximum allowed memory per mpi process

#SBATCH --mail-user=joshua.althueser@tu-dortmund.de     ## replace mail by personal mail address
#SBATCH --mail-type=END ## most relevant options: NONE, BEGIN, END, FAIL

## Optional:
#SBATCH --hint=nomultithread    ## deactivate Hyperthreading (recommended); for Hyperthreading comment out this line
#SBATCH --constraint=CascadeLake    ## chose a specific feature, e.g., only nodes with Haswell-architecture
                                ## Feature-Output by "cat /etc/slurm/slurm.conf | grep Feature"
module purge
loadmkl
cd /home/althueser/phd/cpp/LatticeCUT/        ## go to working directory

date
echo "--- START ---"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
## execute binary using mpirun on the allocated computation resource; the number of cores is $
./build_CascadeLake/latticecut params/cluster/single_peak.config

echo "--- END ---"
date
echo
echo "$(whoami) is leaving from $(hostname) ..."
echo
