#!/bin/bash

## For more information on the following SBATCH commands, please refer to the manual page (man sbatch) in the terminal
## or look up the slurm documentation under https://slurm.schedmed.com/sbatch.html

## Mandatory:
#SBATCH --job-name=sc_cut
#SBATCH --output=/home/althueser/phd/cpp/LatticeCUT/output_sc.txt
#SBATCH --time=48:00:00         ## maximum runtime; hours:minutes:seconds
#SBATCH --partition=long        ## choose queue

##BATCH --ntasks=1              ## sets number of total mpi processes
##BATCH --tasks-per-node=1      ## mpi processes spawned per node
#BATCH --cpus-per-task=4       ## logical cores used per mpi process: should be 1, except if you want to combine OpenMP with$

#SBATCH --mem=60gb              ## sets maximum allowed memory per node
##SBATCH --mem-per-cpu=100mb    ## sets maximum allowed memory per mpi process

#SBATCH --mail-user=joshua.althueser@tu-dortmund.de     ## replace mail by personal mail address
#SBATCH --mail-type=END ## most relevant options: NONE, BEGIN, END, FAIL

## Optional:
#SBATCH --hint=nomultithread    ## deactivate Hyperthreading (recommended); for Hyperthreading comment out this line
#SBATCH --constraint=IceLake    ## chose a specific feature, e.g., only nodes with Haswell-architecture
                                ## Feature-Output by "cat /etc/slurm/slurm.conf | grep Feature"
module purge
cd /home/althueser/phd/cpp/LatticeCUT/        ## go to working directory

date
echo "--- START ---"

## execute binary using mpirun on the allocated computation resource; the number of cores is $
./build_IceLake/latticecut params/cluster/sc.config

echo "--- END ---"
date
echo
echo "$(whoami) is leaving from $(hostname) ..."
echo
