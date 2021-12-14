#!/bin/bash
#SBATCH --job-name=GDarray
#SBATCH --workdir=/scratch/kj375/GD
#SBATCH --time=0:15:00  	# make this time as short as possible, so your jobs will start faster
#SBATCH --partition=all			
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12000
#SBATCH --array=1-550

# set directories
GDdir="/scratch/kj375/GD"
dataDir="/scratch/kj375/CDPOP/mmdata"

# run this to make a list of CDPOP output folder names
# ls /scratch/kj375/CDPOP/mmdata/out_* > CDPOPoutlist

# this sed command grabs the nth line from the file and assigns it to the nth element in the array, thus giving you the input file
CDPOPout=$(sed -n "$SLURM_ARRAY_TASK_ID"p CDPOPoutlist)
CDPOPoutDir=$(basename $CDPOPout)

# load R and run script
module load R
srun Rscript GDstats.R $dataDir $GDdir $CDPOPoutDir
