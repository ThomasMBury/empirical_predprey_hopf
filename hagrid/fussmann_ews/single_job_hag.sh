#!/bin/bash

# # Thomas Bury
# # PhD Candidate
# # Bauch computational epidemiology research group
# # Department of Applied Mathematics
# # Faculty of Mathematics
# # University of Waterloo

#SBATCH --mem=1000MB
#SBATCH --time=1-00:00:00
#SBATCH --output=Jobs/output/job-%j.out
#SBATCH --ntasks=1
#SBATCH --mail-user=tbury@uwaterloo.ca
#SBATCH --mail-type=BEGIN,END,FAIL

echo Job $SLURM_JOB_ID released

echo Install python modules
# python3 -m pip install -r requirements.txt


mkdir -p Jobs/job-$SLURM_JOB_ID
cd Jobs/job-$SLURM_JOB_ID

time python3 ../../script_fold_gen.py `cat ../../par_table.txt | head -n $1 | tail -n 1`
time python3 ../../script_flip_gen.py `cat ../../par_table.txt | head -n $1 | tail -n 1`

# Add parameters to a text file
cat ../../par_table.txt | head -n 1 | tail -n 1 >> pars.txt
cat ../../par_table.txt | head -n $1 | tail -n 1 >> pars.txt

cd ../../