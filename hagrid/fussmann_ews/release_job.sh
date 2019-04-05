#!/bin/bash

mkdir -p Jobs
mkdir -p Jobs/output


. par_table_gen.sh

MAX=`cat par_table.txt | wc -l`

for i in `seq 2 $MAX`; do
	sbatch single_job_hag.sh $i
	sleep 1.0
done