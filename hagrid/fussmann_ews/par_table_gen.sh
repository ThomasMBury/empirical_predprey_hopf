#!/bin/bash

rm par_table.txt
touch par_table.txt

declare -a SPAN_VALS=(80);
declare -a RW_VALS=(1);
declare -a HAM_LENGTH_VALS=(40 80);
declare -a HAM_OFFSET_VALS=(0.5);
declare -a W_CUTOFF_VALS=(0.8 1);
declare -a SWEEP_VALS=('true');
declare -a BLOCK_SIZE_VALS=(20 40);
declare -a BS_TYPE_VALS=('Stationary' 'Circular');
declare -a N_SAMPLES_VALS=(100);


echo "span rw ham_length ham_offset w_cutoff sweep block_size bs_type n_samples" >> par_table.txt;


for span in "${SPAN_VALS[@]}"; do
for rw in "${RW_VALS[@]}"; do
for ham_length in "${HAM_LENGTH_VALS[@]}"; do
	for ham_offset in "${HAM_OFFSET_VALS[@]}"; do
	for w_cutoff in "${W_CUTOFF_VALS[@]}"; do
	for sweep in "${SWEEP_VALS[@]}"; do
	for block_size in "${BLOCK_SIZE_VALS[@]}"; do
		for bs_type in "${BS_TYPE_VALS[@]}"; do
		for n_samples in "${N_SAMPLES_VALS[@]}"; do

			echo "$span $rw $ham_length $ham_offset $w_cutoff $sweep $block_size $bs_type $n_samples" >> par_table.txt;

		done
		done
	done
	done
	done
	done
done
done
done

