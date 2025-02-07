#!/bin/sh 
#input user HPC specifics, time=12:00:00 

source activate /.../.conda/envs/rnabioinfopipeline

#Set the path to your HISAT2 index
index_path="/.../hisat/tran_index"

#Set the path to your raw data folder
raw_data_folder="/.../stylophorap_cleanseq_data/00.CleanData/"

#Set the output folder
output_folder="/.../hisat"

#Set the number of threads
threads=4

# Loop through the sets of files
for set_num in {1..32}; do
    # Set the sample name
    sample="SP${set_num}"

    # Build file paths
    file1="${raw_data_folder}/${sample}_1.clean.fq.gz"
    file2="${raw_data_folder}/${sample}_2.clean.fq.gz"

    # Run HISAT2
    hisat2 -p $threads --new-summary --rf --dta -q -x $index_path -1 $file1 -2 $file2 -S "${output_folder}/${sample}_output.sam" --summary-file "${output_folder}/${sample}_output.txt"
done
