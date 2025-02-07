###samtools


##need to convert.sam files to .bam files using samtools do this in HISAT-Sp directory

#!/bin/sh
#user HPC specifics

cd $SLURM_SUBMIT_DIR
conda deactivate
module load conda/py3-latest
source activate /.../.conda/envs/samtools

#set the path to your samtools executable
samtools_exec="/.../.conda/envs/samtools"

#set the path to the folder containing the .sam files
sam_folder="/.../HISAT-Sp/"

#set the output folder for .bam files
bam_output_folder="/.../StringTie/"

#loop through all .sam files in the folder
for sam_file in "${sam_folder}"*.sam; do
    #extract the file name without the path and extension
    base_name=$(basename -s .sam "${sam_file}")
    
    #set the output file path
    bam_output="${bam_output_folder}${base_name}.bam"
    
    #run samtools sort
    $samtools_exec sort -@ 8 -o $bam_output $sam_file

done
