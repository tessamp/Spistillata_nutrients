### Assemble aligned reads and quantify transcripts 

#Create a new directory for StringTie. 

#do this in the already made StringTie directory


#!/bin/sh
#user HPC specifics

cd $SLURM_SUBMIT_DIR
module load conda/py3-latest

#set path to executable
stringtie_exec="/.../StringTie/stringtie-2.1.5.Linux_x86_64/stringtie"

#set the path to the folder containing the .bam files
bam_folder="/.../StringTie/"

#set the output folder for .gtf files
output_folder="/.../StringTie/"

# Set the path to the reference genome GTF file
reference_gtf="/.../StringTie/GCF_002571385.2_Stylophora_pistillata_v1.1_genomic.gff"

# Loop through the .bam files
for i in {1..32}; do
    input_bam="${bam_folder}SP${i}output.bam"
    output_gtf="${output_folder}SP${i}output.gtf"
    gene_abundance="${output_folder}SP${i}_gene_abundance.tab"
    
    # Run StringTie
    $stringtie_exec $input_bam -A $gene_abundance -e -G $reference_gtf --rf -o $output_gtf
done


stringtie HH1output.bam -A HH1_gene_abundance.tab -e -G GCF_002571385.2_Stylophora_pistillata_v1.1_genomic.gff --rf -o HH1output.gtf 

sbatch stringtie-loop.sh
