# Bioinformatic analyses were conducted on the University of Southampton’s High Performance Computer Cluster “Iridis”. 
### FastQ files have been deposited on NCBI at BioProject PRJNA1156634.

## Tools used

[bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) to make reference genome index

[HISAT2](https://daehwankimlab.github.io/hisat2/) for aligning to reference genome of _S. pistillata_

[SAMtools](http://www.htslib.org/) Converting .sam files in preparation for alignment

[StringTie](https://ccb.jhu.edu/software/stringtie/) assemble alignments into transcripts and get counts matrix for expression analysis

## Make new directory in scratch folder
mkdir HISAT-Sp

## Make bowtie indexes for genome work in HISAT-Sp folder

bowtie2-build /.../HISAT-Sp/GCF_002571385.2_Stylophora_pistillata_v1.1_genomic.fna /.../HISAT-Sp/Spistillata_genome_ref.fasta

## Make a conda environment for rna bioinformatics pipeline
module load conda/py3-latest

conda create -n rnabioinfopipeline

conda activate rnabioinfopipeline

# Install necessary programs into rnabioinfopipeline
conda install -c bioconda hisat2

conda install -c bioconda gffcompare

## Get Stringtie
wget -c http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.5.Linux_x86_64.tar.gz

tar xzf stringtie-2.1.5.Linux_x86_64.tar.gz

# Build the gtf map specific for hisat, use a genome
HISAT2-build function http://daehwankimlab.github.io/hisat2/manual/

source activate rnabioinfopipeline

hisat2-build /.../HISAT-Sp/GCF_002571385.2_Stylophora_pistillata_v1.1_genomic.fna /.../HISAT-Sp/genome

# Perform alignment

Run hisat2 using [hisat2-loop.sh](https://github.com/tessamp/Spistillata_nutrients/blob/main/Bioinformatics_workflow/hisat2-loop.sh)

## Need to convert.sam files to .bam files using samtools do this in HISAT directory
Create new conda environment for SAMtools

conda create -p /.../samtools 
conda activate samtools
conda install samtools==1.11

Run sam tools using [samtools_loop.sh](https://github.com/tessamp/Spistillata_nutrients/blob/main/Bioinformatics_workflow/samtools_loop.sh)


# Assemble aligned reads and quantify transcripts 

Create a new directory for StringTie. 

mkdir StringTie

Run StringTie using [stingtie-loop.sh](https://github.com/tessamp/Spistillata_nutrients/blob/main/Bioinformatics_workflow/stringtie-loop.sh)

### Create a .txt "merge-files.txt" file with the names of the GTF files you want to merge 

#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --job-name=stmerge
#SBATCH --time=10:00:00
#SBATCH --nodes=1

cd $SLURM_SUBMIT_DIR
conda deactivate
module load conda/py3-latest
source activate rnabioinfopipeline

### set path to executable
stringtie_exec="/.../StringTie/stringtie-2.1.5.Linux_x86_64/stringtie"

stringtie --merge -p 8 -G GCF_002571385.2_Stylophora_pistillata_v1.1_genomic.gff -o stringtie_merged.gtf [merge-files.txt](https://github.com/tessamp/Spistillata_nutrients/blob/main/Bioinformatics_workflow/merge-files.txt)

# Compare merged GTF files to reference genome

#!/bin/sh
...

source activate rnabioinfopipeline

gffcompare -r GCF_002571385.2_Stylophora_pistillata_v1.1_genomic.gff -o compared stringtie_merged.gtf

# Compile GTF files into gene count matrices
The StringTie program includes the script 'prepDE.py' that compiles the assembly-generated files into gene count matrices. The script requires a .txt file including sample names and file path.  Which needs to be in the StringTie directory samples_for_countmat.txt

source activate rnabioinfopipeline

prepDE.py3 -g ./counts_matrix.csv -i ./[samples_for_countmat.txt](https://github.com/tessamp/Spistillata_nutrients/blob/main/Bioinformatics_workflow/samples_for_countmat.txt)

##make sure to chmod +x prepDE.py3 to make it executable (first time only)

## Now you can export to R and continue with [differential expression analysis](https://github.com/tessamp/Spistillata_nutrients/blob/main/Differential%20expression%20analysis)
