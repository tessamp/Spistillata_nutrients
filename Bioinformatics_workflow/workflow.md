# Bioinformatic analyses were conducted on the University of Southampton’s High Performance Computer Cluster “Iridis”. 

# Tools used

[bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) to make reference genome index
[HISAT2](https://daehwankimlab.github.io/hisat2/) for aligning to reference genome of _S. pistillata_
SAMtools


# Make bowtie indexes for genome

#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie_index
#SBATCH --nodes=1

cd $SLURM_SUBMIT_DIR
conda deactivate
module load conda/py3-latest
source activate bowtie2env

bowtie2-build /scratch/tp1y22/HISAT/GCF_002571385.2_Stylophora_pistillata_v1.1_genomic.fna /scratch/tp1y22/HISAT/Spistillata_genome_ref.fasta


# Below using HISAT make a conda environment for rna

module load conda/py3-latest

conda create -n rnabioinfopipeline
conda activate rnabioinfopipeline

# Install necessary programs into rnabioinfopipeline
conda install -c bioconda hisat2
conda install -c bioconda multiqc  
conda install -c bioconda cutadapt
conda install -c bioconda fastqc 
conda install -c bioconda gffcompare

# Get Stringtie

wget -c http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.5.Linux_x86_64.tar.gz

tar xzf stringtie-2.1.5.Linux_x86_64.tar.gz

# Build the gtf map specific for hisat, use a genome
# HISAT2-build function http://daehwankimlab.github.io/hisat2/manual/

#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --job-name=hisat-build
#SBATCH --nodes=1

cd $SLURM_SUBMIT_DIR
conda deactivate
module load conda/py3-latest
source activate rnabioinfopipeline

hisat2-build /scratch/tp1y22/HISAT/GCF_002571385.2_Stylophora_pistillata_v1.1_genomic.fna /scratch/tp1y22/HISAT/genome

#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --job-name=hisathh1
#SBATCH --time=10:00:00
#SBATCH --nodes=1

cd $SLURM_SUBMIT_DIR
conda deactivate
module load conda/py3-latest
source activate /scratch/tp1y22/.conda/envs/rnabioinfopipeline

hisat2 -p 4 --new-summary --rf --dta -q -x Spistillata_genome_ref -1 /scratch/tp1y22/stylo_rawdat_rna/X204SC23052155-Z01-F003_02/00.CleanData/HH1/HH1_EKRN230027310-1A_HN5K7DSX5_L1_1.clean.fq.gz -2 /scratch/tp1y22/stylo_rawdat_rna/X204SC23052155-Z01-F003_02/00.CleanData/HH1/HH1_EKRN230027310-1A_HN5K7DSX5_L1_2.clean.fq.gz -S output.sam --summary-file output.txt

# Need to convert.sam files to .bam files using samtools do this in HISAT directory
# Create new conda environment for SAMtools

conda create -p /.../samtools 
conda activate samtools
conda install samtools==1.11


#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --job-name=samtools
#SBATCH --time=10:00:00
#SBATCH --nodes=1

cd $SLURM_SUBMIT_DIR
module load conda/py3-latest
conda deactivate
source activate /scratch/tp1y22/.conda/envs/samtools

samtools sort -@ 8 -o /scratch/tp1y22/StringTie/HH1output.bam /scratch/tp1y22/HISAT/HH1output.sam

### Assemble aligned reads and quantify transcripts 

Create a new directory for StringTie. 

mkdir StringTie

sbatch ~/scripts/stringtie-loop.sh

#stringtie HH1output.bam -A HH1_gene_abundance.tab -e -G GCF_002571385.2_Stylophora_pistillata_v1.1_genomic.gff --rf -o HH1output.gtf 

# Create a .txt "merge-files.txt" file with the names of the GTF files you want to merge 

#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --job-name=stmerge
#SBATCH --time=10:00:00
#SBATCH --nodes=1

cd $SLURM_SUBMIT_DIR
conda deactivate
module load conda/py3-latest
source activate rnabioinfopipeline

#set path to executable
stringtie_exec="/scratch/tp1y22/StringTie/stringtie-2.1.5.Linux_x86_64/stringtie"

stringtie --merge -p 8 -G GCF_002571385.2_Stylophora_pistillata_v1.1_genomic.gff -o stringtie_merged.gtf merge-files.txt

#Compare merged GTF files to reference genome

#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --job-name=gffcomp
#SBATCH --time=10:00:00
#SBATCH --nodes=1

cd $SLURM_SUBMIT_DIR
conda deactivate
module load conda/py3-latest
source activate rnabioinfopipeline

gffcompare -r SP_transcriptome_Fake_modified.gtf -o compared stringtie_merged.gtf

# Compile GTF files into gene count matrices. The StringTie program includes the script 'prepDE.py' that compiles the assembly-generated files into gene count matrices. This script

Compiles the GTF files into gene and transcript count matrices. The StringTie program includes the script `prepDE.py` that compiles the assembly-generated files into gene and transcript count matrices. The script requires a .txt file including sample names and file path.  Which needs to be in the StringTie directory samples_for_countmat.txt

#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --job-name=gene2
#SBATCH --time=04:00:00
#SBATCH --nodes=1

cd $SLURM_SUBMIT_DIR
conda deactivate
module load conda/py3-latest
source activate rnabioinfopipeline

prepDE.py3 -g ./counts_matrix.csv -i ./samples_for_countmat.txt

##make sure to chmod +x prepDE.py3 to make it executable (first time only)

# Now you can export to R and continue with differential expression analysis
