#!/bin/bash

# Change partition info to match configuration

#SBATCH --time 48:00:00
#SBATCH --job-name=trim-spades-mask
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=short
#SBATCH --mem=128GB
#SBATCH --mail-type ALL
#SBATCH --mail-user user@domain

echo "SLURM_NODELIST: "$SLURM_NODELIST

# Requires paired-end read datasets

# 1. Place raw reads into a directory
# e.g. cp strainID_1.fq.gz; strainID_2.fq.gz RAW_READS

# 2. Run assembly using command: $ sbatch trim-spades.sh <reads-directory> <reads-prefix>
# e.g. sbatch trim-spades.sh RAW_READS strainID

# 3. Use for loop to submit assemblies in parallel
# e.g. for prefix in `cat reads-prefix-list.txt`; do sbatch trim-spades.sh RAW_READS $prefix; done

readsdir=$1

prefix=$2

mkdir $prefix 

shopt -s extglob

cp "$readsdir/${prefix}"*_1.@(fq|fastq|fq.gz|fastq.gz) $prefix/

cp "$readsdir/${prefix}"*_2.@(fq|fastq|fq.gz|fastq.gz) $prefix/

cd $prefix

if [ $3 == 'yes' ]
then
  singularity run --app trimmomatic039 /share/singularity/images/ccs/conda/amd-conda2-centos8.sinf trimmomatic PE \
  -threads 32 -phred33 \
  "${prefix}"*_1.@(fq|fastq|fq.gz|fastq.gz) "${prefix}"*_2.@(fq|fastq|fq.gz|fastq.gz) \
  ${prefix}_R1_paired.fq.gz ${prefix}_R1_unpaired.fq.gz \
  ${prefix}_R2_paired.fq.gz ${prefix}_R2_unpaired.fq.gz \
  ILLUMINACLIP:/project/farman_uksr/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:20:20 MINLEN:90;
fi

rm ${prefix}_R1_unpaired.fq.gz
rm ${prefix}_R2_unpaired.fq.gz

# now run spades

singularity run --app spades3155 /share/singularity/images/ccs/conda/amd-conda9-rocky8.sinf \
  spades.py --pe1-1 ${prefix}_R1_paired.fq.gz --pe1-2 ${prefix}_R2_paired.fq.gz \
  -o ${prefix}_assembly


# following version uses unpaired read data also - usually produces lower quality assembly

#singularity run --app spades3155 /share/singularity/images/ccs/conda/amd-conda9-rocky8.sinf \
#  spades.py --pe1-1 ${prefix}_R1_paired.fq.gz --pe1-2 ${prefix}_R2_paired.fq.gz \
#  --pe1-s ${prefix}_R1_unpaired.fq --pe1-s ${prefix}_R2_unpaired.fq \
#  -o ${prefix}_assembly

# add genome ID to contig-level assembly filename
mv ${prefix}_assembly/contigs.fasta ${prefix}_assembly/${prefix}.fasta

# add genome ID to scaffold-level assembly filename
mv ${prefix}_assembly/scaffolds.fasta ${prefix}_assembly/${prefix}_scaffolds.fasta

perl /project/farman/SCRIPTs/GenomeProcessFull.pl ${prefix}_assembly/${prefix}_scaffolds.fasta

# singularity calls and path to trimmomatic/spades images will need to be updated for your system

