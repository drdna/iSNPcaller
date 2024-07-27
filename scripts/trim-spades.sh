#!/bin/bash
# runs trimmomatic then spades assembler

#SBATCH --time 48:00:00
#SBATCH --job-name=trim-spades
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal
#SBATCH --mem=500GB
#SBATCH --mail-type ALL
#SBATCH -A coa_farman_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user farman@uky.edu

echo "SLURM_NODELIST: "$SLURM_NODELIST

# Requires paired-end read datasets

#1. Place raw reads into RAW_READS directory

#2. Run assembly using command: $ sbatch trim-spades.sh RAW_READS SRRXXXXXXX yes

#3. Use for loop to submit assemblies in parallel: $ for f in `cat SRRlist.txt`; do sbatch trim-spades.sh RAW_READS $f yes; done


dir=$1

f=$2

mkdir $f 

cp $dir/$f*_1*f*q* $f/

cp $dir/$f*_2*f*q* $f/

cd $f

if [ $3 == 'yes' ]
then
  singularity run --app trimmomatic039 /share/singularity/images/ccs/conda/amd-conda2-centos8.sinf trimmomatic PE \
  -threads 16 -phred33 -trimlog ${f}_errorlog.txt \
  $f*_1*.f*q* $f*_2*.f*q* \
  ${f}_R1_paired.fq ${f}_R1_unpaired.fq \
  ${f}_R2_paired.fq ${f}_R2_unpaired.fq \
  ILLUMINACLIP:/project/farman_uksr/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:20:20 MINLEN:90;    # change path to fasta file of adaptor sequences as appropriate
fi

# now run spades

singularity run --app spades3155 /share/singularity/images/ccs/conda/amd-conda9-rocky8.sinf spades.py --pe1-1 ${f}_R1_paired.fq --pe1-2 ${f}_R2_paired.fq --pe1-s ${f}_R1_unpaired.fq --pe1-s ${f}_R2_unpaired.fq -o ${f}_assembly

# values for --app and path to trimmomatic/spades imags will need to be updated for your system
