#!/usr/bin/bash

# Usage: AlignAndCallVariants.sh <path/to/reference_masked.fasta> <genomes-dir>
# 1. Blast masked references against all files in specified directory
# 2. Call variants in uniquely aligned regions

if [[ -z "$refgenome" || -z "$genomesdir" ]]; then
    echo "Usage: $0 <path/to/reference_masked.fasta> <genomes-dir>" >&2
    exit 1
fi

refgenome=$1
genomesdir=$2
blastdir=${genomesdir}_BLAST
snpdir=${genomesdir}_SNP

refgenomebase=$(basename "$refgenome")
refgenomeID=${refgenomebase/\.*/}
refgenomeID=${refgenomeID%%_*}

mkdir -p "$blastdir"

for f in "$genomesdir"/*masked.fasta; do

    # 1. Get the filename without the path (e.g., path/to/file.fasta -> file.fasta)
    genomeID=$(basename "$f")

    # 2. Extract the prefix before the first underscore (e.g., file.fasta -> file)
    prefix=${genomeID%%_*}

    blastn -query $refgenome \
           -subject "$f" \
           -evalue 1e-20 \
           -max_target_seqs 20000 \
           -outfmt '6 qseqid sseqid qstart qend sstart send btop' \
           > "$blastdir/$refgenomeID.${prefix}.BLAST"

done

perl /project/farman_uksr/SCRIPTs/Run_SU4.pl "$blastdir" "$snpdir"





