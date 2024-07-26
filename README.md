# iSNPcaller
A package that simplifies the SNP calling process by allowing one to deposit genome assemblies in a single (GENOMES) directory and then running the caller. The following operations are performed:
1. Sequence header lines are standardized.
2. Repeat masking is performed using a novel algorithm that outperforms the industry standard tools.
3. All genomes are compared with one another in a pairwise fashion (Note: the program can be run in a way that allows comparison against a single reference genome. Functionality to specify this modality at run time will be implemented soon).
4. Cryptic repeats are masked.
5. SNPs are called only in the uniquely-aligned regions of each genome.

As new genomes are generated, they are added to the analysis by simply re-populating the GENOMES directory and typing "run". 

# Reference-based SNP Calling
The iSNPcaller algorithm can also be used for reference-based variant calls as follows:

1. Ensure that genome assembly filename has the strainID at the beginning and that the ID contains no spaces, underscores or periods (hyphens are recommended as separators). All filename suffixes should be .fasta.
2. Place all genomes in a single directory
3. Change into this directory:
```bash
cd <directory_name>
```
4. Convert fasta headers in all genomes to a standard format:
```bash
for f in `ls`; do perl SimpleFastHeaders_SB.pl $f ${f/\.*/}; done
```
5. Make a new directory to receive BLAST results in the same directory that houses the directory used for the genome fasta files:
```bash
mkdir ../BLASTdata
```
6. BLAST a repeat-masked version of the reference genome against the reformatted genomes:
```bash
for f in `ls *nh.fasta`; do blastn -query Masked_ref.fasta -subject $f -evalue 1e-20 -max_target_seqs 20000 -outfmt '6 qseqid sseqid qstart qend sstart send btop' > ../BLASTdata/Masked-ref.${f/_*/}.BLAST; done
```
7. Change up a directory level:
```bash
cd ..
```
8. Run the SNP caller by pointing to the directory containing the BLAST output files and naming a new directory to receive the SNP calls:
```bash
perl StrictlyUniqueSNPs.pl BLASTdata SNPdata
```
This will create files with an "_out" suffix that house the SNP calls and "SNP_counts_out(date).txt" and SNP_counts_cumulative.txt files that summarize the SNP counts.

9. Create a strain metadata table (tab- or space-delimited) with the format: strainID lineageID 1 (the last column is for use with Chromopainter software).
10. Generate a haplotypes file using the Generate_haplotypes.pl script:
Usage: perl Generate_haplotypes.pl <StrainIDList.txt> <SNP directory> <BLAST directory> <Haplotypes outfile> <path/to/masked reference genome> <sequence type - chromosome/contig/sequence>
```bash
perl Generate_haplotypes.pl StrainIDList.txt SNPdata BLASTdata HaplotypesOutfile.txt Masked_ref.fasta contig
```
This will generate two haplotypes outfiles: one that contains haplotype calls for all variant sites (missing data repesented as 9s) and a "complete" file with no "missing" data.

11. Generate FASTA file from the complete Haplotypes data:
```bash
perl Haplotypes2FASTA.pl HaplotypesOutFile.complete.txt
```
 
