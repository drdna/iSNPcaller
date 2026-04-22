# iSNPcaller
A package that simplifies the SNP calling process by allowing one to deposit genome assemblies in a single (GENOMES) directory and then running the caller. The following operations are performed:
1. Sequence header lines are standardized.
2. Repeat masking is performed using a novel algorithm that outperforms the industry standard tools.
3. All genomes are compared with one another in a pairwise fashion (Note: the program can be run in a way that allows comparison against a single reference genome. Functionality to specify this modality at run time will be implemented soon).
4. Cryptic repeats uncovered in the pairwise alignment are masked.
5. SNPs are called only in the uniquely-aligned regions of each genome.

--
# Dependencies

- [SPAdes](https://github.com/ablab/spades)

- [Trimmomatic](https://github.com/usadellab/trimmomatic)
 
- Parallel:ForkManager

Have system admin install as a standard perl module or install using conda:
```
conda install -c bioconda perl-parallel-forkmanager

```
Perl Modules:
- RM.pm
- UniqueVariants.pm

These packages are in the lib directory and are automatically called as long as the lib directory is inside the same directory as the GenomeProcessFull.pl script when it is called.

--
# Align Genomes and Call Variants

## Preprocess genomes
1. Place all genomes in a single directory (e.g. assembly-dir)
2. uncompress them if gzipped:
```bash
gunzip assembly-dir/*gz
```
4. Activate the conda environment containing Parallel::ForkManager
5. Run the GenomeProcessFull.pl script to:
- standardize assembly names
- standardize sequence headers
- generate repeat-masked genomes
```bash
perl path/to/GenomeProcessFull.pl path/to/assembly-dir
```
## Align masked genomes to reference assembly and call SNPs
```bash
perl AlignAndCallVariants.sh path/to/Reference.fasta path/to/assembly-dir
```
This will create two directories:
- assembly-dir_BLAST which houses the genome alignment files
- assembly-dir_SNP which houses the SNP calls and "SNP_counts_out(date).txt" and SNP_counts_cumulative.txt files summarizing SNP counts

Alignment files have the format:
qseqid sseqid qstart sstart qallele sallele orientation repeated?
