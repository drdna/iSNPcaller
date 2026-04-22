# iSNPcaller
A package that simplifies the SNP calling process by allowing one to deposit genome assemblies in a single (GENOMES) directory and then running the caller. The following operations are performed:
1. Sequence header lines are standardized.
2. Repeat masking is performed using a novel algorithm that outperforms the industry standard tools.
3. All genomes are compared with one another in a pairwise fashion (Note: the program can be run in a way that allows comparison against a single reference genome. Functionality to specify this modality at run time will be implemented soon).
4. Cryptic repeats uncovered in the pairwise alignment are masked.
5. SNPs are called only in the uniquely-aligned regions of each genome.

# Dependencies

Perl Modules:
- RM.pm
- UniqueVariants.pm
- Parallel:ForkManager - install using conda:
```
conda install -c bioconda perl-parallel-forkmanager

```

## Preprocess genomes
1. Place all genomes in a single directory (e.g. assembly-dir)
2. Activate the conda environment containing Parallel::ForkManager
3. Run the GenomeProcessFull.sh script to:
- standardize assembly names
- standardize sequence headers
- generate repeat-masked genomes

- 
```bash
perl GenomeProcessFull.pl assembly-dir
```
## Align masked genomes to reference assembly and call SNPs
```bash
perl AlignAndCallVariants.sh Reference.fasta assembly-dir
```
This will create two directories:
- assembly-dir_BLAST which houses the genome alignment files
- assembly-dir_SNP which houses the SNP calls and "SNP_counts_out(date).txt" and SNP_counts_cumulative.txt files summarizing SNP counts
