# iSNPcaller
A package that simplifies the SNP calling process by allowing one to deposit genome assemblies in a single (GENOMES) directory and then running the caller. The following operations are performed:
1. Sequence headers are standardized
2. Repeat masking is performed using a novel algorithm that outperforms the industry standard tools.
3. All genomes are compared with one another in a pairwise fashion (Note: the program can be run in a way that allows comparison against a single reference genome. Functionality to specify this modality at run time will be implemented soon).
4. Cryptic repeats are masked.
5. SNPs are called only in the uniquely-aligned regions of each genome.
As new genomes are generated, they are added to the analysis by simply re-populating the GENOMES directory. 
