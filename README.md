# iSNPcaller
A package that simplifies the SNP calling process by allowing one to deposit genome assemblies in a single (GENOMES) directory and then running the caller. The following operations are performed:
1. Sequence header lines are standardized.
2. Repeat masking is performed using a novel algorithm that outperforms the industry standard tools.
3. All genomes are compared with one another in a pairwise fashion (Note: the program can be run in a way that allows comparison against a single reference genome. Functionality to specify this modality at run time will be implemented soon).
4. Cryptic repeats are masked.
5. SNPs are called only in the uniquely-aligned regions of each genome.

As new genomes are generated, they are added to the analysis by simply re-populating the GENOMES directory and typing "run". 

 # Installing and Running iSNPcaller
 ## 1. Installation
 Copy AlignAndCallVariants.pl script into a directory of your choice:
 ```bash
 cp AlignAndCallVariants.pl VariantCalling/
 ```
 Install RM.pm module in a directory in your current search paths. Find current paths with the following command:
 ```bash
 perl -e "print join($/,@INC);"
 ```
 or...
 install in a specific directory that you can add to your PERL search path:
 ```bash
 mv RM.pm PERL_MODULES
 ```
 Add the directory to your search path by including the following line in the .bashrc/.bash_profile file:
 ```bash
 export PERL5LIB=$PERL5LIB:~/PERL_MODULES
 ```
 Make the module executable:
 ```bash
 chmod +x PERL_MODULES/RM.pm
 ```
 Repeat steps for all required modules.
 
 ## 2. Initiating an iSNPcaller project
If this is the first run for a given project, initiate a project as follows:
```bash
cd iSNPcaller
perl iSNPcaller_MT.pl
```
You will be prompted for the name of a new project (should not match any file/folder in the current working directory). Type the name of the project, hit return and confirm when prompted).

...alternatively you can initiate a new project by providing the project name at the command line:
```bash
perl iSNPcaller.pl myProject
```
You will receive instructions to copy genome assemblies into the newly-created project directory. 

## 3. Adding genomes to an iSNPcaller project
Open a new terminal window and copy genome assemblies into the GENOMES folder in the newly-created project directory:
```bash
cp path/to/myGenomesLocation/*fasta myProject/GENOMES
```
## 4. Run the project
Return to the original terminal window, type run and hit return. The run dialog should now appear.

## 5. Adding new genomes to a project
After a run has finished, copy the new genomes into the myProject/GENOMES directory and run iSNPcaller again:
```bash
cp path/to/myNewGenomesLocation/*fasta iSNPcaller/myProject/GENOMES
cd iSNPcaller
perl iSNPcaller_MT.pl myProject
```

# Reference-based SNP Calling
The iSNPcaller algorithm can also be used for reference-based variant calls as follows:

## Preprocess genomes
1. Place all genomes in a single directory (e.g. assembly-dir)
2. Run the GenomeProcessFull.sh script to:
- standardize assembly names
- standardize sequence headers
- generate repeat-masked genomes
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
