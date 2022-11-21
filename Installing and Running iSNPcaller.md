 # Installing and Running iSNPcaller
 ## 1. Installation
 Copy iSNPcaller.pl script into a directory of your choice:
 ```bash
 cp iSNPcaller.pl iSNPcaller
 ```
 Install perl modules in a directory in your current search paths. Find current paths with the following commanns:
 ```bash
 perl -e "print join($/,@INC);"
 ```
 or...
 install in a specific directory that you can add to your PERL search path:
 ```bash
 mv Repeatmasterfaster.pm PERL_MODULES
 ```
 Add the directory to your search path by including the following line in the .bashrc/.bash_profile file:
 ```bash
 export PERL5LIB=~/PERL_MODULES
 ```
 Make the module executable:
 ```bash
 chmod a+x PERL_MODULES/Repeatmaskerfast.pm
 ```
 Repeat steps for all modules.
 
 ## 2. Initiating an iSNPcaller project
If this is the first run for a given project, initiate a project as follows:
```bash
cd iSNPcaller
perl iSNPcaller.pl
```
You will be prompted for the name of a new project (should not match any file/folder in the current working directory). Type the name of the project, hit return and confirm when prompted).

...alternatively you can initiate a new project by providing the project name at the command line:
```bash
perl iSNPcaller.pl myProject
```
confirm the project name when prompted. You will receive instructions to copy genome assemblies into the newly-created project directory. 

## 3. Adding genomes to an iSNPcaller project
Open a new terminal window and copy genome assemblies into the GENOMES folder in the newly-created project directory:
```bash
cp path/to/myGenomesLocation/*fasta myProject/GENOMES
```
## 4. Run the project
Return to the original terminal window, type run and hit return. The run dialog should now appear.

## 5. Adding new genomes to a project
Copy the new genomes into the myProject/GENOMES directory and run iSNPcaller:
```bash
cp path/to/myNewGenomesLocation/*fasta iSNPcaller/GENOMES
cd iSNPcaller
perl iSNPcaller.pl myProject run
```



