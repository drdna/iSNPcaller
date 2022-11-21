 # Installing and Running iSNPcaller
 1. Installation
 
 a. Copy iSNPcaller.pl script into a directory of your choice:
 ```bash
 cp iSNPcaller.pl iSNPcaller
 ```
 b.Install perl modules in a directory in your current search paths. Find current paths with the following commanns:
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
 
 2. Running iSNPcaller

a. If this is the first run for a given project, initiate a project as follows:
```bash
perl iSNPcaller.pl
```
You will be prompted for the name of a new project (should not match any file/folder in the current working directory). Type the name of the project, hit return and confirm when prompted).

...alternatively you can initiate a new project by providing the project name at the command line:
```bash
perl iSNPcaller.pl myProject
```
confirm the project name when prompted.
