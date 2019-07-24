# HLAEpitopes_R
HLA epitope association analysis


## Installation


### To update an existing installation

##### 1. Get source code
Simply download all .R files (git clone https://github.com/RonSchuyler/HLAEpitopes_R/) and move to existing R/ directory.

##### 2. Get allele sequence data
Download current allele sequences zip file from: 

ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Alignments_Rel_3370.zip

and rename to Alignments.zip

(Note that AlleleImport2.txt is no longer used.)


### For a new installation

##### Expected directory structure:

The base directory must contain the Alignments.zip file and subdirectories for input data, output data, and source code.

Data/ - Input data Affected/Control allele list goes here.

Output/ - Results are found here.

R/ - All source code goes here.

Alignments.zip - Download current allele sequences zip file here: 
ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/
For example: Alignments_Rel_3370.zip and rename to Alignments.zip



##### Install required packages:
In an R session, type:
```
install.packages("gWidgets")
install.packages("abind")
install.packages("tools")

```

for Windows installations only, also install gWidgetsGTK package:
```
install.packages("gWidgetsRGtk2")
```



