genome information was downloaded from:
```
ftp://ftp.broadinstitute.org:21/pub/assemblies/mammals/elephant/loxAfr4/Chromosomes.v2.map
```
from Palkopoulou et al 2018, PNAS. 

This is a FISH mapped assembly: FISH was conducted on ~200 BAC clones to anchor
85% of the LoxAfr3.0 assembly to chromosomes. FISH mapped assembly is named LoxAfr4.

For our study, we primarily used the LoxAfr3.0 assembly, which had annotation information.
Nonetheless, we mapped scaffolds to chromosomes using the LoxAfr4 assembly.

The script in this directory named `mapFile_2_table.py` was used to create a table that
indicated the mapping between LoxAfr3.0 scaffolds and LoxAfr4.0 chromosomes.

The file `ScaffoldsInChromosomes.txt` is the output of this script and was used in our study.

