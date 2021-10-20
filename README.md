# Repository of code for analyzing genomes of tuskless and tusked African elephants.

Scripts for various analyses are within the corresponding directories. Many of these scripts simply involve running previously published software but using python scripts to parallelize tasks across scaffolds, in order to dramatically speed up analyses. We used Princeton University's high-performance computing cluster "Della", which uses SLURM to manage job submissions.

`survey_data_analyses` contains raw data and R markdown files used in analyzing the survey data. The subdirectory `population_survey_data_analysis` contains the analysis for Figure 1 in which we study the elephant population sizes before and after the war along with the frequency of tusklessness among females. The subdirectory `mother_offspring_data_analysis` contains the analysis of the mother offspring data in which we look at how often tuskless mothers give rise to tuskless and tusked daughters and also how often they give rise to sons and daughters overall. The R markdown files are heavily annotated, so please see these files for more information.

`HWE` contains a python script for parsing a VCF file and assessing deviations in heterozygosity within tuskless samples.

`LD` contains python scripts to run vcftools in a parallel fashion. These analyses quantify genotype linkage disequilibrium (LD) using vcftools (a separate job per scaffold to speed up analyses), and then an additional script parses this vcftools output to quantify how genotype LD decays with inter-SNP distance.

`RAiSD` contains python scripts to run the RAiSD program, parallelized across scaffold. We first split a genome-wide VCF in to scaffold-specific VCFs (see the `VCF_preprocessing` subdirectory), and then run RAiSD on each of these independently to speed up analyses.

`angsd` contains python scripts to run ANGSD software in a parallel fashion, by scaffold. There are separate subdirectories for the code to comput Fst and Dxy. Scripts have alphanumeric prefixes that indicate the order in which they were run.

`loxAfr4_genomeBuild` contains information about how scaffolds from the 3.0 genome build map to chromosomes in the 4.0 build. We mapped paired-end Illumina reads from our African elephant samples to the 3.0 genome build, since this version had annotaion information. However, we used the mapping of scaffolds to chromosomes for plotting purposes and to infer which scaffolds mapped to autosomes and which mapped to the X chromosome.

`overlap` contains scripts that were used for finding genomic intervals that had significant deviations in several genomic summary statistics. We did this by computing various statistics (Fst, Dxy, LD, Heterozygosity) for each locus (in a sliding window analysis), and we then overlapped the results using [pybedtools](https://daler.github.io/pybedtools/), a python package that implements bedtools. Further details may be found in the README in that subdirectory.


Also in this repository are other files necessary for reproducing our analyses. These include 
1. `DentalGenesWithEnsLafgIDs.txt`, a list of genes known to be involved in tooth development.
2. `position_2_gene.bed`, a list of genomic coordinates corresponding genic and intergenic regions (two genes separated by a hyphen '-'). This file was generated from our VCF that was annotated with [SNPeff](http://pcingola.github.io/SnpEff/) using the `loxAfr3.99` database.

