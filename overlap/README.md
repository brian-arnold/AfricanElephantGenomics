Since many analyses were done by scaffold, producing one file per scaffold in a specific directory, script `01_combine_tests.py` navigates into these directories and collects the results across scaffolds. This script then uses the genome-wide results for each statistic to compute a significance threshold, e.g. values that are above the 95% quantile of the genome-wide distribution. Significant windows/loci are printed out to intermediate files, in BED format, in a separate subdirectory named `outliers`.

Scripts `02_overlap_tests.py` and `02b_overlap_tests_forVennDiagram.py` then take the BED results in the `outliers` subdirectory and overlap them using [pybedtools](https://daler.github.io/pybedtools/), depositing results in the subdirectory `overlaps`. The latter 02b script does so in a way that produces intermediate files for making a Venn Diagram with [VennDiagram](https://cran.r-project.org/web/packages/VennDiagram/index.html) R package.

Finally, the script `03_overlap_dental_add_chr.py` takes selected results from the `overlaps` subdirectory and adds chromosome-level information (LoxAfr4 genome build), since many of the analyses were done by scaffold (LoxAfr 3.0 build, which has gene annotations), and prints out windows/loci that involve known tooth genes (based on information in the file `DentalGenesWithEnsLafgIDs.txt`).


