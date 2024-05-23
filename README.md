post GWAS analysis for PTSD project
=============================================

Step 0:data collection.
------------------------
The two cohorts were obtained from the [PGC](https://pgc.unc.edu/for-researchers/download-results/) and [FINNGEN](https://finngen.gitbook.io/documentation/v/r8/data-description) databases.

Step 1: meta analysis
------------------------
Tool: [Metal](https://csg.sph.umich.edu/abecasis/metal/).
References：[METAL: fast and efficient meta-analysis of genomewide association scans](https://pubmed.ncbi.nlm.nih.gov/20616382/).

The METAL software is designed to facilitate meta-analysis of large datasets (such as several whole genome scans) in a convenient, rapid and memory efficient manner. 
```
metal metal.txt
```

The metal.txt file shows the specific parameters. SCHEME STDERR, weights effect size estimates using the inverse of the corresponding standard errors. AVERAGEFREQ ON, tracks the effect allele frequency across all files and report the mean effect allele frequency. TRACKPOSITIONS ON, track SNP genome positions. GENOMICCONTROL ON, automatically correct test statistics to account for small size. automatically correct test statistics to account for small amounts of population stratification or unaccounted for relatedness.

```
SCHEME STDERR
AVERAGEFREQ ON
TRACKPOSITIONS ON
GENOMICCONTROL ON
```
