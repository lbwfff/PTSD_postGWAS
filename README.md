post GWAS analysis for PTSD project
=============================================

Step 0:data collection.
------------------------
The two cohorts were obtained from the [PGC](https://pgc.unc.edu/for-researchers/download-results/) and [FINNGEN](https://finngen.gitbook.io/documentation/v/r8/data-description) databases.

Adjusting the data structure, effect allele frequency was manually calculated for PGC data. The formula is: (Frequency of coded allele in cases * Number of cases + Frequency of coded allele in controls * Number of controls)/(Number of cases + Number of controls).
Example:
```
awk 'BEGIN {FS=OFS="\t"} NR==1 {header=$0; print $0; next} {Nca=$17; FRQ_A_4363=$6; Nco=$18; FRQ_U_10976=$7; if ((Nca + Nco) != 0) {af_alt = (Nca * FRQ_A_4363 + Nco * FRQ_U_10976) / (Nca + Nco)} else {af_alt = 0}; print $0, af_alt}' pts_lat_freeze2_overall.results > pts_lat_freeze2_overall.results.adj
sed -i '1s/$/\taf_alt/' pts_lat_freeze2_overall.results.adj
```
Simple adjustment of column names for FINNGEN data.

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

Step 2:
------------------------
