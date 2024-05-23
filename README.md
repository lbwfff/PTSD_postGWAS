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

