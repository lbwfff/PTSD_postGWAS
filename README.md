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
The following segment defines the parameters for processing genetic marker data from two datasets, which in this case are the FINNGEN and PGC PTSD data. It specifies marker identifiers, alleles, effect sizes, standard errors, p-values, allele frequencies, chromosome locations, and genomic positions for each dataset. Each column was filled in with its corresponding column in the database. Following data processing using these parameters, a meta-analysis is executed to synthesize findings across both datasets.

```
# Describe and process the SardiNIA input files
MARKER   rsids
ALLELE   ref alt
EFFECT   beta
STDERR   sebeta
PVAL     pval
FREQLABEL af_alt
CHROMOSOME chrom
POSITION pos

PROCESS FINNGEN

# Describe and process the SardiNIA input files
MARKER   SNP
ALLELE   A1 A2
EFFECT   log(OR)
STDERR   SE
PVAL     P
FREQLABEL FRQ_A_
CHROMOSOME CHR
POSITION BP

PROCESS PGCPTSD

# Execute meta-analysis
ANALYZE
```

Step 2: Post Meta-analysis GWAS
------------------------
Tool: [gwaslab](https://github.com/Cloufield/gwaslab/).
References：[GWASlab](https://cloufield.github.io/gwaslab/).

gwaslab is a comprehensive tool designed for Genome-Wide Association Studies (GWAS) to perform statistical analyses, visualize genomic data in manhattan and qq plots, and interpret genetic associations.

```
import gwaslab as gl
```
To conduct a GWAS using the meta-analysis file, ensure that its columns correspond precisely to the parameters specified in gwaslab.
```
mysumstats =gl.Sumstats("METAANALYSIS1.TBL",
                        rsid="MarkerName",
                         ea="Allele1",
                         chrom='Chromosome',
                         pos="Position",
                         nea="Allele2",                        
                         beta="Effect",
                         se="StdErr",
                         p="P-value",
                         direction="Direction",
                         build="38")
```
Complete tutorial for next steps is found on https://cloufield.github.io/gwaslab/tutorial_3.4/.
The p-value threshold used is 1e-6.


Step X: TwoSampleMR and colocalisation analysis
------------------------
MR_and_coloc.R allows twosample MR as well as colocalization analysis based on SMR output.
Additional information on the MAF value of the outcome data is required. all_outcome is a data frame with two columns: rsids and af_alt.

```
all_outcome<-data.table::fread('finngen/finngen_PTSD')
all_outcome<-data.frame(rsids=c(all_outcome$rsids),
                        af_alt=c(all_outcome$af_alt))

source('MR_and_coloc.R')

result <- MR_and_coloc('./testExi.ENSG00000182481.10.txt','PTSD','./test',
             2860,400000,37,all_outcome,'./clump/1kg.v3/EUR')
```

Parameters to be provided by MR_and_coloc include: 
smrfile: SMR output, trait_name: outcome name, e.g. "PTSD", plotpath: folder to save plots, n_qtl: number of samples of qtl data, n_ gwas: number of samples in gwas data, g_version: gene version (37 or 38) of qtl data. 


