# Landscape of AID-related mutations
Code for reproducing analysis on the paper "Pan-cancer landscape of AID-related mutations, composite mutations and their potential role in the ICI response".
![alt text](https://github.com/iS4i4S/Landscape-AICDA-mutations/blob/main/Data/Pan-cancer%20landscape%20of%20AID%20mutations.jpeg "Hi there!")

### Abstract

Activation-induced cytidine deaminase, AICDA or AID, is a driver of somatic hypermutation and class-switch recombination in immunoglobulins. In addition, this deaminase belonging to the APOBEC family may have off-target effects genome-wide, but its effects at pan-cancer level are not well elucidated. Here, we used different pan-cancer datasets, totaling more than 50,000 samples analyzed by whole-genome, whole-exome or targeted sequencing. AID mutations are present at pan-cancer level with higher frequency in hematological cancers and higher presence at transcriptionally active TAD domains. AID synergizes initial hotspot mutations by a second composite mutation. AID mutational load was found to be independently associated with a favorable outcome in immune-checkpoint inhibitors (ICI) treated patients across cancers after analyzing 2,000 samples. Finally, we found that AID-related neoepitopes, resulting from mutations at more frequent hotspots if compared to other mutational signatures, enhance CXCL13/CCR5 expression, immunogenicity, and T-cell exhaustion, which may increase ICI sensitivity.

### Citation
If you use any data or code derived from this study, please cite:

- .......(missing citation) 
- DOI: [](insertLink later)


### Instructions
Code in this repository can be used to recreate the essential parts of the main figures in the manuscript and some supplementary figures. Most of the code is in R (>=4.0.5), some other in python (>=3.0) or bash (PRIME program).

  * Install required R packages from CRAN and Bioconductor:
```r
## check for missing required packages, install them
required.packages <- c('data.table','ggplot2','cowplot','RColorBrewer','Seurat',
                       'parallel','ggsignif','binom','scales','forestplot','ggpubr','survminer','sqldf','annotate','matrixStats','IHW','reshape2',
                       'cancereffectsizeR','ggrepel','Hmisc','Rcpp','pheatmap','ComplexHeatmap','PRIME','Lawstat','e1071','ggbeeswarm','stats','corrplot',
                       'ggforestplot','Meta','rmarkdown','karyoploteR','GSVA','survival','clusterProfiler','circlize')

new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)

## install packages from Bioconductor if not installed
if(!c('DESeq2',"TCGAutils","MutationalPatterns","maftools","dndscv","Palimpsest","BiocOncoTK","seqinr","TxDb.Hsapiens.UCSC.hg19.knownGene","org.Hs.eg.db","genefilter","Biobase","DOSE") %in% installed.packages()) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(c("DESeq2","TCGAutils","BiocOncoTK","MutationalPatterns","maftools","dndscv","Palimpsest","seqinr","TxDb.Hsapiens.UCSC.hg19.knownGene","org.Hs.eg.db","genefilter","Biobase","DOSE","enrichplot"))
}
 ```
  * Install ABSOLUTE by logging into Broad cga website: https://software.broadinstitute.org/cancer/cga then download and install the file "ABSOLUTE_1.0.6.tar.gz" Then:
 ```r
  install.packages("ABSOLUTE_1.0.6.tar.gz", repos = NULL, type = "source")
 ```

  * Install DoAbsolute using devtools: 
  ```r
  devtools::install_github("ShixiangWang/DoAbsolute")
```

Each folder has the title regarding each Figure number (e.g.Figure1) within the manuscript and is divided in 3 folders:
  * Data: contains necessary information to load into an R session in order to reproduce analyses. Some might contain a subfolder named "Raw_calculations" which is extra data needed for the analyses.
  * Rmarkdowns: contains Rmd files with the needed code to reproduce the panels within the main figures and some supplemental figures.
  * r: contains manual written codes needed to reproduce the analyses. **This section contains the Rdat objects containing the functions to calculate and extract AID and APOBEC realted mutations**. The file is named [AICDA_enrich_function_WES_WGS.RDat](https://github.com/iS4i4S/Landscape-AICDA-mutations/blob/main/r/AICDA_enrich_function_WES_WGS.RData).  
  

### HTMLs
Visualize HTML files containing graphics for each main figure by clicking on the coresponding links here:

 * [Figure1.html](http://htmlpreview.github.io/?https://github.com/iS4i4S/Landscape-AICDA-mutations/blob/main/Rmarkdowns/Figure1.html).
 * Figure 2a and Figure 2b: Code for graphs was generated in python using the [Figure2ab.py](https://github.com/iS4i4S/Landscape-AICDA-mutations/blob/main/Rmarkdowns/MutationAggregate.py) which uses as input data generated from Figure2c.html (section Figure 2ab). Additional information can be found through this [github page from Akdemir Lab](https://github.com/akdemirlab/MutationalDistribution).
 * [Figure2c.html](http://htmlpreview.github.io/?https://github.com/iS4i4S/Landscape-AICDA-mutations/blob/main/Rmarkdowns/Figure2.html).
 * [Figure3.html](http://htmlpreview.github.io/?https://github.com/iS4i4S/Landscape-AICDA-mutations/blob/main/Rmarkdowns/Figure3.html).
 * [Figure4.html](http://htmlpreview.github.io/?https://github.com/iS4i4S/Landscape-AICDA-mutations/blob/main/Rmarkdowns/Figure4.html).
 * [Figure5.html](http://htmlpreview.github.io/?https://github.com/iS4i4S/Landscape-AICDA-mutations/blob/main/Rmarkdowns/Figure5.html).

### Data adquisition
Due to data regulations, maf, cnv, rna-seq and other data from the TCGA, ICGC, Composite cohorts cannot be included in the git but can be directly downloaded from their databases, you can consult our article with the links for all the datasets.

### Tracking AID-related mutations

As stated in the article, we developed a code to detect c-AID related mutations over *wrCy/rGyw* (+/- strand, where "W" stands to either adenine or thymine, "R" to purine, and "Y" to pyrimidine) motifs, giving a total of 8 motifs per strand (positive-strand = AACC, AACT, AGCC, AGCT, TACC, TACT, TGCC, TGCT; negative-strand = TTGG, TTGA, TCGG, TCGA, ATGG, ATGA, ACGG, ACGA).

The code was developed under R, takes a maf (mutation annotation format) object as input and outputs an S3 class object containing: 
- i) A matrix of the 768 possible tetranucleotide substitutions across the samples
- ii) A data table with all the needed values for enrichment calculation, the enrichment score, Fisher exact-test p-value and FDR for enrichment, the fraction of AID mutations, among others. 
- iii) A maf-like data table, with the same format as the input, containing only the attributed AID mutations. 

You can visualize an html example of using this function [here](link..missing).


### Contact
E-mail any questions to [isaias.hernandez@icm-institute.org].
