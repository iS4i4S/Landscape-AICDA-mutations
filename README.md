# Landscape of AICDA mutations
Code for reproducing analysis on the paper "Pan-cancer landscape of AID somatic related mutations and its potential role in the ICI response".
#Graphical abstract later--> ![alt text](https://github.com/...png "Hi there!")



### Instructions
Code in this repository can be used to recreate the essential parts of the main figures in the manuscript and some supplementary figures. Most of the code is in R (>=4.0.5), some other in python (>=3.0) or bash (PRIME program).

  * Install required R packages from CRAN and Bioconductor:
```r
## check for missing required packages, install them
required.packages <- c('data.table','ggplot2','cowplot','RColorBrewer',
                       'parallel','ggsignif','binom','scales','forestplot','ggpubr','survminer','sqldf','annotate','matrixStats','IHW','reshape2',
                       'cancereffectsizeR','ggrepel','Hmisc','Rcpp','pheatmap','ComplexHeatmap','PRIME','Lawstat','e1071','ggbeeswarm','stats','corrplot',
                       'ggforestplot','Meta','rmarkdown','karyoploteR','GSVA','survival','clusterProfiler','circlize')

new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)

## install packages from Bioconductor if not installed
if(!c('DESeq2',"TCGAutils","MutationalPatterns","maftools","dndscv","DoAbsolute","ABSOLUTE","Palimpsest","BiocOncoTK","seqinr","TxDb.Hsapiens.UCSC.hg19.knownGene","org.Hs.eg.db","genefilter","Biobase","DOSE","Seurat") %in% installed.packages()) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(c("DESeq2","TCGAutils","BiocOncoTK","MutationalPatterns","maftools","dndscv","DoAbsolute","ABSOLUTE","Palimpsest","seqinr","TxDb.Hsapiens.UCSC.hg19.knownGene","org.Hs.eg.db","genefilter","Biobase","DOSE","enrichplot","Seurat"))
}
 ```
  * Install MAESTRO (Model-based AnalysEs of Single-cell Transcriptome and RegulOme) by following the specifications on the original github page from [liulab-dfci](https://github.com/liulab-dfci/MAESTRO).
    
  * Install required python programs (this is for analyzing the single cell datasets), you can have extra information for installing scanpy [here](https://scanpy.readthedocs.io/en/stable/installation.html).
```bash
conda install seaborn scikit-learn statsmodels numba pytables
conda install -c conda-forge python-igraph leidenalg
#Alternately
pip install scanpy

```

Each folder has the title regarding each Figure number (e.g.Figure1) within the manuscript and is divided in 3 folders:
  * Data: contains necessary information to load into an R session in order to reproduce analyses. Some might contain a subfolder named "Raw_calculations" which is extra data needed either for the analyses for preculculated data included to simplify the coding.
  * Rmarkdowns: contains Rmd files with the needed code to reproduce the panels within the main figures and some supplemental figures.
  * r: contains manual written codes needed to reproduce the analyses.


### Data adquisition
Due to data regulations, maf, cnv, rna-seq and other data from the TCGA, ICGC, Composite cohorts cannot be included in the git but can be directly downloaded from their databases, you can consult our article with the links for all the datasets.

### Citation
- URL: 
- DOI: 

### Contact
E-mail any questions to [isaias.hernandez@icm-institute.org].
