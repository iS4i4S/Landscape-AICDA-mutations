# Landscape of AICDA mutations
Code for reproducing analysis on the paper "Pan-cancer landscape of AID somatic related mutations and its potential role in the ICI response".
#Graphical abstract later--> ![alt text](https://github.com/...png "Hi there!")



### Instructions
Code in this repository can be used to recreate the essential parts of the main figures in the manuscript and most supplementary figures. All code is in R (>=4.0.5).

Install required R packages from CRAN and Bioconductor:
```r
## check for missing required packages, install them
required.packages <- c('data.table','ggplot2','cowplot','RColorBrewer',
                       'parallel','ggsignif','binom','scales',
                       'cancereffectsizeR','ggrepel','Hmisc','Rcpp','pheatmap','ComplexHeatmap','PRIME','Lawstat','e1071',
                       'ggforestplot','Meta','Seurat','MAESTRO','rmarkdown','karyoploteR','GSVA','survival','clusterProfiler','circlize')

new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)

## install packages from Bioconductor if not installed
if(!c('DESeq2',"MutationalPatterns","maftools","dndscv","DoAbsolute","ABSOLUTE","Palimpsest") %in% installed.packages()) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(c("DESeq2","MutationalPatterns","maftools","dndscv","DoAbsolute","ABSOLUTE","Palimpsest"))
}
```

Each folder has the title regarding each section/paragraph within the manuscript and is divided in 3 folders:
  * Data: contains necessary information to load into an R session in order to reproduce analyses.
  * Rmarkdowns: contains Rmd files with the needed code to reproduce the panels within the main figures and some supplemental figures.
  * r: contains manual written codes needed to reproduce the analyses.




### Citation
- URL: 
- DOI: 

### Contact
E-mail any questions to [isaias.hernandez@icm-institute.org].
