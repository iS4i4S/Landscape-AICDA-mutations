# Landscape of AID-related mutations
Code for reproducing analysis on the paper "Pan-cancer landscape of AID somatic related mutations and its potential role in the ICI response".
![alt text](https://github.com/iS4i4S/Landscape-AICDA-mutations/blob/main/Data/Pan-cancer%20landscape%20of%20AID%20mutations.jpeg "Hi there!")



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
  * Install MAESTRO (Model-based AnalysEs of Single-cell Transcriptome and RegulOme) by following the specifications on the original github page from [liulab-dfci](https://github.com/liulab-dfci/MAESTRO).
  * Install ABSOLUTE by logging into Broad cga website: https://software.broadinstitute.org/cancer/cga then download and install the file "ABSOLUTE_1.0.6.tar.gz" Then:
 ```r
  install.packages("ABSOLUTE_1.0.6.tar.gz", repos = NULL, type = "source")
 ```

  * Install DoAbsolute using devtools: 
  ```r
  devtools::install_github("ShixiangWang/DoAbsolute")
```

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

### HTMLs
Generate HTML files containing graphics for each main figure from the command line, which can be viewed in your web browser:

```bash
R -e "rmarkdown::render('Rmarkdowns/Figure1.Rmd', output_file = 'figure-1.html')"
R -e "rmarkdown::render('Rmarkdowns/Figure2.Rmd', output_file = 'figure-2.html')"
R -e "rmarkdown::render('Rmarkdowns/Figure3.Rmd', output_file = 'figure-3.html')"
R -e "rmarkdown::render('Rmarkdowns/Figure4.Rmd', output_file = 'figure-4.html')"
R -e "rmarkdown::render('Rmarkdowns/Figure5.Rmd', output_file = 'figure-5.html')"
R -e "rmarkdown::render('Rmarkdowns/Figure6.Rmd', output_file = 'figure-6.html')"
R -e "rmarkdown::render('Rmarkdowns/Figure7.Rmd', output_file = 'figure-7.html')"
```
Alternatively, you can visualize the files by clicking on the coresponding links here:
 * [Figure1.html](http://htmlpreview.github.io/?https://github.com/iS4i4S/Landscape-AICDA-mutations/blob/main/Rmarkdowns/Figure1.html).
 * [Figure2.html](http://htmlpreview.github.io/?https://github.com/iS4i4S/Landscape-AICDA-mutations/blob/main/Rmarkdowns/Figure2.html).
 * [Figure3.html](http://htmlpreview.github.io/?https://github.com/iS4i4S/Landscape-AICDA-mutations/blob/main/Rmarkdowns/Figure3.html).
 * [Figure4.html](http://htmlpreview.github.io/?https://github.com/iS4i4S/Landscape-AICDA-mutations/blob/main/Rmarkdowns/Figure4.html).
 * [Figure5.html](http://htmlpreview.github.io/?https://github.com/iS4i4S/Landscape-AICDA-mutations/blob/main/Rmarkdowns/Figure5.html).
 * [Figure6.html](http://htmlpreview.github.io/?https://github.com/iS4i4S/Landscape-AICDA-mutations/blob/main/Rmarkdowns/Figure6.html).
 * [Figure7C.html](http://htmlpreview.github.io/?https://github.com/iS4i4S/Landscape-AICDA-mutations/blob/main/Rmarkdowns/Figure7C.html).
 * [Figure7.html](http://htmlpreview.github.io/?https://github.com/iS4i4S/Landscape-AICDA-mutations/blob/main/Rmarkdowns/Figure7.html).

### Data adquisition
Due to data regulations, maf, cnv, rna-seq and other data from the TCGA, ICGC, Composite cohorts cannot be included in the git but can be directly downloaded from their databases, you can consult our article with the links for all the datasets.

### Citation
- URL: 
- DOI: 

### Contact
E-mail any questions to [isaias.hernandez@icm-institute.org].
