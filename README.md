# PlotCoverage: a method for visualization of read gene covarage in bulk transcriptomics data

PlotCoverage is R package, included functions calculation and visualization of gene coverage for specified genome coordinates in bulk RNA-seq data.
This method supports multigroup analysis.  


## Installation

Install the GitHub version:

```{r}
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("mamarkevi/plotCoverage")

```

Load package after installation:

```{r}
library(plotCoverage)
```
## Input
The input of PlotCoverage must be sorted bam-files of RNA-Seq read alignment and bam index-files. It allows to combine files into different groups to define differences of gene coverage between conditions.  