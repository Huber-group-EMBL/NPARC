# NPARC
Nonparametric Analysis of Response Curves from Thermal Proteome Profiling Experiments.    
R package implementation of the method described by [Childs, Bach, Franken et al. (2019) Mol. Cell. Proteomics](https://doi.org/10.1074/mcp.TIR119.001481).

## Installation

### Installation from Bioconductor (recommended)

`NPARC` is in the process of being made available through Bioconductor! Therefore, the most reliable way to install it is via BiocManager:

```{R}
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("NPARC")
```

### Installation of the development version

```{R}
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(“Huber-group-EMBL/NPARC”)
```

## Getting started

The easiest way of learning how to use the `NPARC` package is to browse it's vignette:
```{R}
library(NPARC)
browseVignettes("NPARC")
```
