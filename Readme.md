# **FunRes** 

## FunRes is a computational method designed for the identification of functional cell states. FunRes utilises tissue single cell RNA-seq data to reconstruct the functional cell-cell communication network which is leveraged for partitioning each cell type into functional states.

## Citation

This tool was pas published in the following article:
- Sascha Jung, Kartikeya Singh, Antonio del Sol, FunRes: resolving tissue-specific functional cell states based on a cell–cell communication network model, Briefings in Bioinformatics, Volume 22, Issue 4, July 2021, bbaa283, [DOI](https://doi.org/10.1093/bib/bbaa283)

## Quick start

## Installation

 1) Instalar Bioconductor (si aún no lo tienes)
```r
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
```
 2) Instalar las dependencias de Bioconductor
```r
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db", "RSpectra"))
```

 3) Install CRAN  dependencies
```r
install.packages(c("data.table","doParallel","dplyr","ggplot2","gtools",
  "igraph","Matrix","pheatmap","plyr","reshape2","rlist",
  "Seurat","snow","stringr","taRifx","textshape"))
```
 4) Instalar funres desde GitHub
```r
if (!requireNamespace("remotes", quietly=TRUE))
  install.packages("remotes")
remotes::install_github("aitoe96/funres")
```




#### Authors

_Modified and packaged by Aitor Martínez Pérez._ 
- [Aitor Martínez Pérez](https://www.cicbiogune.es/people/amperez)

  FunRes was developed in the [Computational Biology Group](https://wwwen.uni.lu/lcsb/research/computational_biology) by

- [Kartikeya Singh](https://wwwen.uni.lu/lcsb/people/kartikeya_singh)
- [Sascha Jung](https://www.cicbiogune.es/people/sjung)
- [Antonio del Sol](https://wwwfr.uni.lu/lcsb/people/antonio_del_sol_mesa)
