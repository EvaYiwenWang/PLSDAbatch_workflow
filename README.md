# PLSDAbatch
A new multivariate and non-parametric batch effect correction method based on Projection to Latent Structures Discriminant Analysis for microbiome data.

## pkgdown web page

The web page includes all the functions and vignettes within the package.

Link: https://evayiwenwang.github.io/PLSDAbatch/

## Installation

### Installation without vignettes

```r
devtools::install_github("https://github.com/EvaYiwenWang/PLSDAbatch")
```

### Installation with vignettes

First, we need to install the packages necessary for vignettes.

```r
# CRAN
cran.pkgs <- c('pheatmap', 'vegan', 'ruv', 'UpSetR', 'gplots', 
               'ggplot2', 'gridExtra', 'performance', 'BiocManager')
               
for(c in seq_len(length(cran.pkgs))){
if (!requireNamespace(cran.pkgs[c], quietly = TRUE))
    install.packages(cran.pkgs[c])
}
    
# Bioconductor
bioc.pkgs <- c('mixOmics', 'sva', 'limma', 'Biobase', 'metagenomeSeq')

for(b in seq_len(length(bioc.pkgs))){
if (!requireNamespace(bioc.pkgs[b], quietly = TRUE))
    BiocManager::install(bioc.pkgs[b])
}
```

Then, we are able to install *PLSDAbatch* with vignettes using the following command line. 


```r
devtools::install_github("https://github.com/EvaYiwenWang/PLSDAbatch", dependencies = T, build_vignettes = T)
```

## Reference

Wang, Y., & Lê Cao, K. A. (2020). A multivariate method to correct for batch effects in microbiome data. bioRxiv.
https://www.biorxiv.org/content/10.1101/2020.10.27.358283v1



