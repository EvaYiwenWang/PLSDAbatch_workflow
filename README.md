# PLSDAbatch
A new multivariate and non-parametric batch effect correction method based on Partial Least Squares Discriminant Analysis for microbial metagenome data.

## pkgdown web page

The web page includes all the functions and vignettes within the package.

Link: https://evayiwenwang.github.io/PLSDAbatch/

## Installation

### Installation without vignettes

```r
devtools::install_github("https://github.com/EvaYiwenWang/PLSDAbatch")
```
### Installation with vignettes

The listed CRAN and Bioconductor packages need to be installed first.
```r
# CRAN
cran.pkgs <- c('vegan', 'UpSetR', 'gplots', 'prettydoc', 'knitr', 'rmarkdown')
install.packages(cran.pkgs)

# Bioconductor
bioc.pkgs <- c('sva', 'limma')

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install(bioc.pkgs)
```
```r
devtools::install_github("https://github.com/EvaYiwenWang/PLSDAbatch", build_vignettes = T)
```

## Reference

Wang, Y., & Lê Cao, K. A. (2020). A multivariate method to correct for batch effects in microbiome data. bioRxiv.
https://www.biorxiv.org/content/10.1101/2020.10.27.358283v1



