
This vignette provides all reproducible codes for our article: 

# PLSDA-batch: a multivariate framework to correct for batch effects in microbiome data

Yiwen Wang, Kim-Anh Lê Cao

**Abstract:** Microbial communities are highly dynamic and sensitive to changes in the environment. Thus, microbiome data are highly susceptible to batch effects, defined as sources of unwanted variation that are not related to and obscure any factors of interest. Existing batch effect correction methods have been primarily developed for gene expression data. As such, they do not consider the inherent characteristics of microbiome data, including zero inflation, overdispersion and correlation between variables. We introduce new multivariate and non-parametric batch effect correction methods based on Partial Least Squares Discriminant Analysis (PLSDA). PLSDA-batch first estimates treatment and batch variation with latent components, then subtracts batch-associated components from the data. The resulting batch-effect-corrected data can then be input in any downstream statistical analysis. Two variants are proposed to handle unbalanced batch x treatment designs and to avoid overfitting when estimating the components via variable selection. We compare our approaches with popular methods managing batch effects, namely, removeBatchEffect, ComBat and Surrogate Variable Analysis, in simulated and three case studies using various visual and numerical assessments. We show that our three methods lead to competitive performance in removing batch variation while preserving treatment variation, especially for unbalanced batch × treatment designs. Our downstream analyses show selections of biologically relevant taxa. This work demonstrates that batch effect correction methods can improve microbiome research outputs. Reproducible code and vignettes are available on GitHub.

**Keywords:** microbiome data, multivariate, non-parametric, dimension reduction, batch effect correction

This article has been published and available as:

Wang, Y., & Lê Cao, K. A. (2023). PLSDA-batch: a multivariate framework to correct for batch effects in microbiome data. *Briefings in Bioinformatics*, 24(2), bbac622.

Link: <https://academic.oup.com/bib/article/24/2/bbac622/6991121>

We have also implemented a bookdown: <https://evayiwenwang.github.io/PLSDAbatch_workflow/>

The R package "PLSDAbatch" can be installed from <https://github.com/EvaYiwenWang/PLSDAbatch>.

