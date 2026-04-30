# SDDGP 
*A spatial Bayesian MCMC cell type deconvolution method for spatial transcriptomics*  
![](https://github.com/njjyxl/SDDGP/blob/devel/Overview.jpg)  
SDDGP leverages the hierarchical non-linear mapping capabilities of DGP to capture complex spatial dependencies while maintaining probabilistic interpretability through Markov Chain Monte Carlo (MCMC) inference.   
## Installation  
You can install SDDGP on Github with the following code.  
### Dependencies  
+ R version >= 4.3.1.  
+ R packages: Matrix, Rcpp,RANN, methods,nnls,fields,sf,concaveman,ggplot2.
```  
## Installation
install.packages('devtools')  

# install the Spatialsmooth package  
devtools::install_github('njjyxl/SDDGP')  

# load package  
library(SDDGP)  
```
The R package has been installed successfully on Operating systems:  
+ CentOS Linux release 7.5.1804 (Core)  
+ Windows 11  

## Issues  
All feedback, bug reports and suggestions are warmly welcomed! Please make sure to raise issues with a detailed and reproducible exmple and also please provide the output of your sessionInfo() in R!  

## Quick start

```r
library(SDDGP)

# 1. Build reference
ref <- buildReferenceMatrix(sc_count, spatial_count, cell_meta,
                            ct_varname = "celltype")

# 2. Compute alpha prior
alpha <- calculateAlpha(ref$basis, sc_count, spatial_coords,
                        concentration = 20)

# 3. Spatial kernel (optional)
ED <- fields::rdist(as.matrix(spatial_coords[, c("x","y")]))
kernel_mat <- maternKernel(ED, range = 5, smoothness = 30)
diag(kernel_mat) <- 0

# 4. Stage 1 deconvolution
stage1 <- runDeconvolution(spatial_count, ref, spatial_coords, alpha,
                           kernel_mat, mcmc_iterations = 3000,
                           early_stop = TRUE)

# 5. Resolution enhancement
enhanced <- enhanceResolution(stage1, ref, spatial_coords,
                              num_grids = 51000)

# 6. Stage 2 refinement
stage2 <- runStage2Refinement(stage1, ref, spatial_count_full)

# Diagnostics
printConvergenceReport(stage1)
plotConvergenceDiagnostics(stage1, spatial_coords)

# 7. Extract results
prop_matrix <- results$theta_estimates

colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
    "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
    "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
SDDGP.pie(prop_matrix, pos,colors = colors ,radius = 0.52)
SDDGP.celllandscape(prop_matrix, pos)
SDDGP.Cor(prop_matrix)
SDDGP.cellabundance(prop_matrix,pos)
```
