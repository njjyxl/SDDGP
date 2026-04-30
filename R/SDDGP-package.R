#' SDDGP: Spatial Deconvolution via Deep Gaussian Processes
#'
#' A two-stage Bayesian framework for spatial transcriptomics
#' deconvolution with full posterior uncertainty quantification.
#'
#' @useDynLib SDDGP, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Matrix rowSums colSums
#' @importFrom RANN nn2
#' @importFrom fields rdist
#' @importFrom sf st_as_sf st_bbox st_area st_intersects st_geometry
#' @importFrom concaveman concaveman
#' @importFrom nnls nnls
#' @importFrom parallel makeCluster stopCluster parLapplyLB detectCores
#'   clusterExport clusterEvalQ
#' @importFrom ggplot2 ggplot aes geom_point geom_col geom_smooth
#'   scale_color_viridis_c labs theme_minimal coord_fixed theme .data
#' @importFrom stats kmeans rnorm quantile sd var setNames median scale loess
#' @importFrom methods is as
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics hist abline legend plot par
#' @importFrom utils tail
"_PACKAGE"
