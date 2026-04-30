#' Calculate Dirichlet Alpha Prior Parameters
#'
#' Computes information-weighted Dirichlet prior parameters for each cell type
#' by combining gene specificity (entropy-based), expression intensity,
#' spatial-reference similarity (weighted Pearson correlation), and optional
#' spatial smoothing via k-nearest neighbours.
#'
#' @param ref_matrix Normalized reference basis matrix (genes x cell types).
#' @param spatial_count Spatial gene expression count matrix (genes x spots).
#' @param spatial_coords Data frame or matrix of spot spatial coordinates
#'   (used for spatial smoothing). Default \code{NULL}.
#' @param method Character; currently only \code{"information_weighted"} is
#'   supported.
#' @param concentration Base concentration parameter for the Dirichlet prior.
#'   Default 2.0.
#' @param pseudocount Small value added for numerical stability. Default 1e-4.
#' @param min_alpha Minimum allowed alpha value per cell type. Default 1e-3.
#' @param normalize_spatial Logical; whether to CPM-normalize the spatial data.
#'   Default \code{TRUE}.
#' @param use_spatial_smoothing Logical; whether to apply spatial smoothing to
#'   the prior estimates. Default \code{TRUE}.
#' @param spatial_bandwidth Bandwidth for the Gaussian smoothing kernel.
#'   Default 0.2.
#'
#' @return Named numeric vector of alpha parameters (one per cell type).
#'
#' @export
calculateAlpha <- function(ref_matrix, spatial_count,
                           spatial_coords = NULL,
                           method = "information_weighted",
                           concentration = 2.0,
                           pseudocount = 1e-4,
                           min_alpha = 1e-3,
                           normalize_spatial = TRUE,
                           use_spatial_smoothing = TRUE,
                           spatial_bandwidth = 0.2) {

  message("=== Alpha Parameter Calculation (", method, ") ===")

  if (any(ref_matrix < 0)) stop("Reference matrix contains negative values")
  if (any(spatial_count < 0)) stop("Spatial count matrix contains negative values")

  ref_matrix    <- as.matrix(ref_matrix)
  spatial_count <- as.matrix(spatial_count)

  common_genes <- intersect(rownames(ref_matrix), rownames(spatial_count))
  if (length(common_genes) < 10)
    stop("Too few common genes (", length(common_genes), ")")
  message("   Common genes: ", length(common_genes))

  ref_sub <- ref_matrix[common_genes, , drop = FALSE]
  sp_sub  <- spatial_count[common_genes, , drop = FALSE]

  ref_norm <- sweep(ref_sub, 2, colSums(ref_sub) + 1e-10, "/")
  sp_cpm   <- sweep(sp_sub, 2, colSums(sp_sub) + 1e-10, "/") * 1e6
  sp_log   <- log2(sp_cpm + 1)

  if (method != "information_weighted")
    stop("Currently only 'information_weighted' method is supported")

  # 1. Gene specificity via entropy
  message("   Gene specificity weights...")
  ref_safe    <- pmax(ref_norm, 1e-10)
  max_entropy <- log2(ncol(ref_norm))
  entropy     <- -rowSums(ref_safe * log2(ref_safe))
  gene_spec   <- pmax((max_entropy - entropy) / max_entropy, 0)

  # 2. Cell type intensity
  ct_intensity <- colSums(ref_norm * gene_spec)

  # 3. Weighted Pearson correlation
  message("   Spatial-reference similarity...")
  w     <- gene_spec
  w_sum <- sum(w)

  ref_wm  <- colSums(ref_norm * w) / w_sum
  spot_wm <- colSums(sp_log * w) / w_sum

  sqrt_w       <- sqrt(w)
  ref_centered <- sweep(ref_norm, 2, ref_wm, "-") * sqrt_w
  sp_centered  <- sweep(sp_log, 2, spot_wm, "-") * sqrt_w

  ref_var  <- sqrt(colSums(ref_centered^2))
  sp_var   <- sqrt(colSums(sp_centered^2))

  ref_normed <- sweep(ref_centered, 2, pmax(ref_var, 1e-12), "/")
  sp_normed  <- sweep(sp_centered, 2, pmax(sp_var, 1e-12), "/")

  sim_matrix <- crossprod(sp_normed, ref_normed)
  sim_matrix <- pmax((sim_matrix + 1) / 2, 0.01)
  colnames(sim_matrix) <- colnames(ref_norm)
  rownames(sim_matrix) <- colnames(sp_sub)

  # 4. Combine
  int_w       <- ct_intensity / sum(ct_intensity)
  spot_priors <- sweep(sim_matrix, 2, int_w, "*")

  # 5. Spatial smoothing
  if (use_spatial_smoothing && !is.null(spatial_coords) &&
      nrow(spatial_coords) == nrow(spot_priors) && nrow(spatial_coords) > 1) {
    message("   Spatial smoothing...")
    k_smooth <- min(30, nrow(spatial_coords) - 1)
    nn_res   <- RANN::nn2(as.matrix(spatial_coords), k = k_smooth + 1)
    nn_idx   <- nn_res$nn.idx[, -1, drop = FALSE]
    nn_dist  <- nn_res$nn.dists[, -1, drop = FALSE]

    n_sp <- nrow(spatial_coords)
    sp_wt <- Matrix::sparseMatrix(
      i = rep(seq_len(n_sp), each = k_smooth),
      j = as.vector(t(nn_idx)),
      x = as.vector(exp(-(t(nn_dist))^2 / (2 * spatial_bandwidth^2))),
      dims = c(n_sp, n_sp)
    )
    rs     <- Matrix::rowSums(sp_wt) + 1e-10
    sp_wt  <- sp_wt / rs
    nb_avg <- as.matrix(sp_wt %*% spot_priors)
    spot_priors <- 0.7 * spot_priors + 0.3 * nb_avg
  }

  # 6. Global prior
  spot_priors_safe <- pmax(spot_priors, pseudocount)
  base_alpha <- exp(colMeans(log(spot_priors_safe)))

  # 7. Adaptive concentration
  prior_unc  <- apply(spot_priors, 2, stats::var)
  mean_unc   <- mean(prior_unc)
  adapt_conc <- concentration * (1 + mean_unc)
  message("   Adaptive concentration: ", round(adapt_conc, 3))

  # 8. Normalize
  base_alpha  <- pmax(base_alpha, min_alpha)
  base_alpha  <- base_alpha / sum(base_alpha)
  final_alpha <- base_alpha * adapt_conc
  names(final_alpha) <- colnames(ref_matrix)

  message("=== Alpha range: [", round(min(final_alpha), 4),
          ", ", round(max(final_alpha), 4), "] ===")
  return(final_alpha)
}
