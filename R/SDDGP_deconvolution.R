#' @title Matérn Kernel Function
#' @description Computes a Matérn covariance matrix for a given distance matrix.
#' @param D Distance matrix.
#' @param range Range (length-scale) parameter.
#' @param smoothness Smoothness parameter (0.5 = exponential, Inf = squared
#'   exponential).
#' @return Covariance matrix of the same dimension as \code{D}.
#' @export
maternKernel <- function(D, range, smoothness) {
  if (smoothness == 0.5) {
    return(exp(-D / range))
  } else if (smoothness == Inf) {
    return(exp(-(D^2) / (2 * range^2)))
  } else {
    scaled <- sqrt(2 * smoothness) * D / range
    part1  <- 2^(1 - smoothness) / gamma(smoothness)
    part2  <- scaled^smoothness * besselK(scaled, smoothness)
    return(part1 * part2)
  }
}

# ---------- RBF kernel helpers ----------

#' @keywords internal
rbfKernelFast <- function(X1, X2, lengthscales, variance, jitter = 1e-6) {
  if (length(lengthscales) == 1)
    lengthscales <- rep(lengthscales, ncol(X1))

  X1s <- sweep(X1, 2, sqrt(lengthscales), "/")
  X2s <- sweep(X2, 2, sqrt(lengthscales), "/")

  sq_dist <- if (identical(X1, X2)) as.matrix(stats::dist(X1s)^2)
             else fields::rdist(X1s, X2s)^2

  K <- variance * exp(-0.5 * sq_dist)
  if (identical(X1, X2)) diag(K) <- diag(K) + jitter
  K
}

#' @keywords internal
rbfKernelSpatial <- function(X1, X2, lengthscales, variance,
                             kernel_mat = NULL, spatial_weight = 0.3,
                             X1_indices = NULL, X2_indices = NULL,
                             jitter = 1e-6) {
  if (length(lengthscales) == 1)
    lengthscales <- rep(lengthscales, ncol(X1))

  X1s <- sweep(X1, 2, sqrt(lengthscales), "/")
  X2s <- sweep(X2, 2, sqrt(lengthscales), "/")

  sq_dist <- if (identical(X1, X2)) as.matrix(stats::dist(X1s)^2)
             else fields::rdist(X1s, X2s)^2

  K_rbf <- variance * exp(-0.5 * sq_dist)

  if (!is.null(kernel_mat) && !is.null(X1_indices) && !is.null(X2_indices)) {
    vx1 <- X1_indices[X1_indices <= nrow(kernel_mat)]
    vx2 <- X2_indices[X2_indices <= ncol(kernel_mat)]
    if (length(vx1) == nrow(X1) && length(vx2) == nrow(X2)) {
      K_sp <- kernel_mat[vx1, vx2]
      if (nrow(K_sp) == nrow(K_rbf) && ncol(K_sp) == ncol(K_rbf))
        K_rbf <- (1 - spatial_weight) * K_rbf + spatial_weight * K_sp * variance
    }
  }

  if (identical(X1, X2)) diag(K_rbf) <- diag(K_rbf) + jitter
  K_rbf
}

#' @keywords internal
safeCholesky <- function(A, jitter = 1e-6) {
  current <- jitter
  for (i in 1:10) {
    result <- tryCatch(t(chol(A + diag(nrow(A)) * current)),
                       error = function(e) NULL)
    if (!is.null(result)) return(result)
    current <- current * 10
  }
  warning("Cholesky failed; falling back to eigendecomposition")
  eig <- eigen(A, symmetric = TRUE)
  eig$values[eig$values < 1e-12] <- 1e-12
  eig$vectors %*% diag(sqrt(eig$values))
}

# ---------- DGP layer ----------

#' @keywords internal
createDGPLayer <- function(X, n_inducing = 50, lengthscale_init = 1.0,
                           variance_init = 1.0, layer_depth = 1,
                           output_dim = NULL, save_inducing_indices = FALSE) {
  n <- nrow(X); d <- ncol(X)
  if (is.null(output_dim)) output_dim <- d

  inducing_indices <- NULL
  if (n_inducing < n) {
    tryCatch({
      km <- stats::kmeans(X, centers = min(n_inducing, n), nstart = 5,
                          iter.max = 50)
      inducing_pts <- km$centers
      if (save_inducing_indices) {
        nn_res <- RANN::nn2(X, inducing_pts, k = 1)
        inducing_indices <- nn_res$nn.idx[, 1]
      }
    }, error = function(e) {
      sel <- sample(n, min(n_inducing, n))
      inducing_pts <<- X[sel, , drop = FALSE]
      if (save_inducing_indices) inducing_indices <<- sel
    })
  } else {
    sel <- sample(n, min(n_inducing, n))
    inducing_pts <- X[sel, , drop = FALSE]
    if (save_inducing_indices) inducing_indices <- sel
  }
  inducing_pts <- as.matrix(inducing_pts)

  var_means <- lapply(seq_len(output_dim),
                      function(i) matrix(stats::rnorm(nrow(inducing_pts)), ncol = 1))
  var_chols <- lapply(seq_len(output_dim),
                      function(i) diag(nrow(inducing_pts)) * 0.1)

  res <- list(
    inducing_points   = inducing_pts,
    variational_means = var_means,
    variational_chols = var_chols,
    log_lengthscales  = rep(log(lengthscale_init), d),
    log_variances     = rep(log(variance_init), output_dim),
    log_noises        = rep(log(0.01), output_dim),
    layer_depth       = layer_depth,
    input_dim         = d,
    output_dim        = output_dim,
    n_inducing        = nrow(inducing_pts)
  )
  if (save_inducing_indices) res$inducing_indices <- inducing_indices
  res
}

# ---------- DGP forward ----------

#' @keywords internal
dgpForward <- function(layers, X, kernel_mat = NULL, spatial_indices = NULL,
                       spatial_weight = 0.3, n_samples = 5) {
  n_layers <- length(layers)
  current_input   <- X
  current_indices <- spatial_indices
  layer_outputs   <- vector("list", n_layers)

  for (l in seq_len(n_layers)) {
    layer <- layers[[l]]
    use_spatial <- (l == 1 && !is.null(kernel_mat) && !is.null(current_indices))

    if (use_spatial) {
      if (!is.null(layer$inducing_indices)) {
        ind_idx <- layer$inducing_indices
      } else {
        nn_res  <- RANN::nn2(current_input, layer$inducing_points, k = 1)
        ind_idx <- current_indices[nn_res$nn.idx[, 1]]
      }
      K_uu <- rbfKernelSpatial(layer$inducing_points, layer$inducing_points,
                               exp(layer$log_lengthscales), exp(layer$log_variances[1]),
                               kernel_mat, spatial_weight, ind_idx, ind_idx)
      K_uf <- rbfKernelSpatial(layer$inducing_points, current_input,
                               exp(layer$log_lengthscales), exp(layer$log_variances[1]),
                               kernel_mat, spatial_weight, ind_idx, current_indices)
    } else {
      K_uu <- rbfKernelFast(layer$inducing_points, layer$inducing_points,
                            exp(layer$log_lengthscales), exp(layer$log_variances[1]))
      K_uf <- rbfKernelFast(layer$inducing_points, current_input,
                            exp(layer$log_lengthscales), exp(layer$log_variances[1]))
    }

    L_uu     <- safeCholesky(K_uu)
    A        <- solve(L_uu, K_uf)
    L_uu_inv <- solve(L_uu)
    At       <- t(A)
    A_sq_cs  <- colSums(A^2)

    out_means <- matrix(0, nrow(current_input), layer$output_dim)
    out_vars  <- matrix(0, nrow(current_input), layer$output_dim)

    for (d in seq_len(layer$output_dim)) {
      S   <- layer$variational_chols[[d]] %*% t(layer$variational_chols[[d]])
      mu  <- At %*% (L_uu_inv %*% layer$variational_means[[d]])
      K_diag <- rep(exp(layer$log_variances[d]), nrow(current_input))
      M   <- L_uu_inv %*% S %*% t(L_uu_inv)
      vf  <- K_diag - A_sq_cs + colSums(A * (M %*% A))
      out_means[, d] <- as.vector(mu)
      out_vars[, d]  <- pmax(vf, 1e-8)
    }

    if (l < n_layers) {
      samps <- array(0, dim = c(nrow(current_input), layer$output_dim, n_samples))
      for (s in seq_len(n_samples))
        for (d in seq_len(layer$output_dim))
          samps[, d, s] <- out_means[, d] + sqrt(out_vars[, d]) *
                             stats::rnorm(nrow(current_input))
      current_input   <- apply(samps, c(1, 2), mean)
      current_indices <- NULL
      layer_outputs[[l]] <- list(mean = out_means, variance = out_vars,
                                 samples = samps)
    } else {
      layer_outputs[[l]] <- list(mean = out_means, variance = out_vars)
    }
  }
  layer_outputs
}

# ---------- Public DGP API ----------

#' Build a Spatial Deep Gaussian Process Model
#'
#' Constructs a multi-layer DGP model for learning spatial cell-type priors.
#' The first layer optionally incorporates an external spatial kernel (e.g.
#' Matérn) to encode prior knowledge of tissue-level spatial correlation.
#'
#' @param spatial_coords Matrix or data frame of spot coordinates.
#' @param ref_matrix Reference matrix list as returned by
#'   \code{\link{buildReferenceMatrix}}.
#' @param kernel_mat Optional pre-computed spatial kernel matrix.
#' @param n_layers Number of DGP layers. Default 3.
#' @param n_inducing_per_layer Number of inducing points per layer. Default 30.
#'
#' @return A list representing the DGP model.
#' @export
buildSpatialDGP <- function(spatial_coords, ref_matrix, kernel_mat = NULL,
                            n_layers = 3, n_inducing_per_layer = 30) {
  norm_coords <- scale(spatial_coords)
  n_ct <- ncol(ref_matrix$basis)
  hidden_dim <- min(10, max(5, n_ct))

  layers <- vector("list", n_layers)
  layers[[1]] <- createDGPLayer(
    norm_coords, n_inducing_per_layer, 0.5, 1.0, 1, hidden_dim,
    save_inducing_indices = !is.null(kernel_mat))

  if (n_layers > 2) {
    for (i in 2:(n_layers - 1)) {
      dummy <- matrix(stats::rnorm(nrow(norm_coords) * hidden_dim),
                      ncol = hidden_dim)
      layers[[i]] <- createDGPLayer(dummy, n_inducing_per_layer, 1.0, 1.0, i,
                                    hidden_dim, FALSE)
    }
  }

  final_in <- matrix(stats::rnorm(nrow(norm_coords) * hidden_dim),
                     ncol = hidden_dim)
  layers[[n_layers]] <- createDGPLayer(
    final_in, min(n_inducing_per_layer, n_ct * 3), 1.0, 1.0, n_layers, n_ct,
    FALSE)

  list(layers = layers, n_layers = n_layers, spatial_coords = norm_coords,
       n_celltypes = n_ct, hidden_dim = hidden_dim,
       coord_center = attr(norm_coords, "scaled:center"),
       coord_scale  = attr(norm_coords, "scaled:scale"),
       kernel_mat   = kernel_mat)
}

#' Predict Cell-Type Priors from a DGP Model
#'
#' @param dgp_model DGP model as returned by \code{\link{buildSpatialDGP}}.
#' @param spot_coords Spot coordinates (same scale as original training data).
#' @param spatial_weight Weight of the spatial kernel component (0-1).
#'   Default 0.3.
#'
#' @return A list with \code{priors} (spots x cell types) and
#'   \code{uncertainty} (per-spot scalar).
#' @export
predictCelltypePriors <- function(dgp_model, spot_coords,
                                  spatial_weight = 0.3) {
  if (is.null(dgp_model$coord_center) || is.null(dgp_model$coord_scale)) {
    nc <- scale(spot_coords)
  } else {
    nc <- scale(spot_coords, center = dgp_model$coord_center,
                scale = dgp_model$coord_scale)
  }

  outputs <- dgpForward(dgp_model$layers, nc, dgp_model$kernel_mat,
                        seq_len(nrow(spot_coords)), spatial_weight, 3)
  final   <- outputs[[dgp_model$n_layers]]
  logits  <- final$mean

  if (ncol(logits) != dgp_model$n_celltypes) {
    if (ncol(logits) > dgp_model$n_celltypes) {
      logits <- logits[, seq_len(dgp_model$n_celltypes)]
    } else {
      extra  <- dgp_model$n_celltypes - ncol(logits)
      logits <- cbind(logits, matrix(stats::rnorm(nrow(logits) * extra, sd = 0.1),
                                     nrow = nrow(logits)))
    }
  }

  priors <- t(apply(logits, 1, function(x) {
    ex <- exp(x - max(x)); ex / sum(ex)
  }))
  priors <- pmax(priors, 1e-6)
  priors <- priors / rowSums(priors)

  list(priors = priors, uncertainty = sqrt(rowMeans(final$variance)))
}

# =========================================================================
# Stage 1: DGP-Enhanced MCMC Deconvolution
# =========================================================================

#' Run DGP-Enhanced Spatial Deconvolution (Stage 1)
#'
#' Performs Bayesian deconvolution for each spot using a negative binomial
#' likelihood with Dirichlet prior, enhanced by DGP-predicted spatial priors.
#' The sampler uses an ALR-parameterised Metropolis-Hastings algorithm with
#' adaptive proposal variance, thinning, and optional early stopping based on
#' split-chain R-hat and bulk ESS.
#'
#' @param spatial_count Count matrix (genes x spots).
#' @param ref_matrix Reference list from \code{\link{buildReferenceMatrix}}.
#' @param spatial_coords Spot coordinates.
#' @param alpha Dirichlet alpha vector from \code{\link{calculateAlpha}}.
#' @param kernel_mat Pre-computed spatial kernel matrix.
#' @param mcmc_iterations Total MCMC iterations per spot. Default 15000.
#' @param target_acceptance_range Target MH acceptance rate window.
#'   Default \code{c(0.23, 0.35)}.
#' @param spots_per_worker Approximate spots per parallel batch. Default 20.
#' @param early_stop Logical; enable early stopping. Default \code{TRUE}.
#' @param rhat_threshold R-hat convergence threshold. Default 1.05.
#' @param min_post_samples Minimum post-burn-in stored samples before
#'   convergence check. Default 100.
#' @param convergence_check_interval Iterations between convergence checks.
#'   Default 200.
#'
#' @return A list containing theta estimates, uncertainties, credible
#'   intervals, acceptance rates, and convergence diagnostics.
#' @export
runDeconvolution <- function(spatial_count, ref_matrix, spatial_coords,
                             alpha, kernel_mat,
                             mcmc_iterations = 15000,
                             target_acceptance_range = c(0.23, 0.35),
                             spots_per_worker = 20,
                             early_stop = TRUE,
                             rhat_threshold = 1.05,
                             min_post_samples = 100,
                             convergence_check_interval = 200) {

  # Built-in parameters
  initial_theta_sd <- 0.5
  initial_phi_sd   <- 0.3
  mcmc_burnin      <- max(500, floor(mcmc_iterations * 0.4))
  adapt_period     <- mcmc_burnin
  target_accept    <- mean(target_acceptance_range)
  adaptation_rate  <- 0.05
  min_prop_sd      <- 0.01
  max_prop_sd      <- 3.0
  mcmc_thin        <- 3
  dgp_weight       <- 0.4
  spatial_weight   <- 0.3
  n_dgp_layers     <- 3
  dgp_n_inducing   <- min(50, max(20, floor(nrow(spatial_coords) * 0.1)))

  n_spots   <- ncol(spatial_count)
  ref_basis <- ref_matrix$basis
  n_ct      <- ncol(ref_basis)

  message("=== DGP-Enhanced Spatial Deconvolution ===")
  message("   ", nrow(spatial_count), " genes x ", n_spots, " spots, ",
          n_ct, " cell types")
  message("   MCMC: ", mcmc_iterations, " iter (burn-in ", mcmc_burnin,
          ", thin ", mcmc_thin, ")")

  # --- DGP model ---
  message("1. Building DGP model...")
  dgp_model <- buildSpatialDGP(spatial_coords, ref_matrix, kernel_mat,
                                n_dgp_layers, dgp_n_inducing)
  message("2. Predicting spatial priors...")
  dgp_pred    <- predictCelltypePriors(dgp_model, spatial_coords, spatial_weight)
  dgp_priors  <- dgp_pred$priors
  uncertainty <- dgp_pred$uncertainty

  # --- Adaptive proposal helper ---
  adapt_sd <- function(current, acc, tgt_range, rate) {
    if (acc < tgt_range[1])      new_sd <- current * (1 - rate)
    else if (acc > tgt_range[2]) new_sd <- current * (1 + rate)
    else                         new_sd <- current
    pmax(pmin(new_sd, max_prop_sd), min_prop_sd)
  }

  # --- Parallel setup ---
  message("3. Setting up parallel backend...")
  ncores       <- parallel::detectCores() - 1
  max_by_data  <- max(1, floor(n_spots / 30))
  actual_cores <- min(ncores, max_by_data, 12)
  message("   Workers: ", actual_cores)

  n_workers <- min(ceiling(n_spots / spots_per_worker), 2000)
  spot_list <- split(seq_len(n_spots),
                     cut(seq_len(n_spots), n_workers, labels = FALSE))

  uncertainty_q90 <- stats::quantile(uncertainty, 0.9, na.rm = TRUE)

  # --- Build a shared parameter list so workers receive everything explicitly ---
  worker_params <- list(
    ref_basis      = ref_basis,
    alpha          = alpha,
    n_ct           = n_ct,
    uncertainty_q90 = uncertainty_q90,
    mcmc_thin      = mcmc_thin,
    early_stop     = early_stop,
    rhat_threshold = rhat_threshold,
    min_post_samples = min_post_samples,
    convergence_check_interval = convergence_check_interval,
    initial_theta_sd = initial_theta_sd,
    initial_phi_sd   = initial_phi_sd,
    mcmc_iterations  = mcmc_iterations,
    mcmc_burnin      = mcmc_burnin,
    adapt_period     = adapt_period,
    target_accept    = target_accept,
    target_acceptance_range = target_acceptance_range,
    adaptation_rate  = adaptation_rate,
    min_prop_sd      = min_prop_sd,
    max_prop_sd      = max_prop_sd,
    dgp_weight       = dgp_weight,
    spatial_count    = spatial_count,
    dgp_priors       = dgp_priors,
    uncertainty      = uncertainty
  )

  if (.Platform$OS.type == "unix") {
    cl <- parallel::makeCluster(actual_cores, type = "FORK")
    message("   Cluster type: FORK (copy-on-write)")
  } else {
    cl <- parallel::makeCluster(actual_cores)
    message("   Cluster type: PSOCK (Windows)")

    # On Windows (PSOCK), each worker is a blank R session.
    # We must load the SDDGP package so that the compiled C++ DLL
    # (mcmc_deconv_cpp etc.) is available in every worker.
    load_ok <- tryCatch({
      parallel::clusterEvalQ(cl, {
        library(SDDGP)
        # Verify the C++ function is callable
        exists("mcmc_deconv_cpp", mode = "function")
      })
    }, error = function(e) {
      # If SDDGP is not installed (e.g. loaded via devtools::load_all()),
      # fall back to finding and loading the DLL manually
      message("   Note: SDDGP not installed; loading DLL directly in workers...")
      dll_path <- getLoadedDLLs()[["SDDGP"]][["path"]]
      parallel::clusterExport(cl, "dll_path", envir = environment())
      parallel::clusterEvalQ(cl, dyn.load(dll_path))
      NULL
    })

    # Export the R wrapper functions that workers need
    parallel::clusterExport(cl, c(".emptySpotResult"),
                            envir = asNamespace("SDDGP"))
  }
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # Define the worker function in the main session so it captures nothing
  # implicitly — everything comes through the parameter list P.
  .workerFn <- function(spot_idx, P) {
    ref_basis      <- P$ref_basis
    alpha          <- P$alpha
    n_ct           <- P$n_ct
    uncertainty_q90 <- P$uncertainty_q90
    mcmc_thin      <- P$mcmc_thin
    early_stop     <- P$early_stop
    rhat_threshold <- P$rhat_threshold
    min_post_samples <- P$min_post_samples
    convergence_check_interval <- P$convergence_check_interval
    initial_theta_sd <- P$initial_theta_sd
    initial_phi_sd   <- P$initial_phi_sd
    mcmc_iterations  <- P$mcmc_iterations
    mcmc_burnin      <- P$mcmc_burnin
    adapt_period     <- P$adapt_period
    target_accept    <- P$target_accept
    dgp_weight       <- P$dgp_weight
    spatial_count    <- P$spatial_count
    dgp_priors       <- P$dgp_priors
    uncertainty      <- P$uncertainty

    .emptyResult <- function(nk, a, dp, u, tsd, psd, st = "skipped") {
      list(theta = rep(NA, nk), theta_sd = rep(NA, nk),
           theta_ci = matrix(NA, nk, 2), theta_cv = rep(NA, nk),
           phi = NA, phi_sd = NA, dgp_prior = dp, uncertainty = u,
           adaptive_weight = NA, enhanced_alpha = rep(NA, length(a)),
           acceptance_theta = NA, acceptance_phi = NA,
           proposal_theta_sd = tsd, proposal_phi_sd = psd,
           rhat = NA, min_ess = NA, n_stored = NA, actual_iterations = NA,
           converged = FALSE, status = st)
    }

    batch_theta_sd <- initial_theta_sd
    batch_phi_sd   <- initial_phi_sd

    lapply(spot_idx, function(i) {
      y    <- spatial_count[, i]
      dp   <- dgp_priors[i, ]
      u_sc <- uncertainty[i]

      if (sum(y) == 0 || any(!is.finite(dp)))
        return(.emptyResult(n_ct, alpha, dp, u_sc,
                            batch_theta_sd, batch_phi_sd))

      norm_u <- pmin(pmax(u_sc / uncertainty_q90, 0), 1)
      adap_w <- dgp_weight * (1 - norm_u)
      enh_a  <- alpha * (1 - adap_w) + (alpha * dp) * adap_w
      enh_a  <- pmax(enh_a, 1e-4)

      if (any(!is.finite(enh_a)))
        return(.emptyResult(n_ct, enh_a, dp, u_sc,
                            batch_theta_sd, batch_phi_sd))

      tryCatch({
        fit <- mcmc_deconv_cpp(
          y = as.integer(y), ref = ref_basis, alpha = enh_a,
          a_phi = 0.5, b_phi = 50, n_iter = mcmc_iterations,
          burn_in = mcmc_burnin, theta_proposal_sd = batch_theta_sd,
          phi_proposal_sd = batch_phi_sd, adapt_period = adapt_period,
          target_accept = target_accept, thin = mcmc_thin,
          early_stop = early_stop, rhat_threshold = rhat_threshold,
          min_post_samples = min_post_samples,
          convergence_check_interval = convergence_check_interval)

        th_mean <- colMeans(fit$theta_samples)
        th_sd   <- apply(fit$theta_samples, 2, stats::sd)
        th_ci   <- apply(fit$theta_samples, 2,
                         function(x) stats::quantile(x, c(0.025, 0.975)))
        th_cv   <- ifelse(th_mean > 1e-6, th_sd / th_mean, NA_real_)

        res <- list(
          theta = th_mean, theta_sd = th_sd, theta_ci = t(th_ci),
          theta_cv = th_cv,
          phi = mean(fit$phi_samples), phi_sd = stats::sd(fit$phi_samples),
          dgp_prior = dp, uncertainty = u_sc, adaptive_weight = adap_w,
          enhanced_alpha = enh_a,
          acceptance_theta = fit$acceptance_theta,
          acceptance_phi   = fit$acceptance_phi,
          proposal_theta_sd = batch_theta_sd,
          proposal_phi_sd   = batch_phi_sd,
          rhat = fit$rhat, min_ess = fit$min_ess,
          n_stored = fit$n_stored, actual_iterations = fit$actual_iterations,
          converged = fit$converged, status = "success")

        rm(fit); gc(verbose = FALSE)
        res
      }, error = function(e) {
        .emptyResult(n_ct, enh_a, dp, u_sc, batch_theta_sd, batch_phi_sd,
                     paste("error:", e$message))
      })
    })
  }

  message("4. Running MCMC sampling...")
  t0 <- Sys.time()

  results <- parallel::parLapplyLB(cl, spot_list, .workerFn,
                                   P = worker_params)

  t1 <- Sys.time()
  message("   MCMC complete: ",
          round(difftime(t1, t0, units = "mins"), 1), " minutes")

  # --- Consolidate ---
  message("5. Consolidating results...")
  all_res <- unlist(results, recursive = FALSE)

  .consolidateStage1(all_res, ref_basis, spatial_count, n_spots,
                     mcmc_iterations, early_stop, target_acceptance_range)
}

# ---------- internal helpers ----------

#' @keywords internal
.emptySpotResult <- function(n_ct, alpha, dgp_prior, unc,
                             th_sd, ph_sd, status = "skipped") {
  list(theta = rep(NA, n_ct), theta_sd = rep(NA, n_ct),
       theta_ci = matrix(NA, n_ct, 2), theta_cv = rep(NA, n_ct),
       phi = NA, phi_sd = NA, dgp_prior = dgp_prior, uncertainty = unc,
       adaptive_weight = NA, enhanced_alpha = rep(NA, length(alpha)),
       acceptance_theta = NA, acceptance_phi = NA,
       proposal_theta_sd = th_sd, proposal_phi_sd = ph_sd,
       rhat = NA, min_ess = NA, n_stored = NA, actual_iterations = NA,
       converged = FALSE, status = status)
}

#' @keywords internal
.consolidateStage1 <- function(all_res, ref_basis, spatial_count,
                               n_spots, mcmc_iterations, early_stop,
                               target_range) {
  ct_names  <- colnames(ref_basis)
  sp_names  <- colnames(spatial_count)

  theta_est  <- do.call(rbind, lapply(all_res, `[[`, "theta"))
  theta_sd   <- do.call(rbind, lapply(all_res, `[[`, "theta_sd"))
  theta_ci_l <- do.call(rbind, lapply(all_res, function(x) x$theta_ci[, 1]))
  theta_ci_u <- do.call(rbind, lapply(all_res, function(x) x$theta_ci[, 2]))
  theta_cv   <- do.call(rbind, lapply(all_res, `[[`, "theta_cv"))

  phi_est <- vapply(all_res, `[[`, numeric(1), "phi")
  phi_sd  <- vapply(all_res, `[[`, numeric(1), "phi_sd")

  acc_th <- vapply(all_res, `[[`, numeric(1), "acceptance_theta")
  acc_ph <- vapply(all_res, `[[`, numeric(1), "acceptance_phi")

  sp_rhat     <- vapply(all_res, `[[`, numeric(1), "rhat")
  sp_ess      <- vapply(all_res, `[[`, numeric(1), "min_ess")
  sp_stored   <- vapply(all_res, `[[`, numeric(1), "n_stored")
  sp_iter     <- vapply(all_res, `[[`, numeric(1), "actual_iterations")
  sp_conv     <- vapply(all_res, `[[`, logical(1), "converged")

  dimnames(theta_est)  <- list(sp_names, ct_names)
  dimnames(theta_sd)   <- list(sp_names, ct_names)
  dimnames(theta_ci_l) <- list(sp_names, ct_names)
  dimnames(theta_ci_u) <- list(sp_names, ct_names)
  dimnames(theta_cv)   <- list(sp_names, ct_names)
  names(phi_est) <- names(phi_sd) <- sp_names

  ci_width <- theta_ci_u - theta_ci_l
  spot_ui  <- apply(theta_cv, 1, function(r) {
    v <- r[!is.na(r) & is.finite(r)]
    if (length(v) == 0) NA_real_ else mean(v)
  })
  names(spot_ui) <- sp_names

  list(
    theta_estimates      = theta_est,
    theta_uncertainties  = theta_sd,
    theta_ci_lower       = theta_ci_l,
    theta_ci_upper       = theta_ci_u,
    theta_cv             = theta_cv,
    ci_width             = ci_width,
    spot_uncertainty_index = spot_ui,
    phi_estimates        = phi_est,
    phi_uncertainties    = phi_sd,
    acceptance_rates_theta = acc_th,
    acceptance_rates_phi   = acc_ph,
    mean_acceptance_theta  = mean(acc_th, na.rm = TRUE),
    mean_acceptance_phi    = mean(acc_ph, na.rm = TRUE),
    convergence = list(
      rhat              = sp_rhat,
      min_ess           = sp_ess,
      n_stored          = sp_stored,
      actual_iterations = sp_iter,
      converged         = sp_conv,
      pct_converged     = round(mean(sp_conv, na.rm = TRUE) * 100, 1),
      mean_rhat         = mean(sp_rhat, na.rm = TRUE),
      max_rhat          = max(sp_rhat, na.rm = TRUE),
      mean_ess          = mean(sp_ess, na.rm = TRUE),
      min_ess_global    = min(sp_ess, na.rm = TRUE),
      mean_actual_iter  = mean(sp_iter, na.rm = TRUE),
      total_saved_iterations = sum(mcmc_iterations - sp_iter, na.rm = TRUE)
    )
  )
}

# =========================================================================
# Resolution Enhancement
# =========================================================================

#' CARD-Style Resolution Enhancement
#'
#' Interpolates cell-type proportions from observed spot locations onto a
#' dense grid within the tissue boundary, using Matérn-weighted k-NN
#' interpolation implemented in C++.
#'
#' @param deconv_results Stage 1 results from \code{\link{runDeconvolution}}.
#' @param ref_matrix Reference list from \code{\link{buildReferenceMatrix}}.
#' @param spatial_coords Spot coordinates with columns \code{x} and \code{y}.
#' @param num_grids Target number of grid points. Default 2000.
#' @param n_neighbors Number of neighbours for interpolation. Default 10.
#' @param bandwidth Kernel bandwidth; if \code{NULL} auto-computed.
#' @param concavity Concavity parameter for tissue boundary detection.
#'   Default 2.0.
#'
#' @return A list with enhanced coordinates, proportions, imputed expression,
#'   and interpolation confidence.
#' @export
enhanceResolution <- function(deconv_results, ref_matrix, spatial_coords,
                              num_grids = 2000, n_neighbors = 10,
                              bandwidth = NULL, concavity = 2.0) {
  message("=== Resolution Enhancement ===")
  theta_obs  <- deconv_results$theta_estimates
  coords_obs <- as.matrix(spatial_coords[, c("x", "y")])
  n_obs <- nrow(coords_obs)

  # Tissue boundary
  pts_sf  <- sf::st_as_sf(as.data.frame(coords_obs), coords = c("x", "y"))
  hull    <- concaveman::concaveman(pts_sf, concavity = concavity)
  hull_pg <- sf::st_geometry(hull)[[1]]

  # Grid
  bbox    <- sf::st_bbox(hull)
  x_r     <- bbox["xmax"] - bbox["xmin"]
  y_r     <- bbox["ymax"] - bbox["ymin"]
  a_ratio <- as.numeric(sf::st_area(hull)) / (x_r * y_r)
  n_adj   <- ceiling(num_grids / a_ratio)
  nx <- ceiling(sqrt(n_adj * x_r / y_r))
  ny <- ceiling(sqrt(n_adj * y_r / x_r))
  grid_full <- expand.grid(x = seq(bbox["xmin"], bbox["xmax"], length.out = nx),
                           y = seq(bbox["ymin"], bbox["ymax"], length.out = ny))
  inside    <- sf::st_intersects(
    sf::st_as_sf(grid_full, coords = c("x", "y")), hull, sparse = FALSE)[, 1]
  grid_new  <- grid_full[inside, ]
  message("   Grid points inside tissue: ", nrow(grid_new))

  if (is.null(bandwidth)) {
    nn_d      <- RANN::nn2(coords_obs, coords_obs, k = 2)$nn.dists[, 2]
    bandwidth <- stats::median(nn_d) * 1.5
  }

  # Interpolate (C++)
  theta_new <- spatial_interpolate_theta_cpp(
    theta_obs, coords_obs, as.matrix(grid_new), bandwidth, n_neighbors, 1.5)
  colnames(theta_new) <- colnames(theta_obs)

  # Impute expression
  expr_new <- ref_matrix$basis %*% t(theta_new) * ref_matrix$median_lib_size

  nn_info   <- RANN::nn2(coords_obs, as.matrix(grid_new), k = 1)
  conf      <- exp(-nn_info$nn.dists[, 1] / bandwidth)

  message("   ", n_obs, " -> ", nrow(grid_new),
          " (~", round(nrow(grid_new) / n_obs, 1), "x)")

  list(refined_theta = theta_new, refined_expression = expr_new,
       refined_coords = grid_new, refined_confidence = conf,
       original_theta = theta_obs,
       original_coords = as.data.frame(coords_obs),
       n_original = n_obs, n_enhanced = nrow(grid_new),
       bandwidth = bandwidth, cell_types = colnames(theta_obs))
}

# =========================================================================
# Stage 2: Single-Cell Refinement
# =========================================================================

#' Stage 2 Single-Cell Refinement
#'
#' Refines spot-level proportions to single-cell resolution by matching each
#' spot to a set of candidate pseudo-cells via cosine similarity and solving
#' a cell-type-budget-constrained NNLS problem per spot.
#'
#' @param stage1_results Results from \code{\link{runDeconvolution}}.
#' @param ref_matrix Reference list from \code{\link{buildReferenceMatrix}}.
#' @param spatial_count Full spatial count matrix (all genes, not just
#'   informative genes).
#' @param n_candidates Number of candidate pseudo-cells per spot. Default 50.
#' @param lambda_reg Regularisation parameter for NNLS. Default 0.01.
#' @param budget_slack Slack factor for cell-type budget constraints.
#'   Default 1.2.
#'
#' @return A list with cell weights, assignments, reconstruction RMSE, etc.
#' @export
runStage2Refinement <- function(stage1_results, ref_matrix, spatial_count,
                                n_candidates = 50, lambda_reg = 0.01,
                                budget_slack = 1.2) {
  message("=== Stage 2: Single-Cell Refinement ===")
  theta    <- stage1_results$theta_estimates
  n_spots  <- nrow(theta); K <- ncol(theta)
  ct_names <- colnames(theta)

  sc_ref    <- ref_matrix$sc_sampled
  type_labs <- ref_matrix$pseudo_type_labels
  n_pseudo  <- ncol(sc_ref)

  cg <- intersect(rownames(spatial_count), rownames(sc_ref))
  message("   Shared genes: ", length(cg))

  sp_sub <- as.matrix(spatial_count[cg, ])
  sc_sub <- as.matrix(sc_ref[cg, ])

  sp_tpm <- sweep(sp_sub + 0.5, 2, colSums(sp_sub + 0.5), "/") * 1e6
  sp_log <- log2(sp_tpm + 1)
  sc_tpm <- sweep(sc_sub + 0.5, 2, colSums(sc_sub + 0.5), "/") * 1e6
  sc_log <- log2(sc_tpm + 1)

  type_map    <- stats::setNames(0:(K - 1), ct_names)
  type_assign <- as.integer(type_map[type_labs])
  type_assign[is.na(type_assign)] <- -1L

  message("   Computing cosine similarity...")
  cos_sim <- cosine_similarity_matrix_cpp(sp_log, sc_log)

  message("   Constrained NNLS per spot...")
  cell_weights <- matrix(0, n_spots, n_pseudo, dimnames = list(rownames(theta), NULL))

  for (s in seq_len(n_spots)) {
    top_idx <- order(cos_sim[s, ], decreasing = TRUE)[seq_len(min(n_candidates, n_pseudo))]
    budget  <- pmax(as.numeric(theta[s, ]) * budget_slack, 1e-6)
    w <- constrained_nnls_spot_cpp(sp_log[, s], sc_log[, top_idx, drop = FALSE],
                                   budget, type_assign[top_idx], lambda_reg, 200L)
    cell_weights[s, top_idx] <- as.numeric(w)
    if (s %% 200 == 0) message("     ", s, " / ", n_spots)
  }

  # Assignments
  cell_assigns <- lapply(seq_len(n_spots), function(s) {
    w  <- cell_weights[s, ]
    nz <- which(w > 1e-6)
    if (length(nz) == 0) return(data.frame(cell_idx = integer(0),
                                           weight = numeric(0),
                                           cell_type = character(0),
                                           stringsAsFactors = FALSE))
    nz <- nz[order(w[nz], decreasing = TRUE)]
    data.frame(cell_idx = nz, weight = w[nz], cell_type = type_labs[nz],
               stringsAsFactors = FALSE)
  })
  names(cell_assigns) <- rownames(theta)

  # Reconstruction error
  pred  <- sc_log %*% t(cell_weights)
  rmse  <- vapply(seq_len(n_spots),
                  function(s) sqrt(mean((pred[, s] - sp_log[, s])^2)),
                  numeric(1))

  message("   Mean RMSE: ", round(mean(rmse), 4))

  list(cell_weights = cell_weights, cell_assignments = cell_assigns,
       cell_type_labels = type_labs, reconstruction_rmse = rmse,
       n_pseudo_cells = n_pseudo)
}
