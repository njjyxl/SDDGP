#' Print MCMC Convergence Report
#'
#' Displays a formatted summary of convergence diagnostics from Stage 1,
#' including R-hat, ESS, early stopping, and uncertainty statistics.
#'
#' @param stage1_results Results from \code{\link{runDeconvolution}}.
#'
#' @return Invisibly returns a quality label (\code{"GOOD"},
#'   \code{"CAUTION"}, or \code{"POOR"}).
#' @export
printConvergenceReport <- function(stage1_results) {
  conv <- stage1_results$convergence
  cat("\n")
  cat(strrep("=", 50), "\n")
  cat("  MCMC Convergence Diagnostic Report\n")
  cat(strrep("=", 50), "\n")

  cat("\n  Convergence (early stop)\n")
  cat("    Spots converged:    ", sprintf("%d / %d (%.1f%%)",
      sum(conv$converged, na.rm = TRUE), length(conv$converged),
      conv$pct_converged), "\n")
  cat("    Mean actual iters:  ", sprintf("%.0f", conv$mean_actual_iter), "\n")
  cat("    Total saved iters:  ", sprintf("%d", conv$total_saved_iterations), "\n")

  cat("\n  Split-chain R-hat (target < 1.05)\n")
  cat("    Mean:               ", sprintf("%.4f", conv$mean_rhat), "\n")
  cat("    Max:                ", sprintf("%.4f", conv$max_rhat), "\n")
  cat("    Spots R-hat > 1.1:  ", sum(conv$rhat > 1.1, na.rm = TRUE), "\n")

  cat("\n  Effective Sample Size\n")
  cat("    Mean min-ESS:       ", sprintf("%.1f", conv$mean_ess), "\n")
  cat("    Global min-ESS:     ", sprintf("%.1f", conv$min_ess_global), "\n")
  cat("    Spots ESS < 50:     ", sum(conv$min_ess < 50, na.rm = TRUE), "\n")

  cat("\n  Uncertainty\n")
  cat("    Mean 95%% CI width: ", sprintf("%.4f",
      mean(stage1_results$ci_width, na.rm = TRUE)), "\n")
  cat("    Mean spot CV index: ", sprintf("%.4f",
      mean(stage1_results$spot_uncertainty_index, na.rm = TRUE)), "\n")

  quality <- "GOOD"
  if (conv$max_rhat > 1.1 || conv$min_ess_global < 50) quality <- "CAUTION"
  if (conv$max_rhat > 1.2 || conv$min_ess_global < 20) quality <- "POOR"
  cat("\n  Overall quality: ", quality, "\n")
  cat(strrep("=", 50), "\n\n")
  invisible(quality)
}

#' Plot Convergence Diagnostics
#'
#' Produces a multi-panel base-R figure showing R-hat distribution, ESS
#' distribution, iteration savings, and (optionally) a spatial uncertainty
#' map.
#'
#' @param stage1_results Results from \code{\link{runDeconvolution}}.
#' @param spatial_coords Optional spatial coordinates for the uncertainty map.
#'
#' @export
plotConvergenceDiagnostics <- function(stage1_results, spatial_coords = NULL) {
  conv <- stage1_results$convergence

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  # Panel 1: R-hat
  rhat_vals <- conv$rhat[!is.na(conv$rhat)]
  graphics::hist(rhat_vals, breaks = 30, col = "#2ECC71",
       main = "Split-chain R-hat Distribution",
       xlab = "R-hat", ylab = "# Spots", border = "white")
  graphics::abline(v = 1.05, col = "red", lwd = 2, lty = 2)

  # Panel 2: ESS
  ess_vals <- conv$min_ess[!is.na(conv$min_ess)]
  graphics::hist(ess_vals, breaks = 30, col = "#3498DB",
       main = "Minimum ESS Distribution",
       xlab = "Min ESS per spot", ylab = "# Spots", border = "white")
  graphics::abline(v = 50, col = "red", lwd = 2, lty = 2)
  graphics::abline(v = 100, col = "orange", lwd = 2, lty = 2)

  # Panel 3: Iterations
  iter_vals <- conv$actual_iterations[!is.na(conv$actual_iterations)]
  max_iter  <- max(iter_vals)
  graphics::hist(iter_vals, breaks = 30, col = "#9B59B6",
       main = "Actual Iterations per Spot",
       xlab = "Iterations used", ylab = "# Spots", border = "white")
  savings <- round((1 - mean(iter_vals) / max_iter) * 100, 1)
  graphics::legend("topleft", legend = paste0("Mean saving: ", savings, "%"),
         bty = "n", cex = 0.9)

  # Panel 4: Spatial uncertainty
  if (!is.null(spatial_coords)) {
    ui <- stage1_results$spot_uncertainty_index
    ui[is.na(ui)] <- 0
    pal <- grDevices::colorRampPalette(c("#2ECC71", "#F1C40F", "#E74C3C"))(100)
    ui_s <- pmin(ui / max(ui + 1e-10), 1)
    cols <- pal[pmax(1, ceiling(ui_s * 100))]
    graphics::plot(spatial_coords[, 1], spatial_coords[, 2], col = cols,
         pch = 16, cex = 0.8, main = "Spatial Uncertainty Index",
         xlab = "x", ylab = "y", asp = 1)
    graphics::legend("topright", legend = c("Low", "Medium", "High"),
           col = c("#2ECC71", "#F1C40F", "#E74C3C"), pch = 16, bty = "n")
  }
}

#' Plot Cell-Type-Specific Proportion and Uncertainty Maps
#'
#' Generates ggplot2 panels showing estimated proportions, 95\% CI width,
#' and an estimate-vs-uncertainty scatter for a single cell type.
#'
#' @param stage1_results Results from \code{\link{runDeconvolution}}.
#' @param spatial_coords Spatial coordinates.
#' @param cell_type Character; cell type to plot. Defaults to the first one.
#'
#' @return Invisibly, a list of three ggplot objects.
#' @export
plotCelltypeUncertainty <- function(stage1_results, spatial_coords,
                                   cell_type = NULL) {
  if (is.null(cell_type))
    cell_type <- colnames(stage1_results$theta_estimates)[1]

  df <- data.frame(
    x        = spatial_coords[, 1],
    y        = spatial_coords[, 2],
    theta    = stage1_results$theta_estimates[, cell_type],
    ci_width = stage1_results$ci_width[, cell_type],
    sd       = stage1_results$theta_uncertainties[, cell_type]
  )

  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y,
                                          color = .data$theta)) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_color_viridis_c(option = "C") +
    ggplot2::labs(title = paste(cell_type, "- Proportion"),
                  color = "Proportion") +
    ggplot2::theme_minimal() + ggplot2::coord_fixed()

  p2 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y,
                                          color = .data$ci_width)) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_color_viridis_c(option = "A", direction = -1) +
    ggplot2::labs(title = paste(cell_type, "- 95% CI Width"),
                  color = "CI Width") +
    ggplot2::theme_minimal() + ggplot2::coord_fixed()

  p3 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$theta,
                                          y = .data$ci_width)) +
    ggplot2::geom_point(alpha = 0.4, size = 0.8) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "red") +
    ggplot2::labs(title = paste(cell_type, "- Estimate vs Uncertainty"),
                  x = "Proportion", y = "95% CI Width") +
    ggplot2::theme_minimal()

  print(p1); print(p2); print(p3)
  invisible(list(proportion_map = p1, uncertainty_map = p2, calibration = p3))
}

#' Plot Enhanced-Resolution Proportions
#'
#' @param enhanced_results Results from \code{\link{enhanceResolution}}.
#' @param cell_type Cell type to plot. Defaults to the first one.
#'
#' @return A ggplot object.
#' @export
plotEnhancedProportions <- function(enhanced_results, cell_type = NULL) {
  coords <- enhanced_results$refined_coords
  theta  <- enhanced_results$refined_theta
  if (is.null(cell_type)) cell_type <- colnames(theta)[1]

  df <- data.frame(x = coords$x, y = coords$y,
                   proportion = theta[, cell_type],
                   confidence = enhanced_results$refined_confidence)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y,
                                    color = .data$proportion)) +
    ggplot2::geom_point(size = 0.5, alpha = pmax(df$confidence, 0.3)) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(title = paste("Enhanced:", cell_type),
                  subtitle = paste(enhanced_results$n_enhanced, "from",
                                   enhanced_results$n_original, "spots")) +
    ggplot2::theme_minimal() + ggplot2::coord_fixed()
}

#' Plot Single-Cell Composition of a Spot
#'
#' @param stage2_results Results from \code{\link{runStage2Refinement}}.
#' @param spot_name Name of the spot to inspect.
#' @param top_n Maximum number of cells to display. Default 10.
#'
#' @return A ggplot object, or \code{NULL} invisibly if no cells assigned.
#' @export
plotSingleCellComposition <- function(stage2_results, spot_name, top_n = 10) {
  assignment <- stage2_results$cell_assignments[[spot_name]]
  if (is.null(assignment) || nrow(assignment) == 0) {
    message("No cells assigned to spot: ", spot_name)
    return(invisible(NULL))
  }
  n_show <- min(top_n, nrow(assignment))
  df <- assignment[seq_len(n_show), ]
  df$label <- paste0("Cell_", df$cell_idx, " (", df$cell_type, ")")
  df$label <- factor(df$label, levels = rev(df$label))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$weight, y = .data$label,
                                    fill = .data$cell_type)) +
    ggplot2::geom_col() +
    ggplot2::labs(title = paste("Decomposition:", spot_name),
                  x = "Weight", y = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right")
}
