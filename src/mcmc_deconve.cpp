#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// ============================================================================
// ALR transform helpers (unchanged logic)
// ============================================================================

// [[Rcpp::export]]
Rcpp::NumericVector inv_alr_cpp(const Rcpp::NumericVector& alr_coords) {
  int n = alr_coords.size();
  double max_alc = *std::max_element(alr_coords.begin(), alr_coords.end());
  
  Rcpp::NumericVector result(n + 1);
  double sum_exp = 0.0;
  
  for (int i = 0; i < n; i++) {
    double scaled = alr_coords[i] - max_alc;
    double exp_scaled = std::exp(scaled);
    result[i] = exp_scaled;
    sum_exp += exp_scaled;
  }
  
  double denominator = sum_exp + std::exp(-max_alc);
  for (int i = 0; i < n; i++) {
    result[i] /= denominator;
  }
  result[n] = 1.0 / (sum_exp * std::exp(max_alc) + 1.0);
  
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector alr_cpp(const Rcpp::NumericVector& theta) {
  int n = theta.size();
  Rcpp::NumericVector result(n - 1);
  double ref = theta[n - 1];
  
  for (int i = 0; i < n - 1; i++) {
    result[i] = std::log(theta[i] / ref);
  }
  
  return result;
}

// ============================================================================
// Precompute log factorial (unchanged)
// ============================================================================

// [[Rcpp::export]]
std::vector<double> precompute_log_factorial(int max_value) {
  std::vector<double> log_factorial(max_value + 1);
  log_factorial[0] = 0.0;
  for (int i = 1; i <= max_value; i++) {
    log_factorial[i] = log_factorial[i - 1] + std::log(i);
  }
  return log_factorial;
}

// ============================================================================
// Optimized NB log-likelihood: branch-free inner loop
// ============================================================================

// [[Rcpp::export]]
double log_lik_nb_cpp(const Rcpp::IntegerVector& y, 
                      const Rcpp::NumericVector& mu, 
                      double phi,
                      const std::vector<double>& log_factorial_cache) {
  int n = y.size();
  double sum_log_lik = 0.0;
  double log_phi = std::log(phi);
  double lgamma_phi = std::lgamma(phi);   // hoist out of loop
  
  for (int i = 0; i < n; i++) {
    double log_mu_phi = std::log(mu[i] + phi);
    if (y[i] > 0) {
      sum_log_lik += std::lgamma(y[i] + phi) - log_factorial_cache[y[i]] - 
        lgamma_phi + phi * log_phi + 
        y[i] * std::log(mu[i]) - (y[i] + phi) * log_mu_phi;
    } else {
      sum_log_lik += phi * log_phi - phi * log_mu_phi;
    }
  }
  
  return sum_log_lik;
}

// [[Rcpp::export]]
double dirichlet_prior_cpp(const Rcpp::NumericVector& theta, 
                           const Rcpp::NumericVector& alpha) {
  double sum = 0.0;
  double theta_sum = 0.0;
  
  for (int i = 0; i < theta.size(); i++) {
    if (theta[i] <= 0.0) return R_NegInf;
    theta_sum += theta[i];
    sum += (alpha[i] - 1.0) * std::log(theta[i]);
  }
  
  if (std::abs(theta_sum - 1.0) > 1e-6) return R_NegInf;
  
  return sum;
}

// ============================================================================
// MCMC with thinning support
// Key optimization: thin parameter reduces stored samples, allows
// more iterations with same memory footprint, or same iterations faster
// ============================================================================

// ============================================================================
// Split-chain Rhat: split stored samples into two halves, compute R-hat
// per dimension, return max across dimensions.
// ============================================================================
static double compute_split_rhat(const Rcpp::NumericMatrix& samples, int n_stored) {
  if (n_stored < 4) return 99.0;  // not enough samples
  int K = samples.ncol();
  int half = n_stored / 2;

  double max_rhat = 0.0;
  for (int k = 0; k < K; k++) {
    // Chain A = first half, Chain B = second half
    double mean_a = 0.0, mean_b = 0.0;
    for (int i = 0; i < half; i++) {
      mean_a += samples(i, k);
      mean_b += samples(half + i, k);
    }
    mean_a /= half;
    mean_b /= half;

    double var_a = 0.0, var_b = 0.0;
    for (int i = 0; i < half; i++) {
      double da = samples(i, k) - mean_a;
      double db = samples(half + i, k) - mean_b;
      var_a += da * da;
      var_b += db * db;
    }
    var_a /= (half - 1);
    var_b /= (half - 1);

    double W = (var_a + var_b) / 2.0;
    double grand_mean = (mean_a + mean_b) / 2.0;
    double B = half * ((mean_a - grand_mean) * (mean_a - grand_mean) +
                       (mean_b - grand_mean) * (mean_b - grand_mean));

    double var_hat = ((half - 1.0) / half) * W + (1.0 / half) * B;
    double rhat = (W > 1e-16) ? std::sqrt(var_hat / W) : 1.0;
    if (rhat > max_rhat) max_rhat = rhat;
  }
  return max_rhat;
}

// ============================================================================
// Bulk ESS: per-dimension effective sample size (batch means estimator)
// Returns the minimum ESS across dimensions.
// ============================================================================
static double compute_min_ess(const Rcpp::NumericMatrix& samples, int n_stored) {
  if (n_stored < 10) return 0.0;
  int K = samples.ncol();

  // Batch means ESS: use sqrt(n) batch size
  int batch_size = std::max(1, (int)std::sqrt((double)n_stored));
  int n_batches = n_stored / batch_size;
  if (n_batches < 2) return static_cast<double>(n_stored);

  double min_ess = 1e18;
  for (int k = 0; k < K; k++) {
    double grand_mean = 0.0;
    for (int i = 0; i < n_stored; i++) grand_mean += samples(i, k);
    grand_mean /= n_stored;

    // Batch means
    double var_bm = 0.0;
    for (int b = 0; b < n_batches; b++) {
      double bm = 0.0;
      for (int i = b * batch_size; i < (b + 1) * batch_size; i++) {
        bm += samples(i, k);
      }
      bm /= batch_size;
      var_bm += (bm - grand_mean) * (bm - grand_mean);
    }
    var_bm /= (n_batches - 1);

    // Marginal variance
    double var_total = 0.0;
    for (int i = 0; i < n_stored; i++) {
      double d = samples(i, k) - grand_mean;
      var_total += d * d;
    }
    var_total /= (n_stored - 1);

    double tau_hat = (var_bm > 1e-16) ? batch_size * var_bm / (var_total + 1e-16) : 1.0;
    double ess_k = n_stored / std::max(tau_hat, 1.0);
    if (ess_k < min_ess) min_ess = ess_k;
  }
  return min_ess;
}

// [[Rcpp::export]]
Rcpp::List mcmc_deconv_cpp(const Rcpp::IntegerVector& y,
                           const Rcpp::NumericMatrix& ref,
                           const Rcpp::NumericVector& alpha,
                           double a_phi = 1.0,
                           double b_phi = 1.0,
                           int n_iter = 5000,
                           int burn_in = 2000,
                           double theta_proposal_sd = 0.05,
                           double phi_proposal_sd = 0.1,
                           int adapt_period = 500,
                           double target_accept = 0.234,
                           int thin = 1,
                           bool early_stop = false,
                           double rhat_threshold = 1.05,
                           int min_post_samples = 100,
                           int convergence_check_interval = 200) {
  int K = ref.ncol();
  int G = ref.nrow();
  Rcpp::NumericVector theta_current(K);
  for (int i = 0; i < K; i++) {
    theta_current[i] = 1.0 / K;
  }
  
  Rcpp::NumericVector alr_current = alr_cpp(theta_current);
  double phi_current = 10.0;
  
  // Thinning: store only every `thin`-th sample after burn-in
  int post_burn = n_iter - burn_in;
  int n_stored_max = (post_burn + thin - 1) / thin;  // ceiling division
  Rcpp::NumericMatrix final_theta(n_stored_max, K);
  Rcpp::NumericVector final_phi(n_stored_max);
  int store_idx = 0;
  
  // Log-posterior trace (one value per stored sample for diagnostics)
  Rcpp::NumericVector log_post_trace(n_stored_max);
  
  int accept_theta = 0;
  int accept_phi = 0;
  int adapt_count_theta = 0;
  int adapt_count_phi = 0;
  
  int max_y = *std::max_element(y.begin(), y.end());
  std::vector<double> log_factorial_cache = precompute_log_factorial(max_y);
  
  // Precompute lambda = ref %*% theta
  Rcpp::NumericVector lambda_current(G);
  for (int g = 0; g < G; g++) {
    double s = 0.0;
    for (int k = 0; k < K; k++) {
      s += ref(g, k) * theta_current[k];
    }
    lambda_current[g] = s;
  }
  
  double log_prior_phi_const = -std::lgamma(a_phi) + a_phi * std::log(b_phi);
  
  double log_lik_current = log_lik_nb_cpp(y, lambda_current, phi_current, log_factorial_cache);
  double log_prior_theta_current = dirichlet_prior_cpp(theta_current, alpha);
  double log_prior_phi_current = (a_phi - 1.0) * std::log(phi_current) - b_phi * phi_current + log_prior_phi_const;
  double log_post_current = log_lik_current + log_prior_theta_current + log_prior_phi_current;
  
  Rcpp::RNGScope scope;
  
  // Pre-allocate proposal vectors outside loop to avoid repeated allocation
  Rcpp::NumericVector lambda_proposal(G);
  
  // Convergence tracking
  bool converged = false;
  int actual_iter = n_iter;
  double final_rhat = 99.0;
  double final_ess = 0.0;
  int next_conv_check = burn_in + std::max(convergence_check_interval, min_post_samples * thin);
  
  for (int i = 0; i < n_iter; i++) {
    // Adaptive step size
    if ((i + 1) % adapt_period == 0 && i < burn_in) {
      double accept_rate_theta = static_cast<double>(adapt_count_theta) / adapt_period;
      double accept_rate_phi = static_cast<double>(adapt_count_phi) / adapt_period;
      
      if (accept_rate_theta > target_accept) {
        theta_proposal_sd = std::min(theta_proposal_sd * 1.1, 0.5);
      } else {
        theta_proposal_sd = std::max(theta_proposal_sd * 0.9, 0.001);
      }
      
      if (accept_rate_phi > target_accept) {
        phi_proposal_sd = std::min(phi_proposal_sd * 1.1, 0.5);
      } else {
        phi_proposal_sd = std::max(phi_proposal_sd * 0.9, 0.001);
      }
      
      adapt_count_theta = 0;
      adapt_count_phi = 0;
    }
    
    // ===== Update theta =====
    int alr_size = alr_current.size();
    Rcpp::NumericVector alr_proposal(alr_size);
    for (int j = 0; j < alr_size; j++) {
      alr_proposal[j] = alr_current[j] + R::rnorm(0.0, theta_proposal_sd);
    }
    
    Rcpp::NumericVector theta_proposal = inv_alr_cpp(alr_proposal);
    
    // Compute proposed lambda and log_posterior
    bool lambda_valid = true;
    for (int g = 0; g < G; g++) {
      double s = 0.0;
      for (int k = 0; k < K; k++) {
        s += ref(g, k) * theta_proposal[k];
      }
      if (s <= 0.0) { lambda_valid = false; break; }
      lambda_proposal[g] = s;
    }
    
    if (lambda_valid) {
      double log_lik_proposal = log_lik_nb_cpp(y, lambda_proposal, phi_current, log_factorial_cache);
      double log_prior_theta_proposal = dirichlet_prior_cpp(theta_proposal, alpha);
      double log_post_proposal = log_lik_proposal + log_prior_theta_proposal + log_prior_phi_current;
      double log_ratio_theta = log_post_proposal - log_post_current;
      
      if (!std::isnan(log_ratio_theta) && std::log(R::runif(0.0, 1.0)) < log_ratio_theta) {
        theta_current = theta_proposal;
        alr_current = alr_proposal;
        // Copy lambda_proposal -> lambda_current
        for (int g = 0; g < G; g++) lambda_current[g] = lambda_proposal[g];
        log_lik_current = log_lik_proposal;
        log_prior_theta_current = log_prior_theta_proposal;
        log_post_current = log_post_proposal;
        accept_theta++;
        if (i < burn_in) adapt_count_theta++;
      }
    }
    
    // ===== Update phi =====
    double phi_proposal = std::exp(std::log(phi_current) + R::rnorm(0.0, phi_proposal_sd));
    
    if (phi_proposal > 0.0) {
      double log_lik_phi_proposal = log_lik_nb_cpp(y, lambda_current, phi_proposal, log_factorial_cache);
      double log_prior_phi_proposal = (a_phi - 1.0) * std::log(phi_proposal) - b_phi * phi_proposal + log_prior_phi_const;
      double log_post_phi_proposal = log_lik_phi_proposal + log_prior_theta_current + log_prior_phi_proposal;
      double log_ratio_phi = log_post_phi_proposal - log_post_current + 
        std::log(phi_proposal) - std::log(phi_current);
      
      if (!std::isnan(log_ratio_phi) && std::log(R::runif(0.0, 1.0)) < log_ratio_phi) {
        phi_current = phi_proposal;
        log_lik_current = log_lik_phi_proposal;
        log_prior_phi_current = log_prior_phi_proposal;
        log_post_current = log_post_phi_proposal;
        accept_phi++;
        if (i < burn_in) adapt_count_phi++;
      }
    }
    
    // Store with thinning
    if (i >= burn_in && ((i - burn_in) % thin == 0)) {
      if (store_idx < n_stored_max) {
        for (int j = 0; j < K; j++) {
          final_theta(store_idx, j) = theta_current[j];
        }
        final_phi[store_idx] = phi_current;
        log_post_trace[store_idx] = log_post_current;
        store_idx++;
      }
    }
    
    // ===== Early stopping convergence check =====
    if (early_stop && i >= next_conv_check && store_idx >= min_post_samples) {
      double rhat_now = compute_split_rhat(final_theta, store_idx);
      double ess_now = compute_min_ess(final_theta, store_idx);
      
      if (rhat_now < rhat_threshold && ess_now >= min_post_samples * 0.5) {
        converged = true;
        actual_iter = i + 1;
        final_rhat = rhat_now;
        final_ess = ess_now;
        break;
      }
      // Schedule next check
      next_conv_check = i + convergence_check_interval;
    }
    
    if ((i + 1) % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }
  }
  
  // Compute final diagnostics on whatever was stored
  if (!converged) {
    actual_iter = n_iter;
    final_rhat = compute_split_rhat(final_theta, store_idx);
    final_ess = compute_min_ess(final_theta, store_idx);
  }
  
  // Trim to actual stored samples
  Rcpp::NumericMatrix theta_out(store_idx, K);
  Rcpp::NumericVector phi_out(store_idx);
  Rcpp::NumericVector lp_out(store_idx);
  for (int i = 0; i < store_idx; i++) {
    for (int j = 0; j < K; j++) theta_out(i, j) = final_theta(i, j);
    phi_out[i] = final_phi[i];
    lp_out[i] = log_post_trace[i];
  }
  
  double acceptance_theta = static_cast<double>(accept_theta) / actual_iter;
  double acceptance_phi = static_cast<double>(accept_phi) / actual_iter;
  
  return Rcpp::List::create(
    Rcpp::Named("theta_samples") = theta_out,
    Rcpp::Named("phi_samples") = phi_out,
    Rcpp::Named("log_posterior_trace") = lp_out,
    Rcpp::Named("acceptance_theta") = acceptance_theta,
    Rcpp::Named("acceptance_phi") = acceptance_phi,
    Rcpp::Named("final_theta_sd") = theta_proposal_sd,
    Rcpp::Named("final_phi_sd") = phi_proposal_sd,
    Rcpp::Named("rhat") = final_rhat,
    Rcpp::Named("min_ess") = final_ess,
    Rcpp::Named("n_stored") = store_idx,
    Rcpp::Named("actual_iterations") = actual_iter,
    Rcpp::Named("converged") = converged
  );
}

// ============================================================================
// Spatial interpolation (unchanged – already well optimized)
// ============================================================================

// [[Rcpp::export]]
Rcpp::NumericMatrix spatial_interpolate_theta_cpp(
    const Rcpp::NumericMatrix& theta_observed,
    const Rcpp::NumericMatrix& coords_observed,
    const Rcpp::NumericMatrix& coords_new,
    double bandwidth,
    int n_neighbors,
    double nu = 1.5) {
  
  int n_obs = coords_observed.nrow();
  int n_new = coords_new.nrow();
  int K = theta_observed.ncol();
  int k = std::min(n_neighbors, n_obs);
  
  Rcpp::NumericMatrix theta_new(n_new, K);
  
  int matern_type = 1;
  if (std::abs(nu - 0.5) < 1e-6) matern_type = 0;
  else if (std::abs(nu - 2.5) < 1e-6) matern_type = 2;
  else if (nu > 10.0) matern_type = 3;
  
  double inv_bw = 1.0 / bandwidth;
  double sqrt3 = std::sqrt(3.0);
  double sqrt5 = std::sqrt(5.0);
  
  std::vector<std::pair<double, int>> dist_idx(n_obs);
  std::vector<double> weights(k);
  
  for (int i = 0; i < n_new; i++) {
    double xi = coords_new(i, 0);
    double yi = coords_new(i, 1);
    
    for (int j = 0; j < n_obs; j++) {
      double dx = xi - coords_observed(j, 0);
      double dy = yi - coords_observed(j, 1);
      dist_idx[j] = std::make_pair(dx * dx + dy * dy, j);
    }
    
    std::partial_sort(dist_idx.begin(), dist_idx.begin() + k, dist_idx.end());
    
    double sum_w = 0.0;
    for (int j = 0; j < k; j++) {
      double d = std::sqrt(dist_idx[j].first);
      double sd = d * inv_bw;
      double w;
      
      switch (matern_type) {
        case 0: w = std::exp(-sd); break;
        case 1: { double v = sqrt3 * sd; w = (1.0 + v) * std::exp(-v); break; }
        case 2: { double v = sqrt5 * sd; w = (1.0 + v + v * v / 3.0) * std::exp(-v); break; }
        default: w = std::exp(-0.5 * sd * sd); break;
      }
      
      weights[j] = w;
      sum_w += w;
    }
    
    if (sum_w < 1e-12) sum_w = 1e-12;
    double inv_sum_w = 1.0 / sum_w;
    
    for (int j = 0; j < k; j++) {
      int obs_idx = dist_idx[j].second;
      double w = weights[j] * inv_sum_w;
      for (int c = 0; c < K; c++) {
        theta_new(i, c) += w * theta_observed(obs_idx, c);
      }
    }
    
    double row_sum = 0.0;
    for (int c = 0; c < K; c++) {
      if (theta_new(i, c) < 0.0) theta_new(i, c) = 0.0;
      row_sum += theta_new(i, c);
    }
    if (row_sum > 0.0) {
      double inv_rs = 1.0 / row_sum;
      for (int c = 0; c < K; c++) theta_new(i, c) *= inv_rs;
    }
  }
  
  return theta_new;
}

// ============================================================================
// Cosine similarity: spots(cols of A) vs cells(cols of B)
// Optimized: precompute norms, column-major access pattern
// Input:  A = genes × spots,  B = genes × cells
// Output: spots × cells
// ============================================================================

// [[Rcpp::export]]
Rcpp::NumericMatrix cosine_similarity_matrix_cpp(
    const Rcpp::NumericMatrix& A,
    const Rcpp::NumericMatrix& B) {
  
  int G = A.nrow();
  int n = A.ncol();  // spots
  int m = B.ncol();  // cells
  
  // Precompute inverse column norms
  std::vector<double> inv_norm_A(n);
  std::vector<double> inv_norm_B(m);
  
  for (int j = 0; j < n; j++) {
    double s = 0.0;
    for (int g = 0; g < G; g++) s += A(g, j) * A(g, j);
    inv_norm_A[j] = 1.0 / std::sqrt(s + 1e-24);
  }
  
  for (int j = 0; j < m; j++) {
    double s = 0.0;
    for (int g = 0; g < G; g++) s += B(g, j) * B(g, j);
    inv_norm_B[j] = 1.0 / std::sqrt(s + 1e-24);
  }
  
  Rcpp::NumericMatrix result(n, m);
  
  // Column-major friendly: iterate over cells in outer loop
  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      double dot = 0.0;
      for (int g = 0; g < G; g++) {
        dot += A(g, i) * B(g, j);
      }
      result(i, j) = dot * inv_norm_A[i] * inv_norm_B[j];
    }
  }
  
  return result;
}

// ============================================================================
// Constrained NNLS for a single spot
// Projected gradient descent with cell-type budget constraints
// ============================================================================

// [[Rcpp::export]]
Rcpp::NumericVector constrained_nnls_spot_cpp(
    const Rcpp::NumericVector& y_spot,
    const Rcpp::NumericMatrix& X_cells,
    const Rcpp::NumericVector& budget,
    const Rcpp::IntegerVector& type_assign,
    double lambda_reg = 0.01,
    int max_iter = 200) {
  
  int G = X_cells.nrow();
  int C = X_cells.ncol();
  int K = budget.size();
  
  Rcpp::NumericVector w(C, 0.0);
  
  // Precompute X^T y
  std::vector<double> XtY(C);
  for (int j = 0; j < C; j++) {
    double s = 0.0;
    for (int g = 0; g < G; g++) s += X_cells(g, j) * y_spot[g];
    XtY[j] = s;
  }
  
  // Precompute column squared norms for step size
  std::vector<double> col_sq_norms(C);
  double max_csn = 0.0;
  for (int j = 0; j < C; j++) {
    double s = 0.0;
    for (int g = 0; g < G; g++) s += X_cells(g, j) * X_cells(g, j);
    col_sq_norms[j] = s;
    if (s > max_csn) max_csn = s;
  }
  double step_size = 1.0 / (max_csn + lambda_reg + 1e-8);
  
  std::vector<double> residual(G);
  
  for (int iter = 0; iter < max_iter; iter++) {
    // residual = y - X * w
    for (int g = 0; g < G; g++) {
      double xw = 0.0;
      for (int j = 0; j < C; j++) {
        if (w[j] > 0.0) xw += X_cells(g, j) * w[j];  // skip zeros
      }
      residual[g] = y_spot[g] - xw;
    }
    
    // Gradient step
    for (int j = 0; j < C; j++) {
      double grad_pos = 0.0;
      for (int g = 0; g < G; g++) grad_pos += X_cells(g, j) * residual[g];
      w[j] += step_size * (grad_pos - lambda_reg * w[j]);
      if (w[j] < 0.0) w[j] = 0.0;  // project non-negative inline
    }
    
    // Project: cell-type budget constraints
    for (int k = 0; k < K; k++) {
      double type_sum = 0.0;
      for (int j = 0; j < C; j++) {
        if (type_assign[j] == k) type_sum += w[j];
      }
      if (type_sum > budget[k] && type_sum > 1e-12) {
        double scale = budget[k] / type_sum;
        for (int j = 0; j < C; j++) {
          if (type_assign[j] == k) w[j] *= scale;
        }
      }
    }
  }
  
  // Final normalization
  double total = 0.0;
  for (int j = 0; j < C; j++) total += w[j];
  if (total > 1e-12) {
    double inv_total = 1.0 / total;
    for (int j = 0; j < C; j++) w[j] *= inv_total;
  }
  
  return w;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//