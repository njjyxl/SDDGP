#' Build Reference Expression Matrix from Single-Cell Data
#'
#' Constructs a cell type-specific reference expression matrix from
#' single-cell RNA-seq data, performs quality control, selects informative
#' marker genes, and samples pseudo-cells for Stage 2 refinement.
#'
#' @param sc_count Gene-by-cell count matrix (dense or sparse).
#' @param spatial_count Gene-by-spot count matrix (dense or sparse).
#' @param cell_type A data frame containing cell type annotations.
#' @param ct_varname Character string specifying the column name in
#'   \code{cell_type} that holds cell type labels. Default \code{"celltype"}.
#' @param min_count_gene Minimum number of cells a gene must be detected in
#'   to pass quality control. Default 100.
#' @param min_count_spot Minimum number of spots a gene must be detected in.
#'   Default 3.
#' @param n_pseudo_per_type Number of pseudo-cells to sample per cell type
#'   for Stage 2. Default 200.
#' @param pseudo_seed Random seed for pseudo-cell sampling reproducibility.
#'   Default 42.
#'
#' @return A list with components:
#'   \describe{
#'     \item{basis}{Normalized reference expression matrix
#'       (informative genes x cell types).}
#'     \item{selected_genes}{Character vector of selected informative genes.}
#'     \item{cell_types}{Character vector of retained cell types.}
#'     \item{median_lib_size}{Median library size of the spatial data
#'       (informative genes only).}
#'     \item{sc_sampled}{Count matrix of sampled pseudo-cells
#'       (all filtered genes x pseudo-cells).}
#'     \item{pseudo_type_labels}{Cell type label for each pseudo-cell.}
#'   }
#'
#' @export
buildReferenceMatrix <- function(sc_count, spatial_count, cell_type,
                                 ct_varname = "celltype",
                                 min_count_gene = 100,
                                 min_count_spot = 3,
                                 n_pseudo_per_type = 200,
                                 pseudo_seed = 42) {

  message("=== Reference Expression Matrix Construction ===")

  # --- Step 1: Ensure sparse format & intersect genes ---
  message("1. Data preprocessing...")
  if (!methods::is(sc_count, "dgCMatrix"))
    sc_count <- methods::as(sc_count, "dgCMatrix")
  if (!methods::is(spatial_count, "dgCMatrix"))
    spatial_count <- methods::as(spatial_count, "dgCMatrix")

  common_genes <- intersect(rownames(sc_count), rownames(spatial_count))
  message("   Common genes: ", length(common_genes))

  # Remove mitochondrial and ribosomal genes
  common_genes <- common_genes[!grepl("^mt-|^MT-", common_genes)]
  common_genes <- common_genes[!grepl("^rp[sl]|^RP[SL]", common_genes)]
  message("   After filtering mt/rp: ", length(common_genes))

  sc_count      <- sc_count[common_genes, ]
  spatial_count <- spatial_count[common_genes, ]

  # --- Step 2: Quality control ---
  message("2. Quality control...")
  gene_det   <- Matrix::rowSums(sc_count > 0)
  valid_genes <- names(gene_det)[gene_det >= min_count_gene]

  sp_gene_det <- Matrix::rowSums(spatial_count > 0)
  valid_genes <- intersect(valid_genes,
                           names(sp_gene_det)[sp_gene_det >= min_count_spot])
  message("   Genes after QC: ", length(valid_genes))

  sc_count      <- sc_count[valid_genes, ]
  spatial_count <- spatial_count[valid_genes, ]

  # --- Step 3: Cell type expression profiles ---
  message("3. Building cell type expression profiles...")
  cell_types <- unique(cell_type[[ct_varname]])
  cell_types <- cell_types[!is.na(cell_types)]

  basis_matrix <- matrix(0, nrow = length(valid_genes), ncol = length(cell_types),
                         dimnames = list(valid_genes, cell_types))

  ct_cell_indices <- list()

  for (ct in cell_types) {
    idx <- which(cell_type[[ct_varname]] == ct)
    ct_cell_indices[[ct]] <- idx
    if (length(idx) >= 3) {
      ct_counts <- sc_count[, idx, drop = FALSE]
      cs <- Matrix::colSums(ct_counts)
      ct_cpm <- sweep(as.matrix(ct_counts), 2, cs, "/") * 1e6
      basis_matrix[, ct] <- rowMeans(log2(ct_cpm + 1))
    } else {
      message("   Warning: ", ct, " has only ", length(idx), " cells, skipped")
    }
  }

  valid_ct <- colSums(basis_matrix) > 0
  basis_matrix    <- basis_matrix[, valid_ct, drop = FALSE]
  cell_types      <- cell_types[valid_ct]
  ct_cell_indices <- ct_cell_indices[cell_types]

  # --- Step 4: Informative gene selection ---
  message("4. Informative gene selection...")
  info_list <- vector("list", length(cell_types))
  for (i in seq_along(cell_types)) {
    ct <- cell_types[i]
    current <- basis_matrix[, ct]
    others  <- rowMeans(basis_matrix[, setdiff(cell_types, ct), drop = FALSE])
    fc <- (current + 0.1) / (others + 0.1)
    info_list[[i]] <- names(fc)[fc > 2 & current > 1]
  }
  info_genes <- sort(unique(unlist(info_list)))
  message("   Informative genes: ", length(info_genes))

  # --- Step 5: Final reference matrix ---
  message("5. Building final reference matrix...")
  final_basis <- basis_matrix[info_genes, , drop = FALSE]
  final_basis <- sweep(final_basis, 2, colSums(final_basis), "/")

  # --- Step 6: Sample pseudo-cells for Stage 2 ---
  message("6. Sampling pseudo-cells...")
  set.seed(pseudo_seed)
  sampled_idx   <- c()
  sampled_types <- c()

  for (ct in cell_types) {
    avail  <- ct_cell_indices[[ct]]
    n_samp <- min(n_pseudo_per_type, length(avail))
    if (length(avail) > 0) {
      sel <- sample(avail, n_samp, replace = (n_samp > length(avail)))
      sampled_idx   <- c(sampled_idx, sel)
      sampled_types <- c(sampled_types, rep(ct, n_samp))
    }
  }

  sc_sampled <- as.matrix(sc_count[, sampled_idx, drop = FALSE])
  colnames(sc_sampled) <- paste0(sampled_types, "_", seq_along(sampled_types))
  message("   Pseudo-cells sampled: ", ncol(sc_sampled))

  median_lib <- stats::median(Matrix::colSums(spatial_count[info_genes, ]))

  ref <- list(
    basis             = final_basis,
    selected_genes    = info_genes,
    cell_types        = cell_types,
    median_lib_size   = median_lib,
    sc_sampled        = sc_sampled,
    pseudo_type_labels = sampled_types
  )

  message("=== Reference matrix complete: ",
          nrow(final_basis), " genes x ", ncol(final_basis), " cell types ===")
  return(ref)
}
