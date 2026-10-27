TPDc_fast <- function(TPDs, sampUnit) {
  sampUnit <- as.matrix(sampUnit)
  
  results      <- list()
  results$data <- TPDs$data
  results$data$sampUnit <- sampUnit
  
  eval_grid_nrow <- nrow(results$data$evaluation_grid)
  TPDs_mat <- do.call(cbind, TPDs$TPDs)  # 预先合并为矩阵，避免反复索引list
  
  samp_names  <- rownames(sampUnit)
  TPDc_list   <- vector("list", nrow(sampUnit))
  Mcount_list <- vector("list", nrow(sampUnit))  # 替代RTPDs，只存M
  abund_list  <- vector("list", nrow(sampUnit))
  species_list <- vector("list", nrow(sampUnit))
  
  for (samp in seq_len(nrow(sampUnit))) {
    present_idx <- which(sampUnit[samp, ] > 0)
    n_spp <- length(present_idx)
    if (n_spp == 0) next
    
    spp_mat <- TPDs_mat[, present_idx, drop = FALSE]  # 直接矩阵切片
    
    # TPDc = 行均值（等权重）
    TPDc_list[[samp]]   <- rowMeans(spp_mat)
    
    # M = 每格点物种数（FRed用），用逻辑矩阵行求和
    Mcount_list[[samp]] <- rowSums(spp_mat > 0)
    
    present_spp <- colnames(sampUnit)[present_idx]
    abund_list[[samp]]   <- setNames(rep(1/n_spp, n_spp), present_spp)
    species_list[[samp]] <- present_spp
  }
  
  names(TPDc_list) <- names(Mcount_list) <- 
    names(abund_list) <- names(species_list) <- samp_names
  
  results$TPDc <- list(
    TPDc       = TPDc_list,
    Mcount     = Mcount_list,  # 替代RTPDs
    abundances = abund_list,
    species    = species_list
  )
  
  class(results) <- "TPDcomm"
  return(results)
}



fast_FRic <- function(species, means, sds, abundance_mat,
                            covar        = FALSE,
                            alpha        = 0.95,
                            trait_ranges = NULL,
                            n_divisions  = 70,
                            n_threads    = 2) {
  means         <- as.matrix(means)
  sds           <- as.matrix(sds)
  abundance_mat <- as.matrix(abundance_mat)
  D             <- ncol(means)
  
  if (is.null(rownames(means)) && !missing(species)) {
    if (length(species) != nrow(means)) {
      stop("length(species) must equal nrow(means) when means has no rownames.")
    }
    rownames(means) <- species
  }
  if (is.null(rownames(means))) {
    stop("means must have rownames, or species must be supplied.")
  }
  if (is.null(rownames(abundance_mat))) {
    rownames(abundance_mat) <- seq_len(nrow(abundance_mat))
  }
  if (is.null(colnames(abundance_mat))) {
    stop("abundance_mat must have species names as colnames.")
  }
  if (is.null(n_divisions) || length(n_divisions) != 1 || n_divisions < 2) {
    stop("n_divisions must be a single number >= 2.")
  }
  n_divisions <- as.integer(n_divisions)
  if (is.null(n_threads) || length(n_threads) != 1 || !is.finite(n_threads) || n_threads < 0) {
    stop("n_threads must be a single non-negative number; use 0 for all available OpenMP threads.")
  }
  n_threads <- as.integer(n_threads)
  
  if (is.null(trait_ranges)) trait_ranges <- rep(5, D)
  if (!is.list(trait_ranges)) {
    r <- trait_ranges
    if (length(r) == 1) r <- rep(r, D)
    if (length(r) != D) {
      stop("trait_ranges must be a list of length D, a scalar, or a numeric vector of length D.")
    }
    trait_ranges <- lapply(seq_len(D), function(d)
      c(min(means[, d] - r[d] * sds[, d]),
        max(means[, d] + r[d] * sds[, d])))
  }
  
  grid_list   <- lapply(seq_len(D), function(d)
    seq(trait_ranges[[d]][1], trait_ranges[[d]][2], length.out = n_divisions))
  cell_volume <- prod(vapply(grid_list, function(g) g[2] - g[1], numeric(1)))
  
  spp_in_means <- rownames(means)           # means有rownames
  spp_in_abd   <- colnames(abundance_mat)
  common_spp   <- intersect(spp_in_means, spp_in_abd)
  common_idx   <- match(common_spp, spp_in_means)  # 用数字索引，不依赖dimnames
  
  if (length(common_spp) == 0) {
    out <- rep(0, nrow(abundance_mat))
    names(out) <- rownames(abundance_mat)
    return(out)
  }
  
  means_sub <- means[common_idx, , drop = FALSE]
  sds_sub   <- sds[common_idx,   , drop = FALSE]   # 改为数字索引
  abd_sub   <- abundance_mat[, common_spp, drop = FALSE]
  
  rownames(means_sub) <- common_spp   # 确保有rownames
  rownames(sds_sub)   <- common_spp
  
  message(sprintf("Step 1: Building Sparse TPDs for %d species (n_divisions=%d)...",
                  nrow(means_sub), n_divisions))
  tictoc::tic("  TPD build")
  sparse_TPDs <- build_sparse_tpds_cpp(means_sub, sds_sub, grid_list,
                                       alpha, as.integer(n_divisions),
                                       n_threads = n_threads)
  tictoc::toc()
  
  message(sprintf("  Non-zero cells: %.0f avg / species",
                  mean(vapply(sparse_TPDs, length, integer(1)))))
  
  message("Step 2: Computing FRic via C++...")
  tictoc::tic("  FRic")
  result <- calc_fric_cpp(
    sparse_TPDs   = sparse_TPDs,
    abundance_mat = abd_sub,
    cell_volume   = cell_volume,
    basin_names   = rownames(abd_sub)
  )
  tictoc::toc()
  
  result
}


