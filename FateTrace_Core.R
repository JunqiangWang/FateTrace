




#'
#'
#'
DownsampleSparseMatrix <- function(mat, depth = 3000, ncores = 8, seed=777, verbose = TRUE) {
  
  set.seed(seed)
  
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package required")
  if (!inherits(mat, "dgCMatrix"))
  {
    message("Input must be dgCMatrix")
    message("Converting matrix to dgCMatrix")
    mat <- as(mat, "dgCMatrix")
  } 
  
  message("Data quality checking...")

  n_genes <- nrow(mat)
  n_cells <- ncol(mat)
  orig_rownames <- rownames(mat)
  
  
  
  depth <- rep(depth, length.out = n_cells)
  valid_cols <- which(Matrix::colSums(mat) > 0 & depth > 0)
  
  if (verbose) {
    message("Processing ", length(valid_cols), "/", n_cells, " non-zero columns")
    if (!is.null(orig_rownames)) {
      message("Row names detected: ", length(orig_rownames), " vs ", n_genes, " rows")
    }
  }
  
  
  # Sparse matrix components
  p <- mat@p
  i <- mat@i  # 0-based indices
  x <- mat@x
  
  
  # Parallel processing function
  process_column <- function(j) {
    if (j %in% valid_cols) {
      col_start <- p[j] + 1
      col_end <- p[j + 1]
      idx <- i[col_start:col_end]  # 0-based indices
      vals <- x[col_start:col_end]
      
      sampled <- rmultinom(1, depth[j], prob = vals/sum(vals))[, 1]
      nonzero <- sampled > 0
      

      
      Matrix::sparseVector(
        x = sampled[nonzero],
        i = idx[nonzero] + 1L,  # Convert to 1-based indexing
        length = n_genes
      )
    } else {
      Matrix::sparseVector(integer(), integer(), n_genes)
    }
  }
  
  # Parallel execution
  if (verbose) message("Downsampling with ", ncores, " cores...")
  
  result <- if (requireNamespace("pbmcapply", quietly = TRUE) & n_cells<100001) {
    pbmcapply::pbmclapply(seq_len(n_cells), process_column, 
                          mc.cores = ncores, ignore.interactive = !verbose)
  } else {
    parallel::mclapply(seq_len(n_cells), process_column, mc.cores = ncores)
  }
  
  
  # Assemble final matrix
  if (verbose) message("Constructing sparse matrix...")
  # Combine 
  i_list <- lapply(seq_along(result), function(j) result[[j]]@i)
  x_list <- lapply(seq_along(result), function(j) result[[j]]@x)
  col_ptrs <- lapply(seq_along(result), function(j) rep(j, length(result[[j]]@i)))
  
  i_all <- unlist(i_list)
  x_all <- unlist(x_list)
  j_all <- unlist(col_ptrs)
  
  
  
  out <- as(Matrix::sparseMatrix(
    i = i_all,
    j = j_all,
    x = x_all,
    dims = c(n_genes, n_cells),
    dimnames = list(orig_rownames, colnames(mat))
  ),  "dgCMatrix")
  

  if (!is.null(orig_rownames) && length(orig_rownames) == n_genes) {
    rownames(out) <- orig_rownames
  } else {
    warning("Row names discarded due to length mismatch", call. = FALSE)
    rownames(out) <- NULL
  }
  colnames(out) <- colnames(mat)
  
  if (verbose) {
    message("Final matrix: ", 
            format(utils::object.size(out), units = "auto"),
            "\nDimensions: ", paste(dim(out), collapse = " x "))
  }
  
  return(out)
}



#'
#'
#'
DownsampleSparseMatrixBatched <- function(mat, depth = 3000, batch_size = 10000, ...) {
  message("Constructing final merged sparse matrix...")
  n_cells <- ncol(mat)
  batches <- split(seq_len(n_cells), ceiling(seq_len(n_cells)/batch_size))
  
  results <- lapply(batches, function(batch) {
    DownsampleSparseMatrix(mat[, batch], depth = depth, ...)
  })
  
  m<-do.call(cbind, results)
  return(m)
  message("The final merged sparse matrix is constructed")
}





#' @param seu seurat object
#' @param depth downsampling depth
#' @param ncores number of cores used
#' @param batch_size batch size 
#' @export
DownsampleMultinomialFastSeu<-function(seu,
                                       depth=3000,
                                       ncores=4,
                                       batch_size=10000){
  
  seu[['RNA']] <-as(seu[['RNA']], Class="Assay")
  
  counts <- GetAssayData(seu, slot = "counts")
  
  meta.data<-seu@meta.data
  
  counts <- DownsampleSparseMatrixBatched(mat=counts, depth = depth, ncores=ncores, batch_size = batch_size)
  
  seu<-CreateSeuratObject(counts=counts, meta.data = meta.data)
  
  seu@meta.data$nCount_RNA<-colSums(counts)
  
  return(seu)
}






#' @param seu seurat object
#' @param var variable to downsample
#' @param size.n downsampled size
#' @param seed seed for random sampling
#' @export
DownsampleCells <- function(seu, var = NULL, size.n = 1000, seed = 777) {
  library(dplyr)
  
  if (!is.null(seed)) set.seed(seed)
  
  meta <- seu@meta.data
  meta$cell_id <- rownames(meta)
  
  if (is.null(var)) {
    picked <- meta %>%
      mutate(.rand = runif(n())) %>%
      arrange(.rand) %>%
      slice_head(n = min(size.n, n())) %>%
      pull(cell_id)
  } else {
    missing <- setdiff(var, colnames(meta))
    if (length(missing)) stop("Columns not found in metadata: ", paste(missing, collapse = ", "))
    
    picked <- meta %>%
      group_by(across(all_of(var))) %>%
      mutate(.rand = runif(n())) %>%
      arrange(.rand, .by_group = TRUE) %>%
      slice_head(n = size.n) %>%  
      ungroup() %>%
      pull(cell_id)
  }
  
  subset(seu, cells = picked)
}




#' Downsample cells with optional per–cell-type targets
#'
#' @param seu Seurat object
#' @param var Character scalar or vector. Column(s) in meta.data to group by.
#'            If using `cell_number_vector`, `var` must be a single column (the cell-type column).
#' @param cell_number_vector Named numeric vector giving the number of cells to sample
#'            for each specified level of `var`. Names must match the entries in `seu@meta.data[[var]]`.
#'            Levels not listed will not be sampled (i.e., 0) unless you use the fallback `size.n`.
#' @param size.n Integer. Fallback sample size (global if `var` is NULL, or per-group if `var` is given).
#'               Ignored when `cell_number_vector` is provided (except to break ties if you combine strategies).
#' @param seed Random seed for reproducibility.
#' @export
DownsampleCellTypes <- function(seu,
                                var = NULL,
                                cell_number_vector = NULL,
                                size.n = 1000,
                                seed = 777) {
  library(dplyr)
  if (!is.null(seed)) set.seed(seed)
  
  meta <- seu@meta.data
  meta$cell_id <- rownames(meta)
  
  # 
  if (is.null(var) && is.null(cell_number_vector)) {
    picked <- meta %>%
      mutate(.rand = runif(n())) %>%
      arrange(.rand) %>%
      slice_head(n = min(size.n, n())) %>%
      pull(cell_id)
    
    return(subset(seu, cells = picked))
  }
  
  # Validate var columns
  if (!is.null(var)) {
    missing <- setdiff(var, colnames(meta))
    if (length(missing)) stop("Columns not found in metadata: ", paste(missing, collapse = ", "))
  }
  
  #
  if (!is.null(cell_number_vector)) {
    if (is.null(var) || length(var) != 1) {
      stop("When using 'cell_number_vector', 'var' must be a single column representing the cell type.")
    }
    
    ct_col <- var
    levels_in_data <- unique(as.character(meta[[ct_col]]))
    requested_types <- names(cell_number_vector)
    
    if (is.null(requested_types) || anyNA(requested_types)) {
      stop("'cell_number_vector' must be a *named* numeric vector (e.g., c(B_cell = 300, T_cell = 500)).")
    }
    
    # 
    not_found <- setdiff(requested_types, levels_in_data)
    if (length(not_found)) {
      warning("These requested cell types were not found in '", ct_col, "': ",
              paste(not_found, collapse = ", "))
    }
    
    # 
    valid_types <- intersect(requested_types, levels_in_data)
    if (!length(valid_types)) {
      stop("None of the names in 'cell_number_vector' match levels in '", ct_col, "'.")
    }
    
    # 
    sampled_list <- lapply(valid_types, function(tp) {
      n_req <- as.integer(cell_number_vector[[tp]])
      if (is.na(n_req) || n_req < 0) {
        stop("Invalid requested count for type '", tp, "'. Must be a non-negative integer.")
      }
      
      df <- meta %>% filter(.data[[ct_col]] == tp)
      n_take <- min(n_req, nrow(df))
      
      if (n_take < n_req) {
        warning("Requested ", n_req, " cells for '", tp,
                "' but only ", nrow(df), " available. Taking ", n_take, ".")
      }
      
      # random sample without replacement
      if (n_take > 0) {
        df %>% slice_sample(n = n_take) %>% pull(cell_id)
      } else {
        character(0)
      }
    })
    
    picked <- unlist(sampled_list, use.names = FALSE)
    
    return(subset(seu, cells = picked))
  }
  
  # 
  picked <- meta %>%
    group_by(across(all_of(var))) %>%
    mutate(.rand = runif(n())) %>%
    arrange(.rand, .by_group = TRUE) %>%
    slice_head(n = size.n) %>%           # up to size.n per group
    ungroup() %>%
    pull(cell_id)
  
  subset(seu, cells = picked)
}








