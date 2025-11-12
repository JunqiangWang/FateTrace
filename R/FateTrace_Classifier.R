GetCorKNN <- function(expr_matrix, k = 20) {
  # expr_matrix: cells as rows, genes as columns
  # k: number of nearest neighbors
  
  if (!is.matrix(expr_matrix) && !is.data.frame(expr_matrix)) {
    stop("expr_matrix must be a matrix or data frame with cells as rows and genes as columns.")
  }
  
  if (k >= nrow(expr_matrix)) {
    stop("k must be less than the number of cells (rows in expr_matrix).")
  }
  
  cor_matrix <- cor(t(expr_matrix), method = "pearson", use = "pairwise.complete.obs")
  
  cor_dist <- 1 - cor_matrix
  
  diag(cor_dist) <- Inf
  
  nn.index <- t(apply(cor_dist, 1, function(x) order(x)[1:k]))
  nn.dist <- t(apply(cor_dist, 1, function(x) sort(x)[1:k]))
  
  return(list(nn.index = nn.index, nn.dist = nn.dist))
}



GetCosineKNN <- function(expr_matrix, k = 20) {
  # Normalize rows (cells) to unit vectors
  normalize <- function(x) x / sqrt(sum(x^2) + 1e-10)  # Add epsilon to avoid division by zero
  norm_expr <- t(apply(expr_matrix, 1, normalize))
  
  cosine_sim <- norm_expr %*% t(norm_expr)
  
  cosine_dist <- 1 - cosine_sim
  
  nn.index <- matrix(0, nrow = nrow(cosine_dist), ncol = k)
  nn.dist <- matrix(0, nrow = nrow(cosine_dist), ncol = k)
  
  for (i in 1:nrow(cosine_dist)) {
    # Exclude self by setting distance to Inf
    dists <- cosine_dist[i, ]
    dists[i] <- Inf
    
    top_k <- order(dists)[1:k]
    nn.index[i, ] <- top_k
    nn.dist[i, ] <- dists[top_k]
  }
  
  return(list(nn.index = nn.index, nn.dist = nn.dist))
}






#'
#'
#'
#'
KNN <- function(pc_matrix, 
                      dims = NULL,
                      label, 
                      k = 20, 
                      max_iter = 1, 
                      tol = 1e-4) {
  
  library(FNN)          # For get.knn
  library(matrixStats)  # For any future statistical operations
  
  
  if (is.null(dims)) {
    dims <- 1:ncol(pc_matrix)
  }
  
  selected_pcs <- pc_matrix[, dims]
  
  
  
  
  message("Compute KNN graph")
  knn_result <- get.knn(selected_pcs, k = k)
  neighbor_indices <- knn_result$nn.index
  message("The KNN graph is computated")
  
  

  
  current_stage <- label_vector
  for (i in 1:max_iter) {
    new_stage <- sapply(1:nrow(neighbor_indices), function(cell) {
      neighbor_stages <- current_stage[neighbor_indices[cell, ]]
      new<-names(sort(table(neighbor_stages), decreasing = TRUE))[1]
      return(new)
    })
    
    
    change_ratio <- mean(new_stage != current_stage)
    message(sprintf("Iteration %d: %.4f of labels changed.", i, change_ratio))
    
    if (change_ratio < tol) {
      message(sprintf("Converged after %d iterations.", i))
      break
    }
    
    current_stage <- new_stage
  }
  
  if (change_ratio >= tol) {
    message("Reached maximum iterations without full convergence.")
  }
  
  return(current_stage)
}




# example
# reference_pc_matrix<-Embeddings(seu, reduction="ref.pca")
# pc.used<-c(1:5) 
# reference_label_vector<-seu@meta.data$age 
# query_pc_matrix<-Embeddings(seu, reduction="ref.pca")
# updated_label_vector <-FateTraceKNNImputation(
#   reference_pc_matrix=reference_pc_matrix,
#   query_pc_matrix=query_pc_matrix,
#   dims = pc.used,
#   reference_stage_vector=reference_stage_vector, 
#   k = 10, 
#   max_iter = 1)

#'
#'
FateTraceKNNImputation <- function(
    reference_pc_matrix,
    query_pc_matrix,
    dims = NULL,
    reference_label_vector, 
    query_label_vector = NULL, 
    k = 10, 
    max_iter = 1){
  
  # The reference-based analysis
  library(FNN)          # For get.knn
  library(matrixStats)  # For any future statistical operations
  
  # Step 1: Filter out PCs
  if (is.null(dims)) {
    dims <- 1:ncol(reference_pc_matrix)
  }
  
  if(ncol(reference_pc_matrix)!=ncol(query_pc_matrix)
  ){
    stop("Error: Dimensions of PCs from the reference_pc_matrix and query_pc_matrix do not match.")
  }
  
  # 
  reference_pc_matrix <- reference_pc_matrix[, dims]
  query_pc_matrix <- query_pc_matrix[, dims]
  
  # Step 2: Compute KNN graph
  message("Compute KNN graph")
  knn_result <- get.knnx(data=reference_pc_matrix, query=query_pc_matrix, k = k)
  
  neighbor_indices <- knn_result$nn.index
  message("The KNN graph is computated")
  
  # Step 3: Iterative label updating
  new_stage <- sapply(1:nrow(neighbor_indices), function(cell) {
    neighbor_stages <- reference_label_vector[neighbor_indices[cell, ]]
    new<-names(sort(table(neighbor_stages), decreasing = TRUE))[1]
    return(new)
  })
  
  # Check for convergence
  if(!is.null(query_label_vector)){
    original_stage <- query_label_vector
    change_ratio <- mean(new_stage != original_stage)
    message(sprintf("Iteration: %.4f of labels changed.", change_ratio))
  }
  
  # if (change_ratio >= tol) {
  #  message("Reached maximum iterations without full convergence.")
  # }
  
  current_stage <- new_stage
  return(current_stage)
}




FateTraceKNNCosine <- function(expr_matrix, 
                            dims = NULL,
                            label_vector, 
                            k = 20, 
                            max_iter = 100, 
                            tol = 1e-4) {
  
  
  knn_result <- GetCosineKNN(expr_matrix, k = k)
  neighbor_indices <- knn_result$nn.index
  

  current_stage <- label_vector
  
  for (i in 1:max_iter) {
    new_stage <- sapply(1:nrow(neighbor_indices), function(cell) {
      neighbor_stages <- label_vector[neighbor_indices[cell, ]]
      new<-names(sort(table(neighbor_stages), decreasing = TRUE))[1]
      return(new)
    })
    
    change_ratio <- mean(new_stage != current_stage)
    message(sprintf("Iteration %d: %.4f of labels changed.", i, change_ratio))
    
    if (change_ratio < tol) {
      message(sprintf("Converged after %d iterations.", i))
      break
    }
    
    current_stage <- new_stage
  }
  
  if (change_ratio >= tol) {
    message("Reached maximum iterations without full convergence.")
  }

  
  return(current_stage)
}














