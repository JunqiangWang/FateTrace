


#'
#'
#'
#'
ConnectComponentsKNN <- function(adj, X, dist.method = "euclidean", k.bridges = k.nn) {
  stopifnot(nrow(adj) == ncol(adj))
  labs <- rownames(adj); stopifnot(identical(labs, colnames(adj)))
  stopifnot(k.bridges >= 1L)
  
  # 
  X <- X[labs, , drop = FALSE] 
  
  # 
  adj <- methods::as(adj, "dgCMatrix")
  
  # 
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE)
  comps <- igraph::components(g)
  if (comps$no <= 1L) {
    adj2 <- (adj + t(adj)) / 2
    Matrix::diag(adj2) <- 0
    return(methods::as(adj2, "dgCMatrix"))
  }
  
  mem <- comps$membership; names(mem) <- igraph::V(g)$name
  comp_list <- split(names(mem), mem)
  
  giant_id  <- as.integer(names(which.max(vapply(comp_list, length, integer(1)))))
  connected <- comp_list[[giant_id]]
  others    <- setdiff(seq_along(comp_list), giant_id)
  
  bridge_i <- integer(0); bridge_j <- integer(0); bridge_x <- numeric(0)
  use_rann <- (dist.method == "euclidean") && requireNamespace("RANN", quietly = TRUE)
  
  for (cid in others) {
    S <- comp_list[[cid]]
    if (length(S) == 0L || length(connected) == 0L) next
    
    eff_k <- min(k.bridges, length(S) * length(connected))  # cannot add more unique pairs than this
    
    if (use_rann) {
      # 
      kq <- min(k.bridges, length(connected))
      nn <- RANN::nn2(
        data  = X[connected, , drop = FALSE],
        query = X[S, , drop = FALSE],
        k = kq
      )
      # 
      ds  <- nn$nn.dists
      idx <- nn$nn.idx
      ord <- order(as.vector(ds), na.last = NA)
      take <- head(ord, eff_k)
      rc <- arrayInd(take, .dim = dim(ds))
      # 
      u <- S[ rc[, 1] ]
      v <- connected[ idx[cbind(rc[, 1], rc[, 2])] ]
      d <- ds[ cbind(rc[, 1], rc[, 2]) ]
    } else {
      # 
      D <- as.matrix(proxyC::dist(X[S, , drop = FALSE],
                                  X[connected, , drop = FALSE],
                                  method = dist.method))
      ord <- order(as.vector(D), na.last = NA)
      take <- head(ord, eff_k)
      rc <- arrayInd(take, .dim = dim(D))
      u <- S[ rc[, 1] ]
      v <- connected[ rc[, 2] ]
      d <- D[ cbind(rc[, 1], rc[, 2]) ]
    }
    
    iu <- match(u, labs); iv <- match(v, labs)
    # 
    bridge_i <- c(bridge_i, iu, iv)
    bridge_j <- c(bridge_j, iv, iu)
    bridge_x <- c(bridge_x, d,  d)
    
    # 
    connected <- c(connected, S)
  }
  
  # 
  if (length(bridge_i)) {
    B <- Matrix::sparseMatrix(
      i = bridge_i, j = bridge_j, x = bridge_x,
      dims = dim(adj), dimnames = dimnames(adj), symmetric = FALSE
    )
    adj2 <- adj + B
  } else {
    adj2 <- adj
  }
  adj2 <- (adj2 + t(adj2)) / 2
  Matrix::diag(adj2) <- 0
  methods::as(adj2, "dgCMatrix")
}





#' @param seu seurat object
#' @param return_pseudotime whether to return pseudotime
#' @export
InferPseudotime <- function(seu, 
                               reduction=NULL,
                               dims=NULL,
                               dist.method = "euclidean",
                               k.nn = 30,
                               block.size=1e6L,
                               seed=777,
                               k.diffuse= 10, # 5-20
                               return_pseudotime = FALSE) {
  #
  require(Matrix)
  require(igraph)
  require(dplyr)
  require(tibble)
  require(RSpectra)
  require(FNN)
  set.seed(seed)
  
  # 
  if (!requireNamespace("igraph", quietly = TRUE)) stop("Package 'igraph' is required.")
  library(igraph)
  
  if(is.null(reduction)){
    stop("Error: please provide the reduction used for calculating the distance")
  }
  

  
  embedding <- Embeddings(seu, reduction = reduction)
  cell_names <- rownames(embedding)
  
  if(!is.null(dims)){
    embedding<-embedding[, dims]
  }
  
  root_cell <- seu@misc$root_cell
  
  

  
  message("calculate distance matrix...")
  
  k.nn <- 50
  knn <- get.knn(embedding, k = k.nn)
  nn_index <- knn$nn.index
  nn_dist  <- knn$nn.dist
  # Build a sparse kNN graph if needed:
  library(Matrix)
  i <- rep(seq_len(nrow(embedding)), each = k.nn)
  j <- as.vector(t(nn_index))
  x <- as.vector(t(nn_dist))
  adj <- sparseMatrix(i = i, j = j, x = x, dims = c(nrow(embedding), nrow(embedding)))
  
  rownames(adj) <- cell_names
  colnames(adj) <- cell_names
  
  message("distance matrix is calculated")
  

  
  
  adj <- (adj + t(adj)) / 2
  diag(adj) <- 0
  
  message("create connected graph...")
  
  
  adj <- ConnectComponentsKNN(adj, X=embedding, dist.method = "euclidean", k.bridges = k.nn)
  
  message("connected graph is created...")
  

  
  
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE)
  comps <- components(g)
  giant <- which(comps$membership == which.max(comps$csize))
  keep <- names(V(g))[giant]
  
  adj <- adj[keep, keep, drop = FALSE]
  
  g <- induced_subgraph(g, keep)
  
 
  
  
  if (!root_cell %in% keep) {
    stop(sprintf("root_id '%s' is not in the largest connected component.", root_cell))
  }
  
  n <- nrow(adj)

  d <- Matrix::rowSums(adj)
  eps <- 1e-12
  d[d < eps] <- eps
  Dinv <- Diagonal(x = 1 / d)
  P <- Dinv %*% adj  # sparse
  
  
  k_eff <- min(k.diffuse + 1, n - 1)  # safety
  eig <- RSpectra::eigs(P, k = k_eff, which = "LR")  # largest real parts
  vals <- Re(eig$values)
  vecs <- Re(eig$vectors)
  
  #
  o <- order(vals, decreasing = TRUE)
  vals <- vals[o]
  vecs <- vecs[, o, drop = FALSE]
  
  # 
  vals_nt <- vals[-1]
  vecs_nt <- vecs[, -1, drop = FALSE]
  
  k_use <- min(k.diffuse, ncol(vecs_nt))
  

  diff_time <- 1L
  

  
  vals_k <- as.numeric(Re(vals_nt[seq_len(k_use)]))
  vecs_k <- as.matrix(Re(vecs_nt[, seq_len(k_use), drop = FALSE]))
  
  psi <- sweep(vecs_k, 2, vals_k^diff_time, `*`)  # n x k
  rownames(psi) <- rownames(adj)
  
  
  root_idx <- match(root_cell, rownames(psi))
  root_vec <- psi[root_idx, , drop = FALSE]
  pt <- sqrt(rowSums((psi - matrix(root_vec, nrow(psi), ncol(psi), byrow = TRUE))^2))
  

  pt <- (pt - min(pt)) / (max(pt) - min(pt) + 1e-12)

  
  names(pt) <- rownames(psi)
  

  pseudotime_list<-list(
    pseudotime = pt,          # named numeric vector ∈ [0,1]
    diffusion_coordinates = psi,     # n x k matrix
    eigenvalues = vals_k,            # length k
    kept_cells = rownames(adj)       # cells in the connected component
  )
  
  seu@misc$pseudotime_list<-pseudotime_list
  
  
  #
  pseudotime<-rep(NA, length=nrow(seu@meta.data))
  names(pseudotime) <- rownames(seu@meta.data)
  
  tmp<-match(names(pseudotime_list$pseudotime), names(pseudotime))
  pseudotime[tmp]<-pseudotime_list$pseudotime
  seu@meta.data$pseudotime<-pseudotime
  
  seu@misc$pseudotime<-pseudotime
  
  if(!isTRUE(return_pseudotime)){
    return(seu)
  }else{
    return(pseudotime)
  }
}





#' @param seu seurat object
#' @param reduction reduction
#' @param show_mst show mst path
#' @export
PlotPseudotime <- function(seu,
                           reduction = "umap",
                           show_mst = FALSE,
                           na.value.color = "grey",
                           use.raster =TRUE
                           ) {
  
  library(ggplot2)
  library(viridis)
  
  # 
  umap_coords <- Embeddings(seu, reduction = reduction)
  pseudotime <- seu@misc$pseudotime
  
  # 
  plot_df <- data.frame(
    UMAP_1 = umap_coords[, 1],
    UMAP_2 = umap_coords[, 2],
    pseudotime = pseudotime[rownames(umap_coords)]
  )
  

  
  if(isTRUE(use.raster)){
    
    require(ggrastr)
    
    if (show_mst == FALSE) {
      p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
        geom_point_rast(size = 0.1, alpha = 0.8, raster.dpi = 300) +   # rasterized points
        scale_color_viridis_c(option = "plasma", na.value = "gray") +   # gray for NA
        theme_minimal() +
        labs(title = "", color = "Pseudotime") + gg.theme2
      return(p)
      
    } else {
      mst_graph <- seu@misc$mst_graph
      
      mst_edges <- igraph::as_data_frame(mst_graph, what = "edges")
      
      edge_df <- data.frame(
        x    = umap_coords[mst_edges$from, 1],
        y    = umap_coords[mst_edges$from, 2],
        xend = umap_coords[mst_edges$to,   1],
        yend = umap_coords[mst_edges$to,   2]
      )
      
      p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
        geom_point_rast(size = 0.1, alpha = 0.9, raster.dpi = 300) +    # rasterized points
        geom_segment(
          data = edge_df,
          aes(x = x, y = y, xend = xend, yend = yend),
          inherit.aes = FALSE,
          color = "black", alpha = 0.5, linewidth = 0.5                   # keep MST as vector
        ) +
        scale_color_viridis_c(option = "plasma", na.value = na.value.color) +  # gray for NA
        theme_minimal() +
        labs(title = "", color = "Pseudotime") 
  
    }
    
    
  }else{
  
  
  if (show_mst == FALSE) {
    
    p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
      geom_point(size = 0.1, alpha = 0.8) +
      scale_color_viridis_c(option = "plasma", na.value = "gray") +   ### gray for NA
      theme_minimal() +
      labs(title = "", color = "Pseudotime") + gg.theme2
    
    return(p)
    
  } else {
    
    mst_graph <- seu@misc$mst_graph
    
    mst_edges <- igraph::as_data_frame(mst_graph, what = "edges")
    
    edge_df <- data.frame(
      x = umap_coords[mst_edges$from, 1],
      y = umap_coords[mst_edges$from, 2],
      xend = umap_coords[mst_edges$to, 1],
      yend = umap_coords[mst_edges$to, 2]
    )
    
    p<-ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
      geom_point(size = 0.1, alpha = 0.9) +
      geom_segment(data = edge_df,
                   aes(x = x, y = y, xend = xend, yend = yend),
                   inherit.aes = FALSE,
                   color = "black", alpha = 0.5, linewidth = 0.5) +
      scale_color_viridis_c(option = "plasma", na.value = na.value.color) +   ### gray for NA
      theme_minimal() +
      labs(title = "", color = "Pseudotime") 
  }
    
  }
  
  p<- p+gg.theme2
  return(p)
  
}


