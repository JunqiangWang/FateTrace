#'
#'
#'  Calculate the means and variance of the highly variable genes
CalculateFeatureStats <- function(seu, 
                                features=NULL,
                                assay = "RNA", 
                                slot = "data") {

  if(is.null(features)){
  if (is.null(VariableFeatures(seu))) {
    stop("No variable features found. Please run FindVariableFeatures() first.")
  }else{
    features <- VariableFeatures(seu)
  }
  }
  
  # Extract data matrix from specified slot and assay
  data_matrix <- GetAssayData(seu, assay = assay, slot = slot)
  feature_matrix <- data_matrix[features, ]
  
  feature_means <- Matrix::rowMeans(feature_matrix)
  feature_vars <- apply(feature_matrix, 1, var)
  
  feature_stats <- data.frame(
    gene = features,
    mean = feature_means,
    variance = feaure_vars,
    row.names = NULL
  )
  
  return(feature_stats)
}


#' 
#' 
GetPCLoadings <- function(seu, 
                          features=NULL,
                          assay = "RNA", 
                          reduction = "pca", 
                          dims = 1:50) {
  
  if (!reduction %in% names(seu@reductions)) {
    stop(paste("Reduction", reduction, "not found in the Seurat object. Run RunPCA() first."))
  }
  
  if(is.null(features)){
    if (is.null(VariableFeatures(seu))) {
      stop("No variable features found. Please run FindVariableFeatures() first.")
    }else{
      features <- VariableFeatures(seu)
    }
  
  # Extract PCA loadings (genes x PCs)
  loadings <- Loadings(seu[[reduction]])
  
  # Subset to HVFs and selected PCs
  if (length(dims) != ncol(loadings)) {
    warning("Some specified dimensions (PCs) are not available in the loadings matrix.")
  }
  
  feature_loadings <- loadings[rownames(loadings) %in% features, dims, drop = FALSE]
  
  return(feature_loadings)
}

  
  
#'
#'
#'
FilterQueryCounts <- function(seu.query, 
                              seu.reference, 
                              features_reference=NULL,
                              assay = "RNA", 
                              slot = "data") {
  # 
  if(is.null(features_reference)){
    if (is.null(VariableFeatures(seu.reference))) {
      stop("No variable features found.")
    }else{
      features_reference <- VariableFeatures(seu.reference)
    }

  # Extract the data matrix from the query object
  query_data <- GetAssayData(seu.query, assay = assay, slot = slot)
  
  common_features <- intersect(rownames(query_data), features_reference)
  
  if (length(common_features) == 0) {
    stop("No common variable features found between the reference and query data.")
  }
  
  filtered_data <- query_data[common_features, , drop = FALSE]
  
  return(filtered_data)
}


  
#'
#'
#'
ScaleCountsByReference <- function(query_matrix, 
                                   ref_means, 
                                   ref_vars,
                                   ) {
  # Ensure names match
  common_genes <- intersect(rownames(query_matrix), names(ref_means))
  common_genes <- intersect(common_genes, names(ref_vars))
  
  if (length(common_genes) == 0) {
    stop("No overlapping features found!")
  }
  
  # Subset everything to common genes
  query_matrix <- query_matrix[common_genes, , drop = FALSE]
  ref_means <- ref_means[common_genes]
  ref_vars <- ref_vars[common_genes]
  
  # Scale using reference mean and variance
  scaled_matrix <- sweep(query_matrix, 1, ref_means, "-")
  scaled_matrix <- sweep(scaled_matrix, 1, sqrt(ref_vars), "/")
  
  scaled_matrix <- as.matrix(scaled_matrix)
  
  return(scaled_matrix)
  
}



#' @description
#' 
#' @param seu.query The query seurat object
#' @param seu.reference The reference seurat object
#' @export 
ScaleDataByReference<-function(seu.query,
                               seu.reference,
                               features_query=NULL,
                               features_reference=NULL
                              ){
  # 
  feature_stats<-CalculateFeatureStats(seu=seu.reference,
                                       features=features_query
                                       )
  
  PCloadings_reference<-GetPCLoadings(seu=seu.reference)
  
  # HVF_PCloading_ref
  tmp<-identical(rownames(PCloadings_reference), feature_stats$gene)
  
  if(isFALSE(tmp)){
    warning("PCA of reference dataset was not performed based on highly varible genes.")
  }
  
  
  query.data<-FilterQueryCounts(seu.query=seu.query, 
                                seu.reference=seu.reference, 
                                features_reference=NULL,
                                assay = "RNA", 
                                slot = "data")
  
  # If not Assay5, set as Assay5

  #
  ref_means<-feature_stats$mean
  names(ref_means)<-feature_stats$gene
  
  ref_vars<-feature_stats$variance
  names(ref_vars)<-feature_stats$gene
  
  
  query.data.scaled<-ScaleCountsByReference(query_matrix=query.data, ref_means=ref_means, ref_vars=ref_vars)
  
  seu.query <- SetAssayData(seu.query, assay = "RNA", layer = "scale.data", new.data = query.data.scaled)
  
  return(seu.query)

}




#' @description
#' This function impute the PC scores using the reference dataset
#' @param seu.query The query seurat object
#' @param seu.reference The reference seurat object
#' @export 
ImputePCs<-function(seu.query,
                    seu.reference,
                    reduction="pca",
                    reduction.name="pca.imputed"
                    ){
  
  ref_pca_model <- Embeddings(seu.reference, reduction = reduction)
  ref_pca_loadings <- Loadings(seu.reference[[reduction]])  # genes × PCs
  
  query.data.scaled <- GetAssayData(seu.query, assay = "RNA", slot = "scale.data")
  
  # 
  # Subset to common genes
  common_genes <- intersect(rownames(ref_pca_loadings), rownames(query.data.scaled))
  query.data.scaled <- query.data.scaled[common_genes, , drop = FALSE]
  ref_pca_loadings <- ref_pca_loadings[common_genes, , drop = FALSE]
  
  # Project query onto reference PC loadings
  query_pca_embeddings <- t(ref_pca_loadings) %*% query.data.scaled  # PCs × cells
  query_pca_embeddings <- t(query_pca_embeddings)  # cells × PCs
  
  # Create a DimReduc object with the projected embedding
  query_pca <- CreateDimReducObject(
    embeddings = as.matrix(query_pca_embeddings),
    key = "PC_",
    assay = "RNA"
  )
  
  # Assign to query Seurat object
  seu.query[[reduction.name]] <- query_pca
  
  return(seu.query)
  
}



RunFateTraceUMAP <- function(
    seu,
    reduction = "pca",
    pc.used = NULL,
    umap.name = "umap",
    metric = c("cosine", "euclidean"),
    n_neighbors = 30,
    min_dist = 0.3,
    repulsion.strength = 0.1,      # Seurat default is 1.0
    seed = 777,                    # your seed
    umap_model_file = "umap_model.umap",
    save_umap_model_file = TRUE
){
  
  stopifnot(inherits(seu, "Seurat"))
  metric <- match.arg(metric)
  
  emb <- Embeddings(seu, reduction = reduction)
  if (is.null(pc.used)) {
    pc.used <- seq_len(min(50L, ncol(emb)))  # Seurat default dims 1:50
  }
  
  pcs_avail <- ncol(emb)
  pc.keep <- sort(unique(pc.used[pc.used >= 1 & pc.used <= pcs_avail]))
  if (length(pc.keep) == 0L) stop("No valid PCs in 'pc.used' after sanitization.")
  
  # Construct X and preserve its rownames for later
  X <- as.matrix(emb[, pc.keep, drop = FALSE])
  storage.mode(X) <- "double"
  if (anyNA(X) || any(!is.finite(X))) {
    stop("Input matrix contains NA/Inf. Please remove or impute before UMAP.")
  }
  
  # Also set the R RNG, just in case
  set.seed(seed)
  
  umap_model <- uwot::umap(
    X = X,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = metric,
    spread = 1.0,
    local_connectivity = 1L,
    repulsion_strength = repulsion.strength,
    n_epochs = 500,
    negative_sample_rate = 5L,
    set_op_mix_ratio = 1.0,
    n_components = 2L,
    init = "spectral",
    fast_sgd = FALSE,          # ← deterministic
    n_threads = 1,             # ← deterministic
    n_sgd_threads = 1,         # ← deterministic
    verbose = TRUE,
    ret_model = TRUE,
    y = NULL,
    seed = seed
  )
  
  if (isTRUE(save_umap_model_file)) {
    uwot::save_uwot(umap_model, file = umap_model_file)
  }
  
  umap_coords <- umap_model$embedding
  # Keep the SAME row order/names as X (cells)
  rownames(umap_coords) <- rownames(X)
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  
  seu[[umap.name]] <- CreateDimReducObject(
    embeddings = umap_coords,
    key = "UMAP_",
    assay = DefaultAssay(seu)
  )
  seu@misc$umap_model <- umap_model
  seu
}



#' @description
#' This function imputes UMAP
#' @param seu seurat object
#' @param uamp_model umap model
#' @export
ImputeUMAP<-function(seu,
                     umap_model,
                     pcs,
                     reduction.name="pca",
                     umap_name = "umap.imputed",
                     seed=777
                     ){
  
  require(uwot)
  
  query_pcs <- Embeddings(seu, reduction = reduction.name)
  query_pcs<- query_pcs[, pcs]
  
  umap <- uwot::umap_transform(query_pcs, 
                               umap_model,
                               seed=seed
                               )
  
  umap_coords <- umap
  rownames(umap_coords) <- rownames(seu@meta.data)
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  
  umap_dimred <- CreateDimReducObject(
    embeddings = umap_coords,
    key = "UMAP_",
    assay = DefaultAssay(seu)
  )
  
  seu[[umap_name]] <- umap_dimred
  
  return(seu)
  
}



