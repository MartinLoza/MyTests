
FindNumberPC <- function(object = NULL, dim2Test = NULL, returnBIC = FALSE){
  UseMethod("FindNumberPC")
}



FindNumberPC.Seurat <- function(object = NULL, dim2Test = NULL, returnBIC = FALSE, ...){
  
  features <- VariableFeatures(object)
  
  # We use the Seurat scaling to get comparable results. If we manually scale the data we could get similar but different results.
  scaData <- t(object[["RNA"]]@scale.data)
  
  if(length(scaData) == 0)
    stop("The Seurat object must be scaled.")
  
  # if("pca" %in% names(object@reductions)){
  #   
  #   # Check that the calculated PCA contains all of the eigen values of the covariance matrix
  #   if(length(object[["pca"]]@stdev) != min(length(features), nrow(object[["pca"]]))){
  #     scaData <- t(object@assays$RNA@scale.data)
  #   }else{
  #     eigVal <- (object[["pca"]]@stdev)^2
  #   }
  # }else{
  #   
  #   scaData <- scale(t(as.matrix(GetAssayData(object, assay = "RNA", slot = "data")[features,])), center = TRUE, scale = TRUE)
  # }
  
  if(!exists("eigVal"))
    eigVal <- eigen(cov(scaData))$values
  
  if(is.null(dim2Test))
    dim2Test <- min(dim(scaData))
  
  bicPCA <- -getBicPca(eig = eigVal, N = ncol(scaData), K = dim2Test)
  selPCA <- which.min(bicPCA)
  
  if(returnBIC)
    output <- list(BIC = bicPCA, selPCA = selPCA)
  else
    output <- selPCA
  
  return(output)
  
}