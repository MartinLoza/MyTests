
FindNumberPC <- function(object = NULL, dim2Test = NULL, returnBIC = FALSE){
  UseMethod("FindNumberPC")
}



FindNumberPC.Seurat <- function(object = NULL, dim2Test = NULL, returnBIC = FALSE){
  
  if("pca" %in% names(object@reductions)){
    
    # Check that the calculated PCA contains all of the eigen values of the covariance matrix
    if(length(object[["pca"]]@stdev) != min(length(VariableFeatures(object)), nrow(object[["pca"]]))){
      scaData <- x@assays$RNA@scale.data
    }else{
      eigVal <- (x[["pca"]]@stdev)^2
    }
  }else{
    scaData <- scale(t(GetAssayData(object, assay = "RNA", slot = "data")), center = TRUE, scale = TRUE)
  }
  
  if(!exists("eigVal"))
    eigVal <- eigen(cov(t(scaData)))$values
  
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