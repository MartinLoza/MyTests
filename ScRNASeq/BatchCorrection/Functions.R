


RemoveCelltype <-  function(sim,
                            cells2remove,
                            remB1 = TRUE,
                            remB2 = TRUE, 
                            prob = 1.0){
  
  metaData <- sim$metaData
  data <- as.matrix(sim$data)
  
  celltypes <- metaData$celltype
  cellNames <- levels(celltypes)
  batch <- metaData$batch
  batchNames <- levels(batch)
  
  B1 <- data[,which(batch == batchNames[1])]
  B2 <- data[,which(batch == batchNames[2])]
  typesB1 <- celltypes[which(batch == batchNames[1])]
  batchB1 <- batch[which(batch == batchNames[1])]
  typesB2 <- celltypes[which(batch == batchNames[2])]
  batchB2 <- batch[which(batch == batchNames[2])]
  
  newB1 <- B1
  newTypeB1 <- typesB1
  newBatchB1 <- batchB1
  newB2 <- B2
  newTypeB2 <- typesB2
  newBatchB2 <- batchB2
  
  if(remB1){
    
    idxB1 <- which(typesB1[]==cells2remove)
    
    if(prob != 1.0){
      idxB1 <- sample(x = idxB1, 
                      size = floor(prob*length(idxB1)),
                      replace = FALSE)
    }
    
    newB1 <- B1[,-idxB1]
    newBatchB1 <- batchB1[-idxB1]
    newTypeB1 <- typesB1[-idxB1]
  }
  
  if(remB2){
    
    idxB2 <- which(typesB2[]==cells2remove)
    
    if(prob != 1.0){
      idxB2 <- sample(x = idxB2, 
                      size = floor(prob*length(idxB2)),
                      replace = FALSE)
    }
    
    newB2 <- B2[,-idxB2]
    newBatchB2 <- batchB2[-idxB2]
    newTypeB2 <- typesB2[-idxB2]
  }
  
  droplevels(newBatchB1)
  droplevels(newTypeB1)
  droplevels(newBatchB2)
  droplevels(newTypeB2)
  
  newCelltype <- factor(c(newTypeB1,newTypeB2))
  levels(newCelltype) = cellNames
  newBatch <- factor(c(newBatchB1,newBatchB2))
  levels(newBatch) = batchNames
  
  newMetaData <- data.frame("batch" = newBatch, "celltype" = newCelltype )
  newData <- cbind(newB1,newB2)
  rownames(newMetaData) <- colnames(newData) 
  
  return(list("data" = newData, "metaData" = newMetaData))
}

## Fold analysis and p-values
# IN: logp1 normalized data
fcAnalysis <- function(cell1 = NULL,
                       cell2 = NULL,
                       test = "limma",
                       alternative = "two.sided",
                       exact = TRUE,
                       correct = TRUE,
                       correctNegatives = TRUE){
  
  genes2Analyse <- rownames(cell1)
  
  cell1 <- expm1(cell1)
  cell2 <- expm1(cell2)
   
  if(correctNegatives){
    minCorrection <- min(min(cell1), min(cell2))
    cell1 <- cell1 - minCorrection
    cell2 <- cell2 - minCorrection
  }
  
  cell1 <- log2(cell1 + 1)
  cell2 <- log2(cell2 + 1)
  
  if(test == "limma")
    dataTest <- cbind(cell1, cell2)
  
  #Calculate p values
  pValues <- rep(NA, times = nrow(cell1))
  for(i in seq_along(pValues)){
    if(test == "t.test"){
      if(i == 1){
        print("Performing T test")
      }
      pValues[i] <- t.test(x = cell2[i,],
                           y = cell1[i,],
                           alternative = alternative)$p.value
    }else if(test == "limma"){
      if(i == 1){
        print("Performing rank sum test from Limma package")
      }
      pValues[i] <- min(1,(2*min(rankSumTestWithCorrelation(statistics = dataTest[i,],
                                                             index = seq_len(ncol(cell1))))))
    }
    else{
      if(i == 1){
        print("Performing Wilcox test")
      }
      pValues[i] <- wilcox.test(x = cell2[i,],
                                y = cell1[i,],
                                exact = exact,
                                correct = FALSE,
                                alternative = alternative)$p.value
    }
    
  }
  names(pValues) <- genes2Analyse
  
  #Calc fold change
  
  foldChange <- rowMeans(cell2) - rowMeans(cell1)
  names(foldChange) <- genes2Analyse
  
  pValues <- pValues[genes2Analyse]
  foldChange <- foldChange[genes2Analyse]
  
  results <- data.frame(Genes = genes2Analyse, 
                        FoldChanges = foldChange,
                        pValues = pValues,
                        logPValues = -log10(pValues))
  
  return(results)
}

## Plot fold changes

plotFC <- function(df, 
                   minPValue = 0.01,
                   limFoldChange = 1,
                   lineColor = "red"){
  
  if(!("logPValues" %in% colnames(df))){
    if("p_val" %in% colnames(df)){
      df$logPValues <- -log10(df$p_val)
    }else if("pValue" %in% colnames(df)){
      df$logPValues <- -log10(df$p_val)
    }else{
      print("Necessary data not found")
    }
  }
  
  if ("FoldChanges" %in% colnames(df)){
    p <- ggplot(df, aes(x = FoldChanges, y = logPValues))
  }else if("avg_logFC" %in% colnames(df)){
    p <- ggplot(df, aes(x = avg_logFC, y = logPValues))
  }else{
    print("Necessary data not found")
  }
  
  p <- p + 
    geom_point() +
    geom_hline(yintercept = -log10(minPValue),
               linetype = "dashed",
               size = 0.1,
               color = lineColor) +
    geom_vline(xintercept = -limFoldChange,
               linetype = "dashed",
               size = 0.1,
               color = lineColor ) +
    geom_vline(xintercept = limFoldChange,
               linetype = "dashed",
               size = 0.1,
               color = lineColor) 
  
  return(p)
}

## Wrapper function for ComBat

RunComBat <- function(Uncorrected, group = "batch", runningTime = FALSE){
  
  features <- VariableFeatures(Uncorrected)
  metaData <- Uncorrected[[]]
  
  time <- system.time({
    ComBat <- sva::ComBat(dat = GetAssayData(Uncorrected,
                                             assay = "RNA", 
                                             slot = "data" ), 
                          batch = metaData[,group])
  })
  
  ComBat <- CreateSeuratObject(counts = ComBat,
                               assay = "integrated",
                               meta.data = metaData)
  
  VariableFeatures(ComBat) <- features
  
  if(runningTime == FALSE)
    return(ComBat)
  else
    return(list(ComBat = ComBat, time = time))
}

## Wrapper function for ComBat-seq

RunComBatseq <- function(Uncorrected, group = "batch", runningTime = FALSE){
  features <- VariableFeatures(Uncorrected)
  metaData <- Uncorrected[[]]
  
  time <- system.time({
    ComBatSeq <- sva::ComBat_seq(counts = as.matrix(GetAssayData(Uncorrected, 
                                                                 assay = "RNA",
                                                                 slot = "counts")),
                                 batch = metaData[,group])
  })
  
  ComBatSeq <- CreateSeuratObject(counts = ComBatSeq,
                                  assay = "integrated",
                                  meta.data = metaData)
  
  VariableFeatures(ComBatSeq) <- features
  
  if(runningTime == FALSE)
    return(ComBatSeq)
  else
    return(list(ComBatSeq = ComBatSeq, time = time))
}

## Wrapper function for Scanorama

RunScanorama <- function(Uncorrected, group = "batch", runningTime = FALSE){
  
  library(here)
  scanorama <- reticulate::import('scanorama')
  
  features <- VariableFeatures(Uncorrected)
  order <- unique(Uncorrected[[]][,group])
  metaData <- Uncorrected[[]]
  xl <- SplitObject(Uncorrected, split.by = group)
  
  # Write data as table
  for(i in seq_along(xl)){
    
    data <- t(as.matrix(GetAssayData(object = xl[[order[i]]],
                                     assay = "RNA",
                                     slot = "data")[features,]))
    
    write.table(x = data, file = here(paste0("B", i, ".txt")),sep = "\t")
  }
  
  #Read data as table
  datasets <- list()
  genes_list <- list()
  for (i in seq_along(xl)) {
    datasets[[i]] <- as.matrix(read.table(here(paste0("B", i, ".txt")), sep="\t"))
    genes_list[[i]] <- features
    
    #remove temporary files
    unlink(here(paste0("B", i, ".txt")))
  }
  
  # Correction 
  time <- system.time({
    Scanorama <- scanorama$correct(datasets, genes_list, return_dense=TRUE)
  })
  
  # Create seurat object
  data <- NULL
  for (i in seq_along(xl)) {
    data <- cbind(data,t(as.matrix(Scanorama[[1]][[i]])))
  }
  
  colnames(data) <- colnames(Uncorrected)
  rownames(data) <- Scanorama[[2]]
  
  Scanorama <- CreateSeuratObject(counts = data ,
                                  meta.data = metaData,
                                  assay = "integrated")
  VariableFeatures(Scanorama) <- features
  
  if(runningTime == FALSE)
    return(Scanorama)
  else
    return(list(Scanorama = Scanorama, time = time))
}

## Wrapper function for scMerge

RunScMerge <- function(Uncorrected, group = "batch", runningTime = FALSE){
  
  library(SingleCellExperiment)
  
  features <- VariableFeatures(Uncorrected)
  metaData <- Uncorrected[[]]
  data <- as.matrix(GetAssayData(Uncorrected, 
                                 assay = "RNA",
                                 slot = "data"))[features,]
  nBatches <- length(unique(metaData[,group]))
  k <- rep(3, nBatches)
  
  scMerge <- SingleCellExperiment(assays = list(counts = data, 
                                                logcounts = data),
                                  colData = metaData)
  
  time <- system.time({
    scMerge <- scMerge::scMerge(sce_combine = scMerge,
                                ctl = features,
                                assay_name = "scMerge",
                                kmeansK = k,
                                batch_name = group,
                                plot_igraph = FALSE,
                                verbose = FALSE)
  })
  
  #Create Seurat object
  scMerge <- CreateSeuratObject(counts = as.matrix(assay(scMerge, "scMerge")),
                                assay = "integrated", 
                                meta.data = Uncorrected[[]])
  
  VariableFeatures(scMerge) <- features
  
  if(runningTime == FALSE)
    return(scMerge)
  else
    return(list(scMerge = scMerge, time = time))
}

## MNN's wrapper

RunMNN <- function(Uncorrected,
                   group = "batch",
                   runningTime = FALSE){
  
  features <- VariableFeatures(Uncorrected)
  metaData <- Uncorrected[[]]
  data <- as.matrix(GetAssayData(Uncorrected,
                                 assay = "RNA",
                                 slot = "data")[features,])
  
  time <- system.time({
    MNN <- mnnCorrect(data,
                      batch = metaData[,group],
                      cos.norm.out = FALSE)
  })
  
  #Create seurat object
  MNN <- CreateSeuratObject(counts = as.matrix(assay(MNN)),
                            meta.data = metaData,
                            assay = "integrated")
  VariableFeatures(MNN) <- VariableFeatures(Uncorrected)
  
  if(runningTime == FALSE)
    return(MNN)
  else
    return(list(MNN = MNN, time = time))
}

## Run kBet metric

RunKBET <- function(object = NULL, batch = "batch", 
                    reduction = "pca", nDim = 10, 
                    per = c(0.05, 0.1, 0.15),
                    samplePer = NULL, sampleSize = NULL){
  
  batch <- factor(object[[]][,batch])
  scores <- rep(NA, length(per))
  names(scores) <- per
  
  if(!is.null(sampleSize)){
    sample <- sample(ncol(object), size = (sampleSize), replace = FALSE)
  }else if(!is.null(samplePer)) {
    sample <- sample(ncol(object), size = (ncol(object)*samplePer), replace = FALSE)
  }else{
    sample <- seq(1,ncol(object))
  }
  
  if(reduction == "umap"){
    nDim = ncol(object[[reduction]])
  }
  
  for(p in seq_along(per)){
    k0 = floor(per[p]*(mean(table(batch)))) #neighbourhood size: mean batch size
    scores[p] <- mean(kBET(df = data.frame(Embeddings(object, reduction = reduction)[sample,1:nDim]),
                           batch = batch[sample], k0 = k0,
                           do.pca = FALSE, heuristic = FALSE,
                           plot = FALSE)$stats$kBET.observed)
  }
  return(scores)
}

RunSilhouette <- function(object = NULL, batch = "batch", 
                          reduction = "pca", nDim = 10,
                          samplePer = NULL, sampleSize = NULL){
  
  batch <- factor(object[[]][,batch])
  
  if(!is.null(sampleSize)){
    sample <- sample(ncol(object), size = (sampleSize), replace = FALSE)
  }else if(!is.null(samplePer)) {
    sample <- sample(ncol(object), size = (ncol(object)*samplePer), replace = FALSE)
  }else{
    sample <- seq(1,ncol(object))
  }
  
  if(reduction == "umap"){
    nDim = ncol(object[[reduction]])
  }
  
  score <- kBET::batch_sil(pca.data = list(x = Embeddings(object, reduction = reduction)[sample,]),
                           batch = batch, nPCs = nDim)
  
  return(score)
}

# Wrapper to use RISC on Seurat objects
RunRISC <- function(Uncorrected = NULL, batch = NULL){
  
  if(!is.list(Uncorrected)){
    yl <- SplitObject(Uncorrected, split.by = batch)
  }else{
    yl <- Uncorrected
  }
  
  yl <- lapply(X = yl, FUN = GetAssayData, assay = "RNA", slot = "counts")
  
  # Pre-processing steps using RISC
  yl <- lapply(yl, FUN = function(i){
    tmp <- as.matrix(i)
    genes <- data.frame(rownames(tmp))
    rownames(genes) <- rownames(tmp)
    cell <- data.frame(colnames(tmp))
    rownames(cell) <- colnames(tmp)
    return(RISC::readscdata(count = tmp, cell = cell, gene = genes))
  })
  
  yl <- lapply(yl, scFilter)
  yl <- lapply(yl, scNormalize)
  yl <- lapply(yl, scDisperse)
  
  y <- RISC::scMultiIntegrate(objects = yl)
  rm(yl)
  
  # Re-arrange data to use Seurat objects
  colData <- y@coldata
  rowData <- y@rowdata
  tmp <- as.matrix(y@assay$logcount)
  colnames(tmp) <- colData$scBarcode
  rownames(tmp) <- rowData$Symbol
  rm(y)
  
  if(is.list(Uncorrected)){
    x <- Reduce(merge, Uncorrected)
  }else{
    x <- Uncorrected
  }
  
  x <- x[rowData$Symbol, colData$scBarcode] 
  x[["integrated"]] <- Seurat::CreateAssayObject(counts = tmp)
  DefaultAssay(x) <- "integrated"
  
  return(x)
}

























