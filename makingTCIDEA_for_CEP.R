
#' Create TCIDEA object
#'
#' Initialize the TCIDEA object and adding option
#' @param raw.data Raw input data.
#' @param logNormalized raw data was already logNormalized or not
#' @param clustergroup Using clustering method or not
#' @param saveRawData save raw.data or not
newTCIDEA <- function(
  raw.data,
  logNormalized = FALSE,
  clustergroup = NULL,
  saveRawData = FALSE,
  normal.use = TRUE,
  normal.raw.data,
  normal.logNormalized = FALSE,
  tSNE.data = NULL,
  tSNE.calculate = FALSE,
  label,
  save.loc = NULL
  ){
  obj <- new(Class = "TCIDEA", raw.data = raw.data, label = label)

  ##check logNormalized
  if(logNormalized == FALSE){
    obj@log.data = obj@raw.data
  }
  else{
    library(Seurat)
    obj@log.data = LogNormalize(obj@raw.data, scale.factor = 1e4, display.progress = TRUE)
  }

  PCs = 20
  ##cluster group exist or not
  if(!is.null(clustergroup)){
    clustergroup = checkClusterGroupOrder(colnames(obj@log.data), clustergroup)
    obj@ident = clustergroup
  }else{
    Cluster.list <- runningSeuratForCluster(obj@log.data, label, PCs, save.loc)   #PCs : 20
    if(tSNE.calculate == TRUE){obj@tSNE.data = Cluster.list$tsne; obj@ident = Cluster.list$cluster}
  }

  ##use normal data in inferredCNV
  obj@normal.use = normal.use
  if(normal.use){
    obj@normal.raw.data = normal.raw.data
    ## if normal.raw.data -> NULL : use default dataset
    if(!normal.logNormalized){obj@normal.log.data = normal.raw.data}
    else{
      library(Seurat)
      obj@normal.log.data = LogNormalize(obj@normal.raw.data, scale.factor = le4, display.progress = TRUE)}
  }

  ##
  if(!is.null(tSNE.data)){
    if(tSNE.calculate == FALSE){
      obj@tSNE.data = tSNE.data
    }
    else{
      ##calculate tSNE data plz.
      obj@tSNE.calculate == FALSE
    }
  }

  ##save raw data or not
  if(saveRawData == TRUE){
  }
  else{
    obj@raw.data = matrix()
    obj@normal.raw.data = matrix()
  }

  return (obj)
}


checkClusterGroupOrder <- function(
  original.names,
  reordered.names
  ){
  original.names <- data.frame(barcode = original.names)
  reordered.names <- data.frame(barcode = rownames(reordered.names), cluster = reordered.names[,1])

  if(nrow(original.names) != nrow(reordered.names)){cat("Please check cluster data. Different number of rows\n"); stop()}
  else{
    library(plyr)
    original.names.annot <- join(original.names, reordered.names, by = "barcode")

    return.names <- data.frame(cluster = original.names.annot$cluster)
    rownames(return.names) <- original.names.annot$barcode

    return (return.names)
  }
}

runningSeuratForCluster <- function(
  Data, name, pc, save.loc
){
  library(Seurat)
  ## source code by Nayoung Kim
  ## Cutoff for cells and genes
  min_cells = ceiling(0.001 * dim(Data)[2]) # genes expressed in >= (~0.1% of the data) cells
  min_genes = 0

  ## data input
  sr1 = CreateSeuratObject(raw.data=Data, min.cells=min_cells, min.genes=min_genes, is.expr=0, names.field=1, names.delim="_", project=name)
  sr1

  ## Log Normalized data (TPM-like values)
  sr1 <- NormalizeData(object=sr1, normalization.method="LogNormalize", scale.factor=10000)
  sr1@data <- as.matrix(sr1@raw.data) # change the values into original log2 (TPM+1)

  ## Centering the data (centering)
  sr1 <- ScaleData(object = sr1, do.scale = FALSE, do.center = TRUE, scale.max = 10)

  ## Run RCA on all genes using IRLBA : use total genes
  sr1 <- RunPCA(sr1, pc.genes=rownames(sr1@data), do.print=FALSE) # PCA on all genes

  ## Combiend PCA plot
  pc_Data <- data.frame (sr1@dr$pca@cell.embeddings[, 1:5])

  if(!is.null(save.loc)){OutPutFile_PCA = paste0(save.loc, "/", "PCA_plot_", name, ".pdf")}
  else{OutPutFile_PCA = paste0("PCA_plot_", name, ".pdf")}

  pdf(file=OutPutFile_PCA, width=7, height=6.5)
  plot(pc_Data, pch=16, col=rgb(0,0,0,0.5), cex=0.5)
  dev.off()
  cat("Save PCA plot\n")

  ## PC heatmap
  if(!is.null(save.loc)){OutPutFile_PCmap = paste0(save.loc, "/","PCA_heatmap_", name, ".pdf")}
  else{OutPutFile_PCmap = paste0("PCA_heatmap_", name, ".pdf")}
  pdf(file=OutPutFile_PCmap, width=14, height=12)
  PCHeatmap(object = sr1, pc.use = 1:12, num.genes = 30, do.balanced = TRUE, use.scale = TRUE, label.columns = FALSE)
  dev.off()
  cat("Save PCA heatmap\n")

  ## Determine statistically significant principal components
  if(!is.null(save.loc)){OutPutFile_PCElbow = paste0(save.loc, "/","PCElbow_plot_", name, ".pdf")}
  else{OutPutFile_PCElbow = paste0("PCElbow_plot_", name, ".pdf")}
  pdf(file=OutPutFile_PCElbow, width=3.7, height=3.2)
  PCElbowPlot(sr1)
  dev.off()
  cat("Save PCElbow plot\n")


  ## Cluster the cells
  dims = 1:pc
  sr1 <- FindClusters(sr1, reduction.type = "pca", dims.use = dims, resolution = 1.2, print.output = 0, save.SNN = TRUE)

  ## tSNE
  sr1 = RunTSNE(sr1, dims.use=dims, do.fast=TRUE)

  ident = as.matrix(sr1@ident)
  ident.df = data.frame(cluster = ident)
  rownames(ident.df) = rownames(ident)
  list = list(cluster = ident.df, tsne = data.frame (sr1@dr$tsne@cell.embeddings))

  return(list)
}


