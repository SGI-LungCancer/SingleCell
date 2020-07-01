library(ggplot2)
library(pals)
library(plyr)
library(dplyr)
library(Seurat)
library(gplots)
library(RColorBrewer)

##################################
source("R/calInferredCNA.R") ## calculate inferredCNV value
source("R/makingTCIDEA.R") ## make TCIDEA object
source("R/TCIDEA_obj_for_CEP.R") ## initiate TCIDEA object
source("R/calCNVScore_for_CEP.R") ## calculate CNV score (MS, Corr) with CEP result
#################################

## write your sample names ##
tumor <- c("Tumor_Example")
normal <- c("Normal_Example")

final_cell_all.merged <- readRDS("cell_annotation_tumor_and_normal.Rds")
tumor_example <- readRDS("tumor_example_log2_normalized.Rds")
normal_example <- readRDS("normal_example_log2_normalized.Rds")

load(file = "Normal_Seurat.rds") # Pan-normal log2-normalized(TPM+1) data
load(file = "Normal_info_data.rds") # Cell type information of Pan-normla data
load("final_cell_all_merged.rds") # All of Tumor log2-normalized(TPM+1) data
load(

## count if > 20% of EP -> add
EP_cutoff = 20

for(s in samples){
  sample.info <- subset(final_cell_all.merged, Samples == s)
  #
  tb <- sample.info %>% dplyr::group_by(cell_type) %>% dplyr::summarise(n = n())
  tb$percent <- tb$n / sum(tb$n) * 100

  setwd("") # set working directory
  dir.create(s, showWarnings = TRUE)
  if(tb$percent[tb$cell_type == "EP"] > EP_cutoff){
    runCEP.AddNormal(s, sample.info)
	## Sample list of EP proportion > EP_cutoff (20%)

  }
  runCEP(s, sample.info)
  ## Sample list of EP proportion <= EP_cutoff (20%)
}


## Calculate Chromosomal expression pattern (CEP)
runCEP <- function(sample, sample.info){

  setwd("")
  dir.create(paste0(sample,"/Original"), showWarnings = TRUE)
  setwd("/Original")

  ##1. load each dataset
  log2tpm.name = paste0("log2TPM_QCpass_",sample,"_rmRedu_UMI.txt")
  data <- read.table(file = log2tpm.name, sep = "\t", row.names = 1, header = T)
  ident <- data.frame(celltype = sample.info$cell_type, row.names = sample.info$Barcode)
  head(ident)
  label = paste0(sample, "_original")

  ##2. making TCIDEA obj
  tcidea <- newTCIDEA(raw.data = data, logNormalized = FALSE, clustergroup = ident, saveRawData = FALSE,
                      normal.use = F, normal.raw.data = NULL, normal.logNormalized = F,
                      label = label, tSNE.data = NULL, tSNE.calculate = FALSE)
  head(tcidea@ident)

  ##3, calculate inferredCNV
  ## Proportion of epithelia cells <= EP_CUTOFF
  ## Only use genes expressed > min.cells
  ## Average of 100 genes (Binning)
  memory.limit(500000)
  tcidea <- calInferredCNV.NotNormal(tcidea, min.cells = 10, MYwalk = 100, z.score = TRUE, limit = TRUE,
                                     AnnotationLevel = "GRCh38", log.file = paste0(label, "_log.txt"), use.total = FALSE)

  levels = as.character(unique(tcidea@ident$cluster))
#  levels

  ## Calculate MS (Mean of squares) and CORR (Correlation)
  ## if MS score > cutoff.score or Correlation score > cutoff.corr : Malignant cells
  final_cell_info <- calCNVScore(tcidea@cnv.data, tcidea@ident, tcidea@label, levels, cutoff.score = 0.02, cutoff.corr = 0.2,meta = NULL)

  ## Save calculated info
  save(final_cell_info, file = "") # save final_cell_info (calculated CNV score)

  ## Save Object file
  save(tcidea, file = "")

  ## Remove object
  rm(tcidea)
  rm(final_cell_info)
}

## Calculate Chromosomal expression pattern (CEP) after Adding normal data ##
runCEP.AddNormal <- function(sample, sample.info){

  setwd("")
  dir.create(paste0(sample,"/Add_normal"), showWarnings = TRUE)
  setwd(paste0("", sample, "/Add_normal"))

  ##1. load each dataset
  log2tpm.name = paste0(sample, "/log2TPM_QCpass_",sample,"_rmRedu_UMI.txt")
  data <- read.table(file = log2tpm.name, sep = "\t", row.names = 1, header = T)
  ident <- data.frame(celltype = sample.info$cell_type, row.names = sample.info$Barcode)
  head(ident)
  label = paste0(sample, "_Add_Normal")

  ##1-1. add Normal data
  list <- addNormalDataset(data, ident)
  data <- as.matrix(list$data); ident <- list$ident

  cat(dim(data))

  ##2. making TCIDEA obj
  tcidea <- newTCIDEA(raw.data = data, logNormalized = FALSE, clustergroup = ident, saveRawData = FALSE,
                      normal.use = F, normal.raw.data = NULL, normal.logNormalized = F,
                      label = label, tSNE.data = NULL, tSNE.calculate = FALSE)
  head(tcidea@ident)

  identical(rownames(ident), colnames(data))

  ##3, calculate inferredCNV
  ## Proportion of epithelia cells <= EP_CUTOFF
  ## Only use genes expressed > min.cells
  ## Average of 100 genes (Binning)
  memory.limit(500000)
  tcidea <- calInferredCNV.NotNormal(tcidea, min.cells = 10, MYwalk = 100, z.score = TRUE, limit = TRUE,
                                     AnnotationLevel = "GRCh38", log.file = paste0(label, "_log.txt"), use.total = FALSE)

  levels = as.character(unique(tcidea@ident$cluster))
  #  levels

  ## Calculate MS (Mean of squares) and CORR (Correlation)
  ## if MS score > cutoff.score or Correlation score > cutoff.corr : Malignant cells
  final_cell_info <- calCNVScore(tcidea@cnv.data, tcidea@ident, tcidea@label, levels, cutoff.score = 0.02, cutoff.corr = 0.2,meta = NULL)

  ## Save calculated info
  save(final_cell_info, file = "") # save final_cell_info (calculated CNV score)

  ## Save Object file
  save(tcidea, file = "")

  ## Remove object
  rm(tcidea)
  rm(final_cell_info)
}


## Add normal data ##
## if % of EP > 20, add normal data and change the proportion of EP.
## using Seurat 2.3.4 function
addNormalDataset <- function(tumor.data, tumor.ident){

  # 1. calculate adding normal number
  tb <- tumor.ident %>% dplyr::group_by(celltype) %>% dplyr::summarise(n = n())
  tb$percent <- tb$n / sum(tb$n) * 100

  ep.n <- tb$n[tb$celltype == "EP"]
  needs.normal.n = 5 * ep.n - sum(tb$n)# to adding ep percent == 20

  ##2. add normal dataset
  set.seed(1011)
  random.s <- sample(rownames(normal_info.data), needs.normal.n)
  normal.use <- SubsetData(normal.seurat, cells.use = random.s, do.center = T, do.scale = F)
  normal.use
  normal_info.data$tmp <- rownames(normal_info.data)
  normal_info.data$celltype <- as.character(normal_info.data$celltype)
  str(normal_info.data)
  normal.use.cluster <-  subset(normal_info.data, normal_info.data$tmp %in% random.s)
  normal.use.cluster$tmp <- NULL
  str(normal.use.cluster)
  #  normal.use.cluster$celltype <- as.character(normal.use.cluster$celltype)

  cat("Adding normal data : ", nrow(normal.use.cluster))

  ## calculate normal ratio
  normal.use.cluster$celltype <- paste0(normal.use.cluster$celltype, "_N")
  merged.ident <- rbind(tumor.ident, normal.use.cluster) # tumor - normal rbind
  merged.ident$celltype <- as.character(merged.ident$celltype)

  ## create seurat obj
  tumor.seurat <- CreateSeuratObject(raw.data = tumor.data, min.cells = 0, min.genes = 0,
                                     is.expr=0, names.field=1, names.delim="_", project="Tumor")
  tumor.seurat <- NormalizeData(object=tumor.seurat, normalization.method="LogNormalize", scale.factor=10000)
  tumor.seurat <- ScaleData(object = tumor.seurat, do.scale = FALSE, do.center = TRUE, scale.max = 10)

  ## Merge Tumor and Normal data ##
  merged.seurat <- MergeSeurat(tumor.seurat, normal.use, min.cells = 0, min.genes = 0, do.scale = F, do.center = T)
  merged.seurat

  list <- list(data = merged.seurat@data, ident = merged.ident)

  return (list)
}


