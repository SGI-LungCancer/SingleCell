library(ggplot2)
library(pals)
library(plyr)
library(dplyr)
library(Seurat)
library(gplots)
library(RColorBrewer)

####################################################################################################
source("R/calInferredCNA_for_CEP.R") ## calculate inferredCNV value
source("R/makingTCIDEA_for_CEP.R") ## make TCIDEA object
source("R/TCIDEA_obj_for_CEP.R") ## initiate TCIDEA object
source("R/calCNVScore_for_CEP.R") ## calculate CNV score (MS, Corr) with CEP result

###################################################################################################
## PARAMETERS ##
EP_cutoff = 20 ## count if > 20% of EP -> add
target.celltypes = "EP"
label = "example"
## If user want to change the moving average parameters, please change in runCEP function ##
###################################################################################################
## Example data ##
cell_annotation_with_tumor <- readRDS(file = "example/cell_info_tumor_example.Rds") # celltype = annotation cell types in transcriptome data
tumor_example <- readRDS(file = "example/log2TPM_tumor_example.Rds")
normal_example <- readRDS(file = "example/log2TPM_normal_example.Rds")
ref_genome_example <- readRDS(file = "example/refgenome_example.Rds")
output_dir = paste0(getwd(), "/", "example")
###################################################################################################

## 1. check proportion of epithelial cells in tumor tissues ##
prop <- as.data.frame(table(cell_annotation_with_tumor$celltype))
prop$Percent = prop$Freq / nrow(cell_annotation_with_tumor) * 100

##2. Check the proportion (adding normal cells or not)
if(prop[prop$Var1 %in% target.celltypes,]$Percent > EP_cutoff){
  
  list <- addNormalDataset(tumor.data = tumor_example, tumor.ident = cell_annotation_with_tumor, target.celltypes = target.celltypes,
                           normal.data = normal_example)
  addnormal_example <- as.matrix(list$data); addnormal_annotation <- list$ident
  
  runCEP(target.normalized = addnormal_example,  
         sample.info = addnormal_annotation, label = paste0(label,"_AddNormal"),
         annotationdata = ref_genome_example, target.celltypes = target.celltypes, output.dir = output.dir) ## Sample list of EP proportion > EP_cutoff (20%)
}else{
  runCEP(target.normalized = tumor_example,  
         sample.info = cell_annotation_with_tumor, label = label,
         annotationdata = ref_genome_example,target.celltypes = target.celltypes, output.dir = output.dir) ## Sample list of EP proportion <= EP_cutoff (20%)
}

###################################################################################################


## Calculate Chromosomal expression pattern (CEP) ##
runCEP <- function(target.normalized, sample.info, label, annotationdata, min.cells = 10, MYwalk = 100,
                   target.celltypes, output.dir){

  ## 1. making TCIDEA object only tumor cells ##
  tcidea <- newTCIDEA(log.data = target.normalized, clustergroup = sample.info, label = label)

  ##2. calculate inferredCNV
  ## Proportion of epithelia cells <= EP_CUTOFF
  ## Only use genes expressed > min.cells
  ## Average of 100 genes (Binning)
  tcidea <- calInferredCNV(tcidea, min.cells = min.cells, MYwalk = MYwalk, z.score = TRUE, limit = TRUE,
                           annotationdata = annotationdata, log.file = paste0(label, "_log.txt"), use.total = FALSE)

  ## Calculate MS (Mean of squares) and CORR (Correlation)
  ## if MS score > cutoff.score or Correlation score > cutoff.corr : Malignant cells
  final_cell_info <- calCNVScore(tcidea@cnv.data, tcidea@ident, tcidea@label, levels, cutoff.score = 0.02, cutoff.corr = 0.2,meta = NULL,
                                 target.celltypes)

  ## Save calculated info
  saveRDS(final_cell_info, file = paste0(output.dir,"/", label, "_after_calc_CNV_score.Rds")) # save final_cell_info (calculated CNV score)

  ## Save Object file
  saveRDS(tcidea, file = paste0(output.dir,"/", label, "_after_calc_CNV_score_TCIDEA_obj.Rds"))

  ## Remove object
  rm(tcidea)
  rm(final_cell_info)
}

