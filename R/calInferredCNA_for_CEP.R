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
                                 target.celltypes, output.dir)
  
  ## Save calculated info
  saveRDS(final_cell_info, file = paste0(output.dir,"/", label, "_after_calc_CNV_score.Rds")) # save final_cell_info (calculated CNV score)
  
  ## Save Object file
  saveRDS(tcidea, file = paste0(output.dir,"/", label, "_after_calc_CNV_score_TCIDEA_obj.Rds"))
  
  ## Remove object
  rm(tcidea)
  rm(final_cell_info)
}



calInferredCNV <- function(
  obj, min.cells = 10, MYwalk = 100, z.score = TRUE,  limit = TRUE, 
  annotationdata, log.file = "log_files.txt",
  use.total = FALSE
  ){

  ##
  GTF_uniq <- annotationdata

  ##
  SC = obj@log.data

  cat("Calculate inferred CNV in ", obj@label, "class\n")

  ### 1. Filter out unreliable genes
  SC. = as.matrix(SC[rowMeans(SC) != 0,])
  dim <-  dim(SC.)
  cat("Raw data : ", dim[1], "genes across", dim[2], "samples\n")

  # Filter-out many genes (expr>1, min.cells>10)
  SC.filter = SC.[rowSums(SC.>1) > min.cells,]; dim(SC.filter)
  dim <-  dim(SC.filter)
  cat("QC passed data : ", dim[1], "genes across", dim[2], "samples\n")

  ### 2. rm low-expressed genes
  if(use.total){
    SC.f <- SC.filter
  }else{ SC.f <- SC.filter[rowMeans(SC.filter) > 0.1,]}

  dim <- dim(SC.f)
  cat("2nd QC passed data (rm low-expressed genes) : ", dim[1], "genes across", dim[2], "samples\n")

  SC.o.substract <- SC.f
  ### 4. "Annotation" of gene symbol with chromosomal information
  SC_anno = QAnno(SC.o.substract, GTF_uniq); dim(SC_anno)

  ### 5. Omit X & Y & MT chromosome
  Omit_XYM = c(23,24,25,26) # c("X","Y","M","GL"|"KI")
  SC_auto = SC_anno[!grepl(paste(Omit_XYM,collapse="|") , SC_anno$chromosome_name), ]
  dim(SC_auto)

  dim <-  dim(SC_auto)
  cat("X & Y & MT chromosome removed data : ", dim[1], "genes across", dim[2]-6, "samples\n")

  NormalZ_auto.r <- NULL
  if(z.score == TRUE){
    ### 5. "Z-scoring" by row
    SC_SD = data.matrix(apply(SC_auto[,-c(1:6)], 1, sd))
    SCZ_auto.r = cbind(SC_auto[,1:6], Zscore(SC_auto[,-c(1:6)], SC_SD)); typeof(SCZ_auto.r) # "list"
    cat("Making Z-scoring data\n")
  }
  else{
    SCZ_auto.r = SC_auto
    cat("Making Not centering, Z-scoring data\n")
  }

  ### 6. Limit the relative expression values to [-3,3] # as Tirosh did
  summary(as.numeric(as.matrix(SCZ_auto.r[,-c(1:6)])))
  if(limit){
    SCZ_auto.r2 = SCZ_auto.r[,-c(1:6)]
    SCZ_auto.r2[SCZ_auto.r2 < (-3)] <- (-3)
    SCZ_auto.r2[SCZ_auto.r2 > 3] <- 3
    summary(as.numeric(as.matrix(SCZ_auto.r2)))
    SCZ_auto.r2 = cbind(SCZ_auto.r[,1:6], SCZ_auto.r2)
    cat("Limit the relative expression values : limitation is -3 ~ 3\n")
    #    write.table(SCZ_auto.r2, file = "LUNG_T18_zscore_by_row_after_lim.txt", sep= "\t")
    cat("After limitation, min : ",min(SCZ_auto.r2[,-c(1:6)])," max : ", max(SCZ_auto.r2[,-c(1:6)]), "\n")
  }
  else{
    SCZ_auto.r2 = SCZ_auto.r
    cat("Not Limit the relative expression values\n")
    cat("min : ",min(SCZ_auto.r2[,-c(1:6)])," max : ", max(SCZ_auto.r2[,-c(1:6)]), "\n")
  }


  ### 7. Moving average of Z-score / Centering data
  library(caTools) ;

  SCZ_MV150 = MyMV_Zscore(SCZ_auto.r2, MYwalk)
  SCZ_MV150_centering = t(t(SCZ_MV150)-colMeans(SCZ_MV150)); typeof(SCZ_MV150_centering) # "double"

  f_SCZ_MV150 = round(SCZ_MV150, digits=3);
  f_SCZ_MV150_centering = round(SCZ_MV150_centering, digits=3);  ## centering by column

  ## return cnv values. z-scored values
  obj@cnv.data <- SCZ_MV150_centering
  obj@cnv.sd <- calSD(SCZ_MV150_centering)

  ## add statistical values
  obj@cnv.mean.sq <-  meanSquare(SCZ_MV150_centering)
  obj@cnv.cv <- calCV(SCZ_MV150_centering)
  obj@cnv.abmean <- calMeanAb(SCZ_MV150_centering)

  return (obj)
}

## calculate mean of absolute ##
calMeanAb <- function(matrix){
  mean.ab = matrix(nrow = ncol(matrix), ncol = 1)
  rownames(mean.ab) <- colnames(matrix)

  for(i in 1:ncol(matrix)){
    abs <- abs(x = matrix[,i])
    mean <- mean(x=abs)
    #
    mean.ab[i,1] <- mean
  }
  ##
  colnames(mean.ab) <- "abMean"
  cat(colnames(mean.ab))

  return (mean.ab)
}

## calculate CV ##
calCV <- function(matrix){
  cv = matrix(nrow = ncol(matrix), ncol = 1)
  rownames(cv) <- colnames(matrix)

  for(i in 1:ncol(matrix)){
    sample.SD <- sd(x = matrix[,i])
    abs <- abs(x = matrix[,i])
    mean <- mean(x=abs)
    cv[i,1] <- sample.SD / mean * 100 #150518
  }

  colnames(cv) <- "CV"
  cat(colnames(cv))
  return (cv)
}

## mean of squares : mean(squares per each values)
meanSquare <- function(matrix){

  m.square = matrix(nrow = ncol(matrix), ncol = 1)
  rownames(m.square) <- colnames(matrix)

  for(i in 1:ncol(matrix)){
    square <- 0; n <- 0;
    for(j in 1:nrow(matrix)){
      square <- square + (matrix[j,i])^2
      n <- n+1
    }
    ##
    m.square[i,1] <- square / n
  }

  colnames(m.square) <- "MS"
  cat(colnames(m.square))

  return(m.square)
}


## in this case, we only provide GRCh38 reference genome.
readGTF.addData <- function(AnnotationLevel){

  if(AnnotationLevel == "GRCh38"){load("GRCh38.rda"); return(GTF_uniq)}
  else{
    cat('TCIDEA only provides GRCh38 gtf files.')
  }
}

## using GTF matrix, we make chromosome table along chromosomal location ##
QAnno = function(row_gene_named_matrix, GTF_uniq){
  ANNO_overlap = intersect(rownames(row_gene_named_matrix) , rownames(GTF_uniq))
  ANNO_GTF     = GTF_uniq[ANNO_overlap, ]
  ANNO_input   = row_gene_named_matrix[ANNO_overlap,];
  ANNO_merge   = cbind(ANNO_GTF,ANNO_input)
  MyANNO       = cbind(ANNO_merge[,1],rownames(ANNO_merge),ANNO_merge[,c(2:5 , 7:ncol(ANNO_merge))])
  colnames(MyANNO)=c("ensembl_gene_id","gene_name","description","chromosome_name","start_position","end_position",colnames(ANNO_input))
  MyANNO[,4]=gsub("X",23,MyANNO[,4]) ; MyANNO[,4]=gsub("Y",24,MyANNO[,4]) ; MyANNO[,4]=gsub("MT",25,MyANNO[,4]) ; MyANNO[,4]=gsub("GL",26,MyANNO[,4]) ; MyANNO[,4]=gsub("KI",26,MyANNO[,4])
  MyANNO[,c(4,5,6)] = sapply(MyANNO[,c(4,5,6)], as.numeric)
  MyTABLE = MyANNO[order(MyANNO[,4] , MyANNO[,5]) , ]
  return(MyTABLE)
}

## calculate z-score ##
Zscore = function(Tumor_ExpRatio, SD_of_TumorRatios){(Tumor_ExpRatio-rowMeans(Tumor_ExpRatio))/SD_of_TumorRatios}

## Normalize CEP to Z-score ##
MyMV_Zscore = function(annotate_matrix, MYwalk){
  annotate_matrix$chromosome_name = sapply(annotate_matrix$chromosome_name, as.numeric)
  MV_input = annotate_matrix[order(annotate_matrix$chromosome_name), ]
  MV_input = MV_input[,-c(1:3,5:6)] ;
  rownames(MV_input)=paste("chr",annotate_matrix[,4],":",annotate_matrix[,5],"-",annotate_matrix[,6]," (",annotate_matrix[,2],")",sep="")
  for(i in 1:22){
    MV.chr = MV_input[MV_input[,1] == i, ];
    MV.dat = apply(MV.chr, 2, runmean, MYwalk)
    if(i ==1){MyMV = MV.dat} else {MyMV=rbind(MyMV,MV.dat)}}
  MV_output = MyMV[,-1];
  colnames(MV_output)=colnames(MV_input[,-1]) ; rownames(MV_output) = rownames(MV_input[,-1]);
  return(MV_output)
}

# calculate SD of each single-cells
calSD <- function(matrix){
  SD = matrix(nrow = ncol(matrix), ncol = 1)
  rownames(SD) <- colnames(matrix)

  for(i in 1:ncol(matrix)){
    sample.SD <- sd(x = matrix[,i])
    SD[i,1] <- sample.SD
  }

  colnames(SD) <- "SD"
  cat(colnames(SD))
  return (SD)
}
