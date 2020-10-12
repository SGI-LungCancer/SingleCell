calCNVScore <- function(sh, cell_info, s, levels, cutoff.corr, cutoff.score, meta, target.celltypes){

  sh2 = sh; dat = t(sh2) ## already remove chromosomal location info
  CNV_score <- data.frame(MS = colMeans(sh2^2), SS = apply(sh2^2, 2, sum), SD = apply(sh2, 2, sd))
  ##
  CNV_score$Row.names <- rownames(CNV_score)
  cell_info$Row.names <- rownames(cell_info)

  cell_info2 <- plyr::join(cell_info, CNV_score, by="Row.names") # boxplot for celltype

  ## MS top 5% cells
  top_MS_cells <- arrange(cell_info2, desc(MS))[1:round(dim(cell_info2)[1]*0.05),]$Row.names  # Top 5%

  ## calculate correlation : corr using 1 cell vs. top_MS_cells
  tmp <- data.frame(Ave_tumor = rowMeans(sh2[,top_MS_cells]))
  for(i in 1:dim(cell_info2)[1]){
    cell_info2$COR[i] <-  cor(sh2[,i, drop=FALSE], data.frame(Ave_tumor = rowMeans(sh2[,top_MS_cells])))
  }


  ##
  cell_info3 <- cell_info2
  rownames(cell_info3) <- cell_info3$Row.names

  tumorcells <- filter(cell_info3, ((MS > cutoff.score | COR > cutoff.corr) & celltype %in% target.celltypes))$Row.names
  nontumorcells <- cell_info3$Row.names[!(cell_info3$Row.names %in% c(tumorcells))]
  immunecells <- cell_info3$Row.names[cell_info3$celltype != target.celltypes]

  ## only classified tumor vs. non-tumor ##
  cell_info3$cell_index <- rep("X", dim(cell_info3)[1])
  cell_info3[tumorcells,]$cell_index <- "Tumor"
  cell_info3[nontumorcells,]$cell_index <- "Nontumor"
  cell_info3[immunecells,]$cell_index <- "Immune"

  ## 2D plot of MS score and correlation ##
  expos<-ggplot(cell_info3, aes(x=MS, y= COR)) + geom_point(aes(fill=cell_index), size=5, alpha=.8, shape=21, colour="black") +
    scale_fill_manual(values = c("Tumor"="red","Immune" = "gray70","Nontumor"="dodgerblue1")) +
    geom_vline(xintercept = cutoff.score, colour="black", size=0.5, linetype = "longdash") + geom_hline(yintercept = cutoff.corr, colour="black", size=0.5, linetype = "longdash") +
    xlab("MS score") + ylab("CNV correlation") + theme_bw() +
    theme(axis.title.x = element_text(face="bold", size=16), axis.text.x  = element_text(size=12)) +
    theme(axis.title.y = element_text(face="bold", size=16), axis.text.y  = element_text(size=12)) +
    theme(panel.border=element_rect(fill=NA, colour="black", size=2), legend.position = 'right')
  expos
  ggsave(paste0("final_",s, "_all_cells_CNV_score_vs_cor_classification.pdf"), width = 7, height = 5)

  return (cell_info3)
}

ReadTotalCelltype <- function(FinalCellType){

  total.celltype <- read.table(file = FinalCellType,
                               sep = "\t", header = T)

  total.celltype$Cell <- as.character(total.celltype$Cell)
  total.celltype$NEW <- as.character(total.celltype$NEW)

  for(i in 1:nrow(total.celltype)){
    tmp <- base::strsplit(x = total.celltype$Cell[i], split = "_")
    total.celltype$Sample[i] <- tmp[[1]][1]; total.celltype$Barcode[i] <- tmp[[1]][2]
  }

  return (total.celltype)
}

ReadClusterForSample <- function(s, total){

  if(grep("-A1", s) > 0){s <- gsub("-A1", "", s)}

  ##
  sample.subset <- subset(total, total$Sample == s)
  dim <- dim(sample.subset)
  cat("Sample ",s," Cell : ", dim[1],"\n")

  ##
  df <- data.frame(cluster = sample.subset$NEW)
  rownames(df) <- sample.subset$Barcode

  return (df)
}

