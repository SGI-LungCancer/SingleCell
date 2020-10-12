
#' Create TCIDEA object
#'
#' Initialize the TCIDEA object and adding option
newTCIDEA <- function(
  log.data,
  label,
  clustergroup
  ){
  
  obj <- new(Class = "TCIDEA", log.data = log.data, label = label)
  obj@ident = clustergroup

  return (obj)
}

## Add normal data ##
addNormalDataset <- function(tumor.data, tumor.ident, target.celltypes,
                             normal.data){
  
  # 1. calculate adding normal number #
  tb <- tumor.ident %>% dplyr::group_by(celltype) %>% dplyr::summarise(n = n())
  tb$percent <- tb$n / sum(tb$n) * 100
  
  ep.n <- tb$n[tb$celltype %in% target.celltypes]
  needs.normal.n = 5 * ep.n - sum(tb$n)# to adding ep percent == 20
  
  ##2. Select normal data (random)
  set.seed(1011)
  random.s <- sample(colnames(normal.data), needs.normal.n)
  
  normal.random = normal.data[,random.s]
  normal.ident = data.frame(Index = random.s, celltype = "Normal", stringsAsFactors = F)
  rownames(normal.ident) = normal.ident$Index
  
  ##3. Add normal data
  intersect.gene = intersect(rownames(tumor.data), rownames(normal.random))
  
  addnormal.data = cbind(tumor.data[intersect.gene,], normal.random[intersect.gene,])
  addnormal.cellinfo = rbind(tumor.ident, normal.ident)
  addnormal.cellinfo = addnormal.cellinfo[colnames(addnormal.data),]
  
  list <- list(data =addnormal.data, ident = addnormal.cellinfo)
  
  return (list)
}



