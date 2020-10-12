
tcidea <- setClass(
  "TCIDEA",
  slots = c(
    log.data = "ANY",
    scale.data = "ANY",
    ident = "ANY",
    tSNE.data = "ANY",
    tSNE.calculate = "ANY",
    cnv.data = "ANY",
    cnv.sd = "ANY",
    cnv.mean.sq = "ANY",
    cnv.cv = "ANY",
    cnv.abmean = "ANY",
    cluster.auc = "ANY",

    ##
    normal.use = "ANY",
    normal.raw.data = "ANY",
    normal.log.data = "ANY",

    ##edge index per cluster and gene set name -> T/F
    cluster.es.esp = "ANY",
    cluster.geneset = "ANY",
#    cluster.rank = "ANY",

    ##
    label = "ANY"
  )
)


##Documentation
#raw.data = raw umi count data
#log.data = log data -> log2 normalized
#scale.data = scaled data -> centering log2 normalized data
setMethod(
  f = "show",
  signature = "TCIDEA",
  definition = function(object) {
    cat(
      "An object of class",
      object@label,
      "\n",
      nrow(x = object@log.data),
      "genes across",
      ncol(x = object@log.data),
      "samples.\n"
    )
    if(object@normal.use){
      cat(
        "Normal reference dataset ",
        nrow(x = object@normal.log.data),
        "genes across",
        ncol(x = object@normal.log.data),
        "samples.\n"
      )
    }
    invisible(x = NULL)
  }
)


