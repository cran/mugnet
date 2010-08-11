#########################################################################
# mgNetwork Class Methods 

mgSearchOrder <- function(data, perturbations = NULL,
                          maxCategories = 2, nodeCategories = 0,
                          maxParentSet = 2, maxComplexity = 0,
                          nodeOrder = NULL,
                          parentsPool = NULL, fixedParentsPool = NULL,
                          emIter = 5, stopDelta = 0, emStartIter = 2,
                          selectMode="BIC", model="Gaus",
                          echo = FALSE) {

  t1 <- proc.time()

  ##if(!is.numeric(data))
  ##  stop("data shuld be numeric")
  if(sum(is.na(data))>0)
     stop("can't handle missing data yet")
    
  asframe <- FALSE
  if(is.data.frame(data)) {
    data <- as.matrix(t(data))
    if(!is.null(perturbations))
      perturbations <- as.matrix(t(perturbations))
    asframe <- TRUE
  }
 
  numnodes <- dim(data)[1]
  numsamples <- dim(data)[2]

  for(i in 1:numnodes) {
    if(!is.numeric(data[,i]))
      stop("Data should be numeric")
  }
  
  if(!is.null(nodeCategories)) {
    if(length(nodeCategories) != numnodes)
      stop("Specify the number of categories for all nodes")
    maxCategories <- 0
    for(cc in nodeCategories)
      if(maxCategories < cc)
        maxCategories <- cc
  }
  if(maxCategories < 2) {
    warning("Set maxCategories = 2")
    maxCategories <- 2
  }
  if(is.null(nodeCategories))
    nodeCategories <- rep(maxCategories, numnodes)

  if(maxParentSet < 1) {
    warning("Set maxParentSet = 1")
    maxParentSet <- 1
  }

  for(i in 1:numnodes)
    data[i,] <- as.numeric(data[i,])

  if(maxComplexity <= 0) {
    maxComplexity <- as.integer(numnodes * (exp(log(maxCategories)*maxParentSet) * (maxCategories-1) + (maxCategories+1)))
    if(echo)
      cat("Set maxComplexity to ", maxComplexity, "\n")
  }

  if(is.null(nodeOrder))
    nodeOrder <- 1:numnodes
  else {
    if(length(nodeOrder) != numnodes) {
      nodeOrder <- 1:numnodes
      warning("nodeOrder set to ", nodeOrder)
    }
  }
  
  if(numnodes < 1 || numsamples < 1)
    stop("No valid data is specified.")

  if(selectMode == "AIC")
    selectMode <- 1
  else {
    if(selectMode == "BIC")
      selectMode <- 2
    else
      selectMode <- 3
  }

  bestnets <- .Call("mgSearchOrderC",
                    data, perturbations,
                    nodeCategories, maxParentSet,
                    maxComplexity, nodeOrder,
                    parentsPool, fixedParentsPool,
                    emIter, stopDelta, emStartIter, 
                    as.integer(selectMode),
                    as.character(model), 
                    echo, 
                    PACKAGE="mugnet")

  nodenames <- seq(1,numnodes)
  if(length(dimnames(data)) == 2) {
    nodenames <- dimnames(data)[[1]]
  }
  if(length(nodenames) == numnodes && length(bestnets) > 0)
      for(i in 1:length(bestnets))
        bestnets[[i]]@nodes <- nodenames[nodeOrder]

  eval <- new("catNetworkEvaluate", numnodes, numsamples, 1)

  for(i in 1:length(bestnets)) {
    eval@complexity[i] <- bestnets[[i]]@complexity
    eval@loglik[i] <- bestnets[[i]]@likelihood

    ## reorder eval@nets[[nn]]'s nodes to match data's nodes
    enetnodes <- bestnets[[i]]@nodes
    if(length(nodenames) == numnodes) {
      ord <- sapply(nodenames, function(c) {
        id <- which(enetnodes==c)
        if(length(id)>0)
          return(id[1])
	stop("nodes do not match")
        })
      if(sum(ord != 1:numnodes) > 0) {
        bestnets[[i]] <- mgReorderNodes(bestnets[[i]], ord)
      }
    }
  }
  eval@nets <- bestnets

  t2 <- proc.time()
  eval@time <- as.numeric(t2[3] - t1[3])
  
  return(eval)
}

