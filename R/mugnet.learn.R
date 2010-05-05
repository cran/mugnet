#########################################################################
# mgNetwork Class Methods 

mgSearchOrder <- function(data, perturbations = NULL,
                          maxCategories = 2, nodeCategories = 0,
                          maxParentSet = 2, maxComplexity = 0,
                          nodeOrder = NULL,
                          parentsPool = NULL, fixedParentsPool = NULL,
                          emIterations = 5, stopDelta = 0, 
                          selectMode="BIC", echo = FALSE) {

  t1 <- proc.time()

  ##if(!is.numeric(data))
  ##  stop("data shuld be numeric")
  if(sum(is.na(data))>0)
     stop("can't handle missing data yet")
    
  asframe <- FALSE
  if(is.data.frame(data)) {
    data <- as.matrix(t(data))
    asframe <- TRUE
  }
 
  numnodes <- dim(data)[1]
  numsamples <- dim(data)[2]

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
    cat("Set maxComplexity to ", maxComplexity, "\n")
  }

  if(is.null(nodeOrder))
    nodeOrder <- 1:numnodes
  
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
                    emIterations, stopDelta,
                    as.integer(selectMode), echo, 
                    PACKAGE="mugnet")

  nodenames <- seq(1,numnodes)
  if(length(dimnames(data)) == 2) {
    nodenames <- dimnames(data)[[1]]
  }
  if(length(nodenames) == numnodes && length(bestnets) > 0)
      for(i in 1:length(bestnets))
        bestnets[[i]]@nodes <- nodenames[nodeOrder]

  eval <- new("catNetworkEvaluate", numnodes, numsamples, 1)

  eval@nets <- bestnets
  for(i in 1:length(bestnets)) {
    eval@complexity[i] <- bestnets[[i]]@complexity
    eval@loglik[i] <- bestnets[[i]]@likelihood
  }

  t2 <- proc.time()
  eval@time <- as.numeric(t2[3] - t1[3])
  
  return(eval)
}

