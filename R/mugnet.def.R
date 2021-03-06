#########################################################################
# mgNetwork Class Methods 

## generate a mgNetwork object
setMethod("initialize", "mgNetwork", 
          function(.Object,...) {
            if(!validObject(.Object))
              stop("Error creating the object.")
            .Object@objectName <- "mgNetwork"
            return(.Object)
          })

setMethod("initialize", "mgNetwork", 
          function(.Object, arg1, arg2, arg3, arg4) {

            .Object@objectName <- "mgNetwork"
            
            if(nargs() >= 4 && is(arg1, "catNetwork") && (arg2=="Gaus" || arg2=="Gauss" || arg2=="Gaussian"))
              return(mgGaus(.Object, arg1, arg3, arg4))

            if(nargs() >= 3 && is(arg1, "catNetwork") && (arg2=="Pois" || arg2=="Poisson"))
              return(mgPois(.Object, arg1, arg3))

            if(nargs() >= 3 && is(arg1, "catNetwork") && (arg2=="Exp" || arg2=="Exponential"))
              return(mgExp(.Object, arg1, arg3))
            
            stop("Not appropriate initialization method is found.\n")
          })


setMethod("show", "mgNetwork",
          function(object) {
            if(is(object, "mgNetwork")) {
              cat("A ", object@model, " mgNetwork with ", object@numnodes, " nodes, ",
                  object@maxParents, " parents and ", object@maxCategories, " categories.\n Likelihood = ", object@likelihood,
                  ", Complexity = ", object@complexity, ".\n")
            }
            })

setMethod("mgSigma", "mgNetwork", 
	function(object) {
		return(object@sigmas)
	})

setMethod("mgBeta", "mgNetwork", 
	function(object) {
		return(object@betas)
	})

valid.mgNetwork <- function(obj, quietly=FALSE) {  
  res = TRUE  
  return(res)

  if(!is(obj, "mgNetwork")) {  
    res <- FALSE  
    return(res)  
  }  
  return(res)   
}  

setMethod("mgGaus", "mgNetwork", 
function(object, cnet, betas, sigmas) {
  object@model <- "Gaus"
  object@numnodes <- cnet@numnodes
  object@nodes <- cnet@nodes
  object@meta <- cnet@meta
  object@maxParents <- cnet@maxParents
  object@parents <- cnet@parents
  object@categories <- cnet@categories
  object@maxCategories <- cnet@maxCategories
  object@probabilities <- cnet@probabilities

  if(is.null(betas) || length(betas) != object@numnodes) {
    betas <- vector("list", object@numnodes)
    for(i in 1:object@numnodes) {
      ncats <- length(object@categories[[i]])
      if(ncats < 2)
        stop("A node with less than two categories: ", i)
      betas[[i]] <- seq(0, 1, 1/(ncats-1))[1:ncats]
    }
  }
  else {
    for(i in 1:object@numnodes)
      betas[[i]] <- as.numeric(betas[[i]])
  }
          
  if(is.null(sigmas)) {
    sigmas <- rep(1, object@numnodes)
    for(i in 1:object@numnodes) {
      ncats <- length(object@categories[[i]])
      if(ncats < 2)
        stop("A node with less than two categories: ", i)
      sigmas[i] <- 1/(ncats-1)
    }
  }
  
  object@betas <- betas
  object@sigmas <- sigmas

  pc <- sapply(1:object@numnodes, function(x) nodeComplexity(object, x)) 
  object@complexity <- as.integer(sum(pc))
  object@likelihood <- 0
  
  return(object)
}
)

setMethod("mgPois", "mgNetwork", 
function(object, cnet, lambdas) {
  object@model <- "Pois"
  object@numnodes <- cnet@numnodes
  object@nodes <- cnet@nodes
  object@meta <- cnet@meta
  object@maxParents <- cnet@maxParents
  object@parents <- cnet@parents
  object@categories <- cnet@categories
  object@maxCategories <- cnet@maxCategories
  object@probabilities <- cnet@probabilities

  if(is.null(lambdas) || length(lambdas) != object@numnodes) {
    lambdas <- vector("list", object@numnodes)
    for(i in 1:object@numnodes) {
      ncats <- length(object@categories[[i]])
      if(ncats < 2)
        stop("A node with less than two categories: ", i)
      lambdas[[i]] <- seq(1, ncats)
    }
  }
  else {
    for(i in 1:object@numnodes)
      lambdas[[i]] <- as.numeric(lambdas[[i]])
  }
  
  object@betas <- lambdas
  object@sigmas <- rep(0, object@numnodes)

  pc <- sapply(1:object@numnodes, function(x) nodeComplexity(object, x)) 
  object@complexity <- as.integer(sum(pc))
  object@likelihood <- 0
  
  return(object)
}
)

setMethod("mgExp", "mgNetwork", 
function(object, cnet, lambdas) {
  object@model <- "Exp"
  object@numnodes <- cnet@numnodes
  object@nodes <- cnet@nodes
  object@meta <- cnet@meta
  object@maxParents <- cnet@maxParents
  object@parents <- cnet@parents
  object@categories <- cnet@categories
  object@maxCategories <- cnet@maxCategories
  object@probabilities <- cnet@probabilities

  if(is.null(lambdas) || length(lambdas) != object@numnodes) {
    lambdas <- vector("list", object@numnodes)
    for(i in 1:object@numnodes) {
      ncats <- length(object@categories[[i]])
      if(ncats < 2)
        stop("A node with less than two categories: ", i)
      lambdas[[i]] <- seq(1, ncats)
    }
  }
  else {
    for(i in 1:object@numnodes)
      lambdas[[i]] <- as.numeric(lambdas[[i]])
  }
  
  object@betas <- lambdas
  object@sigmas <- rep(0, object@numnodes)

  pc <- sapply(1:object@numnodes, function(x) nodeComplexity(object, x)) 
  object@complexity <- as.integer(sum(pc))
  object@likelihood <- 0
  
  return(object)
}
)

mgNew <- function(nodes, cats, parents, probs = NULL, model="Gaus",  betas = NULL, sigmas = NULL) {

  cnet <- cnNew(nodes, cats, parents, probs)
  object <- new("mgNetwork", cnet, model, betas, sigmas) 
  
  if(!valid.mgNetwork(object, TRUE)) 
    stop("Incompatible parameters") 
   
  return(object)   
}

setMethod("mgReorderNodes", c("mgNetwork", "vector"),  
function(object, nodeIndices) { 
  nodeIndices <- nodeIndices[nodeIndices<=object@numnodes] 
  if(length(nodeIndices) != object@numnodes) { 
    warning("length(nodeIndices) != object@numnodes") 
    return(NULL) 
  } 
  nodeIndicesInvert <- nodeIndices 
  for(i in 1:object@numnodes) { 
    nodeIndicesInvert[i] = which(nodeIndices == i) 
  } 
 
  ##cat(nodeIndicesInvert, "\n") 
   
  newnodes <- object@nodes[nodeIndices] 
  parents <- vector("list", object@numnodes) 
  categories <- vector("list", object@numnodes) 
  probabilities <- vector("list", object@numnodes) 
  betas <- vector("list", object@numnodes)
  sigmas <- vector("list", object@numnodes)
  
  for(i in 1:object@numnodes) { 
    if(length(object@parents[[nodeIndices[i]]]) > 0) 
      parents[[i]] <- nodeIndicesInvert[object@parents[[nodeIndices[i]]]] 
    categories[[i]] <- object@categories[[nodeIndices[i]]]
    probabilities[[i]] <- object@probabilities[[nodeIndices[i]]] 
    betas[[i]] <- object@betas[[nodeIndices[i]]]
    sigmas[[i]] <- object@sigmas[[nodeIndices[i]]]
  } 
 
  object@nodes <- newnodes 
  object@categories <- categories 
  object@parents <- parents 
  object@probabilities <- probabilities 
  object@betas <- betas
  object@sigmas <- sigmas
  
  return(object) 
}) 

 
nodeComplexity <- function(object, nnode) { 
  ll <- sapply(object@parents[[nnode]], function(i) length(object@categories[[i]]))
  if(object@model == "Gaus") {
    if(length(ll)>0) 
      return((1+length(object@betas[[nnode]])) + prod(ll)*(length(object@categories[[nnode]])-1)) 
    else 
      return(2*length(object@betas[[nnode]]))
  }
  if(object@model == "Pois" || object@model == "Exp") {
    if(length(ll)>0) 
      return(length(object@betas[[nnode]]) + prod(ll)*(length(object@categories[[nnode]])-1)) 
    else 
      return(2*length(object@betas[[nnode]])-1)
  }
} 
 
setMethod("cnComplexity", signature("mgNetwork"), function(object, node) { 
  if(missing(node) || !is.numeric(node)) { 
    pc <- sapply(1:object@numnodes, function(x) nodeComplexity(object, x)) 
    return(as.integer(sum(pc))) 
  } 
  else 
    return(as.integer(nodeComplexity(object, as.integer(node))))   
  }) 
