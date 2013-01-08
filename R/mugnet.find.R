#########################################################################
# Categorical Network Class Methods
# Find-members

setMethod("mgFind", "catNetworkEvaluate",
          function(object, complexity=0, alpha=0, factor=1) {
            if(!is(object, "catNetworkEvaluate"))
               stop("catNetworkEvaluate object should be specified.")
            if(length(object@nets) < 1)
              return(NULL)
            if(complexity>0) {
              complx <- sapply(object@nets, function(pnet) pnet@complexity)
              idx <- which(complx == complexity)
              if(length(idx) == 0) {
                for(i in 1:length(complx))
                  if(complx[i] > complexity)
                    break
                idx <- i
              }
              return(object@nets[[idx]])
            }
            if(alpha == "BIC") 
              alpha <- -1
            if(alpha == "AIC") 
              alpha <- -2
            if(factor <= 0) {
              factor <- 1
              warning("factor set to 1")
            }
            if(alpha==-2) {##AIC
              score <- object@loglik - object@complexity
            }
            if(alpha==-1) {##BIC
              score <- object@loglik - 0.5*log(object@numsamples)*object@complexity
            }
            idx <- which(score == max(score))
            if(length(idx)!=1)
              return(NULL)
            return(object@nets[[idx]])
            })

setMethod("mgFindAIC", "catNetworkEvaluate", function(object) {
  objectlist <- object@nets
  if(length(objectlist) < 1)
    return(NULL)
  liststr <- ""
  maxobj <- objectlist[[1]]
  numsamples <- object@numsamples
  maxaic <- -Inf
  for(object in objectlist) {
    if(!is(object, "mgNetwork"))
      next
    curaic <- object@likelihood - object@complexity
    if(maxaic < curaic) {
      maxaic <- curaic
      maxobj <- object
    }
  }
  return(maxobj)
})

setMethod("mgFindBIC", "catNetworkEvaluate", function(object) {
  objectlist <- object@nets
  if(length(objectlist) < 1)
    return(NULL)
  liststr <- ""
  maxobj <- objectlist[[1]]
  maxbic <- -Inf
  numsamples <- object@numsamples
  for(object in objectlist) {
    if(!is(object, "mgNetwork"))
      next
    curbic <- object@likelihood - 0.5*object@complexity*log(numsamples)
    if(maxbic < curbic) {
      maxbic <- curbic
      maxobj <- object
    }
  }
  return(maxobj)
})

