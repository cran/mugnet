#########################################################################
# mgNetwork Class Methods 

setMethod("mgSamples", "mgNetwork", function(object, numsamples = 1, output="frame") {
  if(missing(numsamples))
    numsamples <- 1
  if(numsamples < 1)	
    numsamples <- 1	
  data <- .Call("mgSampleC",
                object, as.character(object@model), as.integer(numsamples),
                PACKAGE="mugnet")
  if(is.null(data))
    return(NULL)
  data <- matrix(data, nrow = object@numnodes)	
  rownames(data)<-object@nodes
  if(!missing(output) && output=="matrix")
    return(data)
  return(as.data.frame(t(data)))
})

setMethod("mgLoglik", c("mgNetwork"), 
          function(object, data, bysample=FALSE) {

            if(!is.matrix(data) && !is.data.frame(data))
              stop("'data' should be a matrix or data frame")

            if(bysample==TRUE)
              warning("bysample=TRUE is ignored")
            
	    if(is.data.frame(data))
              data <- as.matrix(t(data))

            if(!is.numeric(data))
              stop("'data' should be numeric")
            
            if(length(dim(data)) == 2 && dim(data)[1] != object@numnodes)
              stop("The number of nodes in  the object and data should be equal")
            
            rownames <- rownames(data)
            if(length(rownames) != object@numnodes)
              stop("The data rows should be named after the nodes of the object.")
            
            if(prod(tolower(rownames) == tolower(object@nodes)) == 0) {
              norder <- order(rownames)
              data <- data[norder,]
              rownames <- rownames(data)
              norder <- order(object@nodes)
              object <- cnReorderNodes(object, norder)
            }
            
            if(prod(tolower(rownames) == tolower(object@nodes)) == 0)
              stop("The row names should correspond to the object nodes.")
            
            numnodes <- dim(data)[1]
            numsamples <- dim(data)[2]
            
            vloglik <- .Call("mgNodeLoglikC", 
                             object, as.character(object@model),
                             1:numnodes, data, NULL, 
                             PACKAGE="mugnet")
            if(length(vloglik)<1)
              return(-Inf)
            return(sum(vloglik))
          })


setMethod("mgNodeLoglik", c("mgNetwork"), 
          function(object, nodes, data) {

            if(!is.matrix(data) && !is.data.frame(data))
              stop("'data' should be a matrix or data frame")
            
            if(is.data.frame(data))
              data <- as.matrix(t(data))

            if(!is.numeric(data))
              stop("'data' should be numeric")

            numnodes <- dim(data)[1]
            numsamples <- dim(data)[2]  
            if(numsamples < 1)
              stop("No samples\n")
            
            rownames <- rownames(data)
            if(length(rownames) != object@numnodes)
              stop("The data rows should be named after the nodes of the object.")
            
            if(prod(tolower(rownames) == tolower(object@nodes)) == 0) {
              norder <- order(rownames)
              data <- matrix(data[norder,], nrow=numnodes)
              rownames <- rownames(data)
              norder <- order(object@nodes)
              object <- cnReorderNodes(object, norder)
            }
            
            if(prod(tolower(rownames) == tolower(object@nodes)) == 0)
              stop("The row names should correspond to the object nodes.")

            if(is.character(nodes)) {
              nodes <- sapply(nodes, function(cnode) which(object@nodes == cnode))
            }
            nodes <- as.integer(nodes)
            if(!is.vector(nodes) || length(nodes) < 1)
              stop("'nodes' should be a vector of node indices")

            id <- which(nodes > numnodes)
            if(length(id) > 0)
              nodes <- nodes[-id]

            if(length(nodes) < 1)
              return(NULL)
            vloglik <- .Call("mgNodeLoglikC", 
                             object, as.character(object@model),
                             nodes, data, NULL, 
                             PACKAGE="mugnet")
            return(vloglik)
          }
          )
