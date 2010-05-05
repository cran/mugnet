#########################################################################
# mgNetwork Class Methods 

setMethod("mgPredict", "mgNetwork", function(object, data) {

  if(!is.matrix(data) && !is.data.frame(data))
    stop("data should be a matrix or data frame of node categories")
  asframe <- FALSE
  if(is.data.frame(data)) {
    data <- as.matrix(t(data))
    asframe <- TRUE
  }
  
  if(length(dim(data)) == 2 && dim(data)[1] != object@numnodes)
    stop("the number of nodes in 'object' and 'data' should be equal")
  
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
  
  if(dim(data)[1] != object@numnodes)
    stop("Incompatible sample dimensions.\n")
   numsamples <- dim(data)[2]
   if(numsamples < 1)
     stop("No samples\n")

   newdata <- .Call("mgPredictC",
                    object, data,
                    PACKAGE="mugnet")
   if(is.null(newdata))
      return(NULL)

  newdata <- matrix(newdata, nrow = object@numnodes)	
  rownames(newdata)<-object@nodes
  if(asframe)
    return(as.data.frame(t(newdata)))
  return(newdata)
})

