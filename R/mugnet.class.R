#########################################################################
# Mixture of Gussians Network Class

setClass("mgNetwork",  
         representation(
                        betas="list",
                        sigmas="vector"
	 ), 
	 contains="catNetwork", 
         validity = function(object) valid.mgNetwork(object),
         package = "mugnet"
	)

setGeneric("mgBeta", 
          function(object)
           standardGeneric("mgBeta")
	)

setGeneric("mgSigma", 
          function(object)
           standardGeneric("mgSigma")
	)

setGeneric("mgFromCatnet", 
          function(object, cnet, betas, sigmas)
           standardGeneric("mgFromCatnet")
	)

setGeneric("mgSamples", 
          function(object, numsamples=1, output="frame")
           standardGeneric("mgSamples")
	)

setGeneric("mgPredict", 
          function(object, data)
           standardGeneric("mgPredict")
	)

setGeneric("mgLoglik", 
          function(object, data, bysample=FALSE)
           standardGeneric("mgLoglik")
	)

setGeneric("mgNodeLoglik", 
          function(object, nodes, data)
           standardGeneric("mgNodeLoglik")
	)

setGeneric("mgReorderNodes", 
          function(object, nodeIndices)
           standardGeneric("mgReorderNodes")
	)
