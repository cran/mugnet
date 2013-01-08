#########################################################################
# Mixture of Gussians Network Class

setClass("mgNetwork",  
         representation(
                        model="character",
                        betas="list",
                        sigmas="vector"
	 ), 
	 contains="catNetwork", 
         validity = function(object) valid.mgNetwork(object),
         package = "mugnet"
	)

setGeneric("mgModel", 
          function(object)
           standardGeneric("mgModel")
	)

setGeneric("mgBeta", 
          function(object)
           standardGeneric("mgBeta")
	)

setGeneric("mgSigma", 
          function(object)
           standardGeneric("mgSigma")
	)

setGeneric("mgGaus", 
          function(object, cnet, betas, sigmas)
           standardGeneric("mgGaus")
	)

setGeneric("mgPois", 
          function(object, cnet, lambdas)
           standardGeneric("mgPois")
	)

setGeneric("mgExp", 
          function(object, cnet, lambdas)
           standardGeneric("mgExp")
	)

setGeneric("mgSamples", 
          function(object, numsamples=1, output="frame")
           standardGeneric("mgSamples")
	)

setGeneric("mgPredict", 
          function(object, data)
           standardGeneric("mgPredict")
	)

setGeneric("mgSetProb", 
          function(object, data)
           standardGeneric("mgSetProb")
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

setGeneric("mgFind", 
          function(object, complexity=0, alpha=0, factor=1)
           standardGeneric("mgFind")
	)
setGeneric("mgFindAIC", 
          function(object)
           standardGeneric("mgFindAIC")
	)
setGeneric("mgFindBIC", 
          function(object)
           standardGeneric("mgFindBIC")
	)
