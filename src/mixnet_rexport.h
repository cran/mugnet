
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/RConverters.h>
#include <R_ext/Rdynload.h>

/* these are called by R-functions directly */

SEXP marginalProb(SEXP cnet, SEXP rnode);
SEXP getProbSlot(SEXP cnet, SEXP rnode, SEXP rcats);
SEXP jointProb(SEXP cnet, SEXP rnode);
SEXP findParentPool(SEXP cnet, SEXP rnode);
SEXP parentSetProb(SEXP cnet, SEXP rnode, SEXP rparents, SEXP rSamples, SEXP rPerturbations);
SEXP searchOrder(
			SEXP rSamples, SEXP rPerturbations,
			SEXP rNodeCategories, SEXP rMaxParents, SEXP rMaxComplexity,
			SEXP rOrder, SEXP rParentsPool, SEXP rFixedParentsPool,
			SEXP rEmIterations, SEXP rStopDelta, 
			SEXP rNetSelection, SEXP rEcho);

SEXP quickSort(SEXP rlist);

SEXP sample(SEXP cnet, SEXP rNumSamples);
SEXP predict(SEXP cnet, SEXP rSamples);
SEXP nodeLoglik(SEXP cnet, SEXP rNodes, SEXP rSamples, SEXP rPerturbations);
