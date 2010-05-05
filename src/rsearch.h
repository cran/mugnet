/*
 * rsearch.h
 *
 *  Created on:Nov 16, 2009
 *      Author: Nikolay Balov
 */

#ifndef RMIX_SEARCH_H
#define RMIX_SEARCH_H

#include "emsearch.h"

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

class RMixSearch : public EMSEARCH<double> {
public:
	RMixSearch();

	SEXP estimateNetworks(SEXP rSamples, SEXP rPerturbations,
                       SEXP rNodeCategories, SEXP rMaxParents, SEXP rMaxComplexity, 
		       SEXP rOrder,
                       SEXP rParentsPool, SEXP rFixedParentsPool, 
		       SEXP rEmIterations, SEXP rStopDelta, 
		       SEXP rNetSelection, SEXP rEcho);
};

#endif /* RMIX_SEARCH_H */
