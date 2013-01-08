/*
 *  mugnet : Mixed Categorical Bayesian networks
 *  Copyright (C) 2009--2010  Nikolay Balov
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.gnu.org/licenses/gpl-2.0.html
 */

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
		       SEXP rEmIterations, SEXP rStopDelta, SEXP rEmStartIterations, 
		       SEXP rNetSelection, SEXP rModel, SEXP rEcho);
};

#endif /* RMIX_SEARCH_H */
