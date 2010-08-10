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
			SEXP rEmIterations, SEXP rStopDelta, SEXP rEmStartIterations, 
			SEXP rNetSelection, SEXP rEcho);

SEXP quickSort(SEXP rlist);

SEXP sample(SEXP cnet, SEXP rModel, SEXP rNumSamples);
SEXP predict(SEXP cnet, SEXP rModel, SEXP rSamples);
SEXP nodeLoglik(SEXP cnet, SEXP rModel, SEXP rNodes, SEXP rSamples, SEXP rPerturbations);
