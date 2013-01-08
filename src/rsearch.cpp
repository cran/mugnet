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
 * rsearch.cpp
 *
 *  Created on: Nov 16, 2009
 *      Author: Nikolay Balov
 */

#include "utils.h"
#include "rsearch.h"
#include "rmixnet.h"
#include "rpoissonmix.h"
#include "rexpmix.h"

RMixSearch::RMixSearch() {
}

SEXP RMixSearch::estimateNetworks(SEXP rSamples, SEXP rPerturbations, 
                       SEXP rNodeCategories, SEXP rMaxParents, SEXP rMaxComplexity, 
		       SEXP rOrder,
                       SEXP rParentsPool, SEXP rFixedParentsPool, 
		       SEXP rEmIterations, SEXP rStopDelta, SEXP rEmStartIterations, 
		       SEXP rNetSelection, SEXP rModel, SEXP rEcho) {

	int i, j, k, len, numnodes, numsamples, maxParentSet, 
		maxComplexity, emIterations, emStartIterations, nnode, 
		numnets, inet, netSelection, echo;
	int maxCategories, *pNodeCategories;
 	int *pRperturbations, *pPerturbations, **parentsPool, **fixedParentsPool, *pRorder, *pPool;
	double *pRsamples, *pSamples, **pBetas, *pSigmas, *pNodeSamples, range, stopDelta;
	char model[256];

	SEXP dim, rparpool, cnetlist, cnetnode;

	if(!isMatrix(rSamples))
		error("Samples is not a real matrix");

	PROTECT(rSamples = AS_NUMERIC(rSamples));

	dim = GET_DIM(rSamples);
	numnodes = INTEGER(dim)[0];
	numsamples = INTEGER(dim)[1];

	if(numnodes <= 1 || numsamples < 1) {
		UNPROTECT(1);
		error("Insufficient sample");
		return R_NilValue;
	}

	PROTECT(rMaxParents = AS_INTEGER(rMaxParents));
	PROTECT(rMaxComplexity = AS_INTEGER(rMaxComplexity));
	PROTECT(rEmIterations = AS_INTEGER(rEmIterations));
	PROTECT(rStopDelta = AS_NUMERIC(rStopDelta));
	PROTECT(rEmStartIterations = AS_INTEGER(rEmStartIterations));
	PROTECT(rNetSelection = AS_INTEGER(rNetSelection));
	PROTECT(rEcho = AS_LOGICAL(rEcho));

	maxParentSet = INTEGER_POINTER(rMaxParents)[0];
	maxComplexity = INTEGER_POINTER(rMaxComplexity)[0];

	emIterations = INTEGER_POINTER(rEmIterations)[0];
	if(emIterations < 2)
		emIterations = 2;
	stopDelta = NUMERIC_POINTER(rStopDelta)[0];
	if(stopDelta < 0)
		stopDelta = 0;
	emStartIterations = INTEGER_POINTER(rEmStartIterations)[0];
	if(emStartIterations < 2)
		emStartIterations = 2;

	netSelection = INTEGER_POINTER(rNetSelection)[0];

	echo = LOGICAL(rEcho)[0];

	UNPROTECT(7);

	PROTECT(rNodeCategories = AS_INTEGER(rNodeCategories));
	if(length(rNodeCategories) != numnodes) {
		UNPROTECT(1);
		error("length(NodeCategories) != numnodes");
	}

	maxCategories = 0;
	pNodeCategories = (int*)CATNET_MALLOC(numnodes*sizeof(int));
	// rNodeCategories is ordered according to the data
	memcpy(pNodeCategories, INTEGER_POINTER(rNodeCategories), numnodes*sizeof(int));
	for(i = 0; i < numnodes; i++) {
		if(pNodeCategories[i] < 2) {
			pNodeCategories[i] = 2;
			warning("Set pNodeCategories = 2\n");
		}
		if(maxCategories < pNodeCategories[i])
			maxCategories = pNodeCategories[i];
 	}

	pSamples = (double*)CATNET_MALLOC(numnodes*numsamples*sizeof(double));
	pRsamples = NUMERIC_POINTER(rSamples);

	PROTECT(rOrder = AS_INTEGER(rOrder));
	pRorder = (int*)CATNET_MALLOC(numnodes*sizeof(int));
	if(length(rOrder) != numnodes) {
		warning("Invalid rOrder parameter. Reset to default node order");
		for(i = 0; i < numnodes; i++)
			pRorder[i] = i + 1;
	}
	else
		memcpy(pRorder, INTEGER(rOrder), numnodes*sizeof(int));
	// rOrder
	UNPROTECT(1);

	pBetas = (double**)CATNET_MALLOC(numnodes*sizeof(double*));
	pSigmas = (double*)CATNET_MALLOC(numnodes*sizeof(double));

	pNodeSamples = (double*)CATNET_MALLOC(numsamples*sizeof(double));

	if(echo)
		Rprintf("Starting parameters:\n");

	for(i = 0; i < numnodes; i++) {
		nnode = pRorder[i] - 1;
		// reorder
		pNodeCategories[i] = INTEGER_POINTER(rNodeCategories)[nnode];
		for(j = 0; j < numsamples; j++) {
			pNodeSamples[j] = pRsamples[j*numnodes + nnode];
			pSamples[j*numnodes + i] = pNodeSamples[j];
		}
		_quick_sort<double>(pNodeSamples, numsamples);
		range = pNodeSamples[numsamples-1] - pNodeSamples[0];
		if(echo)
			Rprintf("  beta[%d] = ", i);
		
		pBetas[i] = (double*)CATNET_MALLOC(pNodeCategories[i]*sizeof(double));
		for(j = 0; j < pNodeCategories[i]; j++) {
			/* uniform */
			pBetas[i][j] = pNodeSamples[0] + range*(j+1)/(pNodeCategories[i]+1);
			/* quantiles */
			//k = (2*j+1)*numsamples / (2*pNodeCategories[i]);
			//pBetas[i][j] = pNodeSamples[(int)k];
			if(echo)
				Rprintf("%.4f  ", pBetas[i][j]);
		}
		if(echo)
			Rprintf("\n");
		pSigmas[i] = range / (2*pNodeCategories[i]);
		pSigmas[i] = pSigmas[i]*pSigmas[i];
		if(echo)
			Rprintf("    sigma[%d] = %.4f\n", i, pSigmas[i]);
	}

	// rSamples & rNodeCategories
	UNPROTECT(2);

	if(pNodeSamples)
		CATNET_FREE(pNodeSamples);

	pPerturbations = 0;
	if(!isNull(rPerturbations)) {
		PROTECT(rPerturbations = AS_INTEGER(rPerturbations));

		pPerturbations = (int*)CATNET_MALLOC(numnodes*numsamples*sizeof(int));
		pRperturbations = INTEGER(rPerturbations);
		for(j = 0; j < numsamples; j++) {
			for(i = 0; i < numnodes; i++) {
				pPerturbations[j*numnodes + i] = 
				pRperturbations[j*numnodes + pRorder[i] - 1];
			}
		}
		UNPROTECT(1);
	}

	parentsPool = 0;
	if(!isNull(rParentsPool)) {
		PROTECT(rParentsPool = AS_LIST(rParentsPool));

		parentsPool = (int**)CATNET_MALLOC(numnodes*sizeof(int*));
		memset(parentsPool, 0, numnodes*sizeof(int*));
		for(i = 0; i < numnodes; i++) {
			rparpool = AS_INTEGER(VECTOR_ELT(rParentsPool, (int)(pRorder[i] - 1)));
			len = length(rparpool);
			if(isVector(rparpool) && len > 0 && len <= numnodes) {
				pPool = INTEGER(rparpool);
				parentsPool[i] = (int*)CATNET_MALLOC(numnodes*sizeof(int));
				for(j = 0; j < len; j++) {
					if(pPool[j] > 0 && pPool[j] <= numnodes) {
						for(k = 0; k < numnodes; k++)
							if(pPool[j] == pRorder[k])
								break;
						if(k < numnodes)
							parentsPool[i][j] = k;
						else
							parentsPool[i][j] = -1;
					}
				}
				for(; j < numnodes; j++)
					parentsPool[i][j] = -1;
			}
		}
		UNPROTECT(1);
	}

	fixedParentsPool = 0;
	if(!isNull(rFixedParentsPool)) {
		PROTECT(rFixedParentsPool = AS_LIST(rFixedParentsPool));

		fixedParentsPool = (int**)CATNET_MALLOC(numnodes*sizeof(int*));
		memset(fixedParentsPool, 0, numnodes*sizeof(int*));
		for(i = 0; i < numnodes; i++) {
			rparpool = AS_INTEGER(VECTOR_ELT(rFixedParentsPool, (int)(pRorder[i] - 1)));
			len = length(rparpool);

			if(isVector(rparpool) && len > 0 && len <= numnodes) {
			 	if(maxParentSet < len)
			    		maxParentSet = len;
				pPool = INTEGER(rparpool);
				fixedParentsPool[i] = (int*)CATNET_MALLOC(numnodes*sizeof(int));
				for(j = 0; j < len; j++) {
					if(pPool[j] > 0 && pPool[j] <= numnodes) {
						for(k = 0; k < numnodes; k++)
							if(pPool[j] == pRorder[k])
								break;
						if(k < numnodes)
							fixedParentsPool[i][j] = k;
						else
							fixedParentsPool[i][j] = -1;
					}
				}
				for(; j < numnodes; j++)
					fixedParentsPool[i][j] = -1;
			}

		}
		UNPROTECT(1);
	}

	PROTECT(rModel = AS_CHARACTER(rModel));
	if(echo)
		Rprintf("model = %s\n", CHARACTER_VALUE(rModel));
	// make an egg
	I_NETPARAMS<double>* pRnormnet = NULL;
	if(!strncmp(CHARACTER_VALUE(rModel), "Gaus", 4)) {
		strcpy(model, "Gaus");
		// make sure maxComplexity is at least the minimum possible complexity
		if(numnodes*(2*maxCategories) > maxComplexity) {
			maxComplexity = numnodes*(2*maxCategories);
			warning("Set maxComplexity = %d\n", maxComplexity);
		}

		pRnormnet = (I_NETPARAMS<double>*)new RMixNet(numnodes, maxParentSet, maxCategories, 
				/*const t_node **nodes =*/ 0, /*const int * pnumpars =*/ 0,
				/*const int **ppars =*/ 0, /*const int *pcats =*/ pNodeCategories);
	}
	else if(!strncmp(CHARACTER_VALUE(rModel), "Pois", 4)) {
		strcpy(model, "Pois");
		// make sure maxComplexity is at least the minimum possible complexity
		if(numnodes*(2*maxCategories-1) > maxComplexity) {
			maxComplexity = numnodes*(2*maxCategories-1);
			warning("Set maxComplexity = %d\n", maxComplexity);
		}

		pRnormnet = (I_NETPARAMS<double>*)new RPoissonMix(numnodes, maxParentSet, maxCategories, 
				/*const t_node **nodes =*/ 0, /*const int * pnumpars =*/ 0,
				/*const int **ppars =*/ 0, /*const int *pcats =*/ pNodeCategories);
	}
	else if(!strncmp(CHARACTER_VALUE(rModel), "Exp", 3)) {
		strcpy(model, "Exp");
		// make sure maxComplexity is at least the minimum possible complexity
		if(numnodes*(2*maxCategories-1) > maxComplexity) {
			maxComplexity = numnodes*(2*maxCategories-1);
			warning("Set maxComplexity = %d\n", maxComplexity);
		}

		pRnormnet = (I_NETPARAMS<double>*)new RExpMix(numnodes, maxParentSet, maxCategories, 
				/*const t_node **nodes =*/ 0, /*const int * pnumpars =*/ 0,
				/*const int **ppars =*/ 0, /*const int *pcats =*/ pNodeCategories);
	}
	else {
		UNPROTECT(1);
		error("Wrong model");
	}

	UNPROTECT(1);
	if(!pRnormnet)
		CATNET_MEM_ERR();

	CATNET_FREE(pNodeCategories);
	pNodeCategories = NULL;

	pRnormnet->setBetas(pBetas, numnodes);
	pRnormnet->setSigmas(pSigmas, numnodes);

	for (i = 0; i < numnodes; i++) {
		if(pBetas[i])
			CATNET_FREE(pBetas[i]);
		pBetas[i] = 0;
	}
	CATNET_FREE(pBetas);
	pBetas = 0;
	CATNET_FREE(pSigmas);
	pSigmas = 0; 

	estimate((I_NETPARAMS<double>*) pRnormnet, numsamples, pSamples, pPerturbations, maxParentSet, 
			maxComplexity, parentsPool, fixedParentsPool, 
			emIterations, stopDelta, emStartIterations, netSelection, echo);

	pRnormnet->releaseRef();
	pRnormnet = 0;

	// create a R-list of catNetworks
	numnets = 0;
	for(i = 0; i < m_nNets; i++) {
		if(m_pCurNets[i]) {
			m_pCurNets[i]->setNodesOrder(pRorder);
			//loglik = 0;
			//for(k = 0; k < numnodes; k++)
			//	loglik += m_pCurNets[i]->findNodeLoglik(k, pSamples, numsamples);
			//m_pCurNets[i]->setLoglik(loglik);
			numnets++;
		}
	}
	
	if(pSamples)
		CATNET_FREE(pSamples);
	if(pPerturbations)
		CATNET_FREE(pPerturbations);
	if(parentsPool) {
		for(i = 0; i < numnodes; i++)
			if(parentsPool[i])
				CATNET_FREE(parentsPool[i]);
		CATNET_FREE(parentsPool);
	}
	if(fixedParentsPool) {
		for(i = 0; i < numnodes; i++)
			if(fixedParentsPool[i])
				CATNET_FREE(fixedParentsPool[i]);
		CATNET_FREE(fixedParentsPool);
	}

	if(pRorder)
		CATNET_FREE(pRorder);

	if(!numnets || !m_pCurNets) {
		return R_NilValue;
	}

	if(echo)
		Rprintf("Found %d networks\n", numnets);

	PROTECT(cnetlist = allocVector(VECSXP, numnets));

	inet = 0;
	cnetnode = 0;
	for(i = 0; i < m_nNets; i++) {
		if(!m_pCurNets[i])
			continue;
		if(!strncmp(model, "Gaus", 4)) 
			PROTECT(cnetnode = ((RMixNet*)m_pCurNets[i])->genRmixnet("mgNetwork"));
		if(!strncmp(model, "Pois", 4)) 
			PROTECT(cnetnode = ((RPoissonMix*)m_pCurNets[i])->genRmixnet("mgNetwork"));
		if(!strncmp(model, "Exp", 3)) 
			PROTECT(cnetnode = ((RExpMix*)m_pCurNets[i])->genRmixnet("mgNetwork"));
		if(cnetnode)
			SET_VECTOR_ELT(cnetlist, inet, cnetnode);
		UNPROTECT(1);
		inet++;
	}

	UNPROTECT(1);

	return cnetlist;
}

