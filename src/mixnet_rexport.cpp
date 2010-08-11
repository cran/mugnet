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
 * mixnet_rexport.cpp
 *
 *  Created on: Nov 16, 2009
 *      Author: Nikolay Balov
 */

#include "utils.h"
#include "rmixnet.h"
#include "rpoissonmix.h"
#include "rexpmix.h"
#include "rsearch.h"

extern "C" {

extern size_t g_memcounter;

SEXP searchOrder(
			SEXP rSamples, SEXP rPerturbations,
			SEXP rNodeCategories, SEXP rMaxParents, SEXP rMaxComplexity,
			SEXP rOrder, SEXP rParentsPool, SEXP rFixedParentsPool,
			SEXP rEmIterations, SEXP rStopDelta, SEXP rEmStartIterations, 
			SEXP rNetSelection, SEXP rModel, SEXP rEcho) {

	if(!isMatrix(rSamples))
		error("Samples should be a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations should be a matrix");
	if(!isInteger(AS_INTEGER(rNodeCategories)))
		error("NodeCategories should be integer");
	if(!isVector(AS_INTEGER(rNodeCategories)))
		error("NodeCategories should be a vector");
	if(!isInteger(AS_INTEGER(rMaxParents)))
		error("MaxParents should be an integer");
	if(!isInteger(AS_INTEGER(rMaxComplexity)))
		error("MaxComplexity should be an integer");
	if(!isVector(rOrder))
		error("rOrder should be a vector");
	if(!isNull(rParentsPool) && !isVector(rParentsPool))
		error("ParentsPool should be a list");
	if(!isNull(rFixedParentsPool) && !isVector(rFixedParentsPool))
		error("FixedParentsPool should be a list");
	if(!isInteger(AS_INTEGER(rEmIterations)))
		error("EmIterations should be an integer");
	if(!isNumeric(AS_NUMERIC(rStopDelta)))
		error("StopDelta should be numeric");

	RMixSearch * pengine = new RMixSearch;
	if(!pengine)
		CATNET_MEM_ERR();
	
	SEXP res = pengine -> estimateNetworks(
		rSamples, rPerturbations, 
		rNodeCategories, rMaxParents, rMaxComplexity, rOrder,
		rParentsPool, rFixedParentsPool, 
		rEmIterations, rStopDelta, rEmStartIterations, 
		rNetSelection, rModel, rEcho);
	
	delete pengine;

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return res;
}

SEXP sample(SEXP cnet, SEXP rModel, SEXP rNumSamples) {

	int res, numnodes, numsamples;
	double *pSamples, *pRsamples;
	SEXP rsamples;

	if(!isInteger(rNumSamples))
		error("rNumSamples should be integer");
	PROTECT(rNumSamples = AS_INTEGER(rNumSamples));
	numsamples = INTEGER_POINTER(rNumSamples)[0];
	UNPROTECT(1);

	rsamples = R_NilValue;
	
	I_NETPARAMS<double>* rnet = NULL;

	PROTECT(cnet);
	PROTECT(rModel = AS_CHARACTER(rModel));
	if(!strncmp(CHARACTER_VALUE(rModel), "Gaus", 4)) {
		rnet = new RMixNet(cnet);
	}
	else if(!strncmp(CHARACTER_VALUE(rModel), "Pois", 4)) {
		rnet = new RPoissonMix(cnet);
	}
	else if(!strncmp(CHARACTER_VALUE(rModel), "Exp", 3)) {
		rnet = new RExpMix(cnet);
	}
	else {
		UNPROTECT(2);
		error("Wrong model");
	}
	UNPROTECT(2);

	if(!rnet) {
		return R_NilValue;
	}

	numnodes = rnet->numNodes();
	
	pSamples = (double*)CATNET_MALLOC(numnodes*numsamples*sizeof(double));

	res = rnet->sample(pSamples, numsamples);

	rnet->releaseRef();

	if(res != ERR_CATNET_OK) {
		if(pSamples)
			CATNET_FREE(pSamples);
		return R_NilValue;
	}

	// output the new matrix
	PROTECT(rsamples = NEW_NUMERIC(numnodes * numsamples));
	pRsamples = NUMERIC_POINTER(rsamples);
	memcpy(pRsamples, pSamples, numnodes*numsamples*sizeof(double));
	UNPROTECT(1); // rsamples

	if(pSamples)
		CATNET_FREE(pSamples);
	return rsamples;
}

SEXP predict(SEXP cnet, SEXP rModel, SEXP rSamples) {

	int numnodes, numsamples;
	double *pRsamples, *pSamples;
	SEXP dim, rsamples;

	rsamples = R_NilValue;
	
	if(!isMatrix(rSamples))
		error("rSamples is not a matrix");
	
	I_NETPARAMS<double>* rnet = NULL;

	PROTECT(cnet);
	PROTECT(rModel = AS_CHARACTER(rModel));
	if(!strncmp(CHARACTER_VALUE(rModel), "Gaus", 4)) {
		rnet = new RMixNet(cnet);
	}
	else if(!strncmp(CHARACTER_VALUE(rModel), "Pois", 4)) {
		rnet = new RPoissonMix(cnet);
	}
	else if(!strncmp(CHARACTER_VALUE(rModel), "Exp", 3)) {
		rnet = new RExpMix(cnet);
	}
	else {
		UNPROTECT(2);
		error("Wrong model");
	}
	UNPROTECT(2);

	if(!rnet) {
		return R_NilValue;
	}

	PROTECT(rSamples = AS_NUMERIC(rSamples));
	pRsamples = NUMERIC_POINTER(rSamples);
	dim = GET_DIM(rSamples);
	numnodes = INTEGER(dim)[0];
	numsamples = INTEGER(dim)[1];
	pSamples = (double*)CATNET_MALLOC(numnodes*numsamples*sizeof(double));
	memcpy(pSamples, pRsamples, numnodes*numsamples*sizeof(double));
	UNPROTECT(1);

	rnet->predict(pSamples, numsamples);

	rnet->releaseRef();

	// output the new matrix
	PROTECT(rsamples = NEW_NUMERIC(numnodes * numsamples));
	pRsamples = NUMERIC_POINTER(rsamples);
	memcpy(pRsamples, pSamples, numnodes*numsamples*sizeof(double));
	UNPROTECT(1); // rsamples

	if(pSamples)
		CATNET_FREE(pSamples);

	return rsamples;
}

SEXP nodeLoglik(SEXP cnet, SEXP rModel, SEXP rNodes, SEXP rSamples, SEXP rPerturbations) {

	int *pPerturbations, *pNodes;
	double *pSamples;
	int numsamples, numnodes, nNodes, nnode, j;
	double *pvec;
	SEXP dim, rvec;

	if(!isMatrix(rSamples))
		error("rSamples should be a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("rPerturbations should be a vector");
	if(!isVector(rNodes))
		error("rNodes should be a vector");

	////////////////////////////////////////
	// Danger Ahead
	// We don's check that sample nodes actually correspond to the cnet's nodes
	// Missmatch of categories possible

	I_NETPARAMS<double>* rnet = NULL;

	PROTECT(cnet);
	PROTECT(rModel = AS_CHARACTER(rModel));
	if(!strncmp(CHARACTER_VALUE(rModel), "Gaus", 4)) {
		rnet = new RMixNet(cnet);
	}
	else if(!strncmp(CHARACTER_VALUE(rModel), "Pois", 4)) {
		rnet = new RPoissonMix(cnet);
	}
	else if(!strncmp(CHARACTER_VALUE(rModel), "Exp", 3)) {
		rnet = new RExpMix(cnet);
	}
	UNPROTECT(2);

	if(!rnet) {
		return R_NilValue;
	}

	PROTECT(rSamples = AS_NUMERIC(rSamples));
	pSamples = NUMERIC_POINTER(rSamples);
	dim = GET_DIM(rSamples);
	numnodes = INTEGER(dim)[0];
	numsamples = INTEGER(dim)[1];

	PROTECT(rNodes = AS_INTEGER(rNodes));
	pNodes = INTEGER(rNodes);
	nNodes = LENGTH(rNodes);

	PROTECT(rvec = NEW_NUMERIC(nNodes));
	pvec = NUMERIC_POINTER(rvec);

	dim = GET_DIM(rSamples);
	numnodes = INTEGER(dim)[0];
	numsamples = INTEGER(dim)[1];

	// nNodes are assumed positive indices
	for(j = 0; j < nNodes; j++) {
		if(pNodes[j] > 0)
			pNodes[j]--;
	}

	pPerturbations = 0;
	if(!isNull(rPerturbations)) {
		PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
		pPerturbations = INTEGER(rPerturbations);
		UNPROTECT(1);
	}

	// work with the sample distribution of pCurNet if necessary
	if(rnet->maxParentSet() > 4 || rnet->maxCategories()*rnet->maxParentSet() > 9 || 
		rnet->numNodes()*rnet->maxCategories()*rnet->maxParentSet() > 120) {
		j = (int)exp(log((double)(rnet->maxCategories()))*(1+rnet->maxParentSet())) * 100;
		//if(becho)
		//	printf("simulate probability with sample of size %d\n", j);
		rnet->catnetSample(j);
	}

	for(nnode = 0; nnode < nNodes; nnode++) {
		pvec[nnode] = rnet->findNodeLoglik(pNodes[nnode], pSamples, numsamples);
	}

	UNPROTECT(3);

	rnet->releaseRef();

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return rvec;
}

SEXP parentSetProb(SEXP cnet, SEXP rnode, SEXP rparents, SEXP rSamples, SEXP rPerturbations) {

	int i, j, node, ncats;
	int *pparents, nparents, numnodes, numsamples;
	double *pRsamples, *psamples;
	double *pvec, **pbetas, *psigmas;
	const double *pmarg;
	SEXP dim, rvec, plist;

	if(!isMatrix(rSamples))
		error("rSamples is not a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("rPerturbations is not a vector");
	if(!isInteger(AS_INTEGER(rnode)))
		error("rnode is not an integer");
	if(!isNull(rparents) && !isVector(rparents))
		error("rparents is not a vector");

	PROTECT(cnet);
	RMixNet *rnet = new RMixNet(cnet);
	UNPROTECT(1);

	PROTECT(rnode = AS_INTEGER(rnode));
	PROTECT(rparents = AS_INTEGER(rparents));
	PROTECT(rSamples = AS_NUMERIC(rSamples));

	node = INTEGER_VALUE(rnode) - 1;

	pparents = 0;
	nparents = length(rparents);
	if(nparents > 0) {
		pparents = (int*)CATNET_MALLOC(nparents*sizeof(int));
		for(i = 0; i < nparents; i++)
			pparents[i] = INTEGER_POINTER(rparents)[i] - 1;
	}

	dim = GET_DIM(rSamples);
	numnodes = INTEGER(dim)[0];
	numsamples = INTEGER(dim)[1];

	pbetas = (double**)CATNET_MALLOC(numnodes*sizeof(double*));
	psigmas = (double*)CATNET_MALLOC(numnodes*sizeof(double));
	for(i = 0; i < numnodes; i++) {
		pbetas[i] = (double*)CATNET_MALLOC(rnet->maxCategories()*sizeof(double));
		for(j = 0; j < rnet->maxCategories(); j++)
			pbetas[i][j] = j + 1;
		psigmas[i] = 0.1;
	}

	psamples = (double*)CATNET_MALLOC(numnodes*numsamples*sizeof(double));
	pRsamples = NUMERIC_POINTER(rSamples);

	for(j = 0; j < numsamples; j++) {
		for(i = 0; i < numnodes; i++) {
			psamples[j*numnodes + i] = pRsamples[j*numnodes + i];
		}
	}

	rnet->setBetas(pbetas, numnodes);
	rnet->setSigmas(psigmas, numnodes);

	for(i = 0; i < numnodes; i++) 
		CATNET_FREE(pbetas[i]);
	CATNET_FREE(pbetas);
	CATNET_FREE(psigmas);


	rnet->findNodeMarginalProb(node, psamples, numsamples);
	rnet->findNodeJointProb(node, pparents, nparents, psamples, numsamples);
	
	UNPROTECT(3);

	CATNET_FREE(pparents);
	CATNET_FREE(psamples);

	PROTECT(plist = allocVector(VECSXP, 4));

	// output pc's
	rvec = R_NilValue;
	ncats = rnet->pc_size();
	pmarg = rnet->pc();
	if(ncats > 0 && pmarg) {
		PROTECT(rvec = NEW_NUMERIC(ncats));
		pvec = NUMERIC_POINTER(rvec);
		for(i = 0; i < ncats; i++) 
			pvec[i] = pmarg[i];
		UNPROTECT(1);
	}
	SET_VECTOR_ELT(plist, 0, rvec);

	// output qc's
	rvec = R_NilValue;
	ncats = rnet->pc_size();
	pmarg = rnet->qc();
	if(ncats > 0 && pmarg) {
		PROTECT(rvec = NEW_NUMERIC(ncats));
		pvec = NUMERIC_POINTER(rvec);
		for(i = 0; i < ncats; i++) 
			pvec[i] = pmarg[i];
		UNPROTECT(1);
	}
	SET_VECTOR_ELT(plist, 1, rvec);

	// output pcC's
	rvec = R_NilValue;
	ncats = rnet->pcC_size();
	pmarg = rnet->pcC();
	if(ncats > 0 && pmarg) {
		PROTECT(rvec = NEW_NUMERIC(ncats));
		pvec = NUMERIC_POINTER(rvec);
		for(i = 0; i < ncats; i++) 
			pvec[i] = pmarg[i];
		UNPROTECT(1);
	}
	SET_VECTOR_ELT(plist, 2, rvec);

	// output qcC's
	rvec = R_NilValue;
	ncats = rnet->pcC_size();
	pmarg = rnet->qcC();
	if(ncats > 0 && pmarg) {
		PROTECT(rvec = NEW_NUMERIC(ncats));
		pvec = NUMERIC_POINTER(rvec);
		for(i = 0; i < ncats; i++) 
			pvec[i] = pmarg[i];
		UNPROTECT(1);
	}
	SET_VECTOR_ELT(plist, 3, rvec);

	UNPROTECT(1); // plist

	delete rnet;

	return plist;
}

SEXP marginalProb(SEXP cnet, SEXP rnode)
{
	int node, i, ncats;
	SEXP rvec = R_NilValue;
	double *pvec;
	int *prnodes;

	PROTECT(cnet);
	RMixNet *rnet = new RMixNet(cnet);
	UNPROTECT(1);

	if(!rnet)
		return rvec;

	PROTECT(rnode = AS_INTEGER(rnode));

	if(length(rnode) == 1) {
		node = INTEGER_VALUE(rnode);
		UNPROTECT(1);
		if(node < 1 || node > rnet->numNodes()) {
			delete rnet;
			return rvec;
		}
		node--;
		if(rnet->marginalProb(node) != ERR_CATNET_OK) {
			delete rnet;
			return rvec;
		}
	}
	else {
		prnodes = INTEGER_POINTER(rnode);
		int *pnodes = (int*)CATNET_MALLOC(length(rnode)*sizeof(int));
		for(i = 0; i < length(rnode); i++) {
			pnodes[i] = prnodes[i] - 1;
		}
		UNPROTECT(1);
		if(rnet->marginalProb(pnodes, length(rnode)) != ERR_CATNET_OK) {
			CATNET_FREE(pnodes);
			delete rnet;
			return rvec;
		}
		CATNET_FREE(pnodes);
	}

	double *pmarg = rnet->margProb();
	ncats = rnet->margProbSize(); //numCategories(node);

	PROTECT(rvec = NEW_NUMERIC(ncats));
	pvec = NUMERIC_POINTER(rvec);
	for(i = 0; i < ncats; i++) {
		pvec[i] = pmarg[i];
	}
	UNPROTECT(1);

	delete rnet;

	return rvec;
}

SEXP getProbSlot(SEXP cnet, SEXP rnode, SEXP rcats)
{
	int node, i;
	SEXP rvec = R_NilValue;
	int *prcats, *pcats;
	double *pvec;
	int *prnodes;

	PROTECT(cnet);
	RMixNet *rnet = new RMixNet(cnet);
	UNPROTECT(1);

	if(!rnet)
		return rvec;

	PROTECT(rnode = AS_INTEGER(rnode));
	if(length(rnode) == 1) {
		node = INTEGER_VALUE(rnode);
		UNPROTECT(1);
		if(node < 1 || node > rnet->numNodes()) {
			delete rnet;
			return rvec;
		}
		node--;
		if(rnet->marginalProb(node) != ERR_CATNET_OK) {
			delete rnet;
			return rvec;
		}
	}
	else {
		prnodes = INTEGER_POINTER(rnode);
		int *pnodes = (int*)CATNET_MALLOC(length(rnode)*sizeof(int));
		for(i = 0; i < length(rnode); i++) {
			pnodes[i] = prnodes[i] - 1;
		}
		UNPROTECT(1);
		if(rnet->marginalProb(pnodes, length(rnode)) != ERR_CATNET_OK) {
			delete rnet;
			return rvec;
		}
		CATNET_FREE(pnodes);
	}

	PROTECT(rcats = AS_INTEGER(rcats));

	prcats = INTEGER_POINTER(rcats);
	pcats = (int*)CATNET_MALLOC(length(rcats)*sizeof(int));
	for(i = 0; i < length(rcats); i++) {
		pcats[i] = prcats[i] - 1;
	}
	UNPROTECT(1);
	double pval = rnet->getCatProb(pcats, length(rcats));
	CATNET_FREE(pcats);

	PROTECT(rvec = NEW_NUMERIC(1));
	pvec = NUMERIC_POINTER(rvec);
	pvec[0] = pval;
	UNPROTECT(1);

	delete rnet;
	return rvec;
}


SEXP jointProb(SEXP cnet, SEXP rnode)
{
	int i, node, jointprobsize;
	SEXP rvec = R_NilValue;
	double *pvec;
	int *prnodes;

	PROTECT(cnet);
	RMixNet *rnet = new RMixNet(cnet);
	UNPROTECT(1);
	if(!rnet)
		return rvec;

	PROTECT(rnode = AS_INTEGER(rnode));

	if(length(rnode) == 1) {
		node = INTEGER_VALUE(rnode);
		UNPROTECT(1);
		if(node < 1 || node > rnet->numNodes()) {
			delete rnet;
			return rvec;
		}
		node--;
		if(rnet->findJointProb(node) != ERR_CATNET_OK) {
			delete rnet;
			return rvec;
		}
	}
	else {
		prnodes = INTEGER_POINTER(rnode);
		int *pnodes = (int*)CATNET_MALLOC(length(rnode)*sizeof(int));
		for(i = 0; i < length(rnode); i++) {
			pnodes[i] = prnodes[i] - 1;
		}
		UNPROTECT(1);
		if(rnet->findJointProb(pnodes, length(rnode)) != ERR_CATNET_OK) {
			delete rnet;
			return rvec;
		}
	}

	jointprobsize = rnet->jointProbSize();
	double *pjoint = rnet->jointProb();

	if(!pjoint) {
		delete rnet;
		return rvec;
	}

	PROTECT(rvec = NEW_NUMERIC(jointprobsize));
	pvec = NUMERIC_POINTER(rvec);
	memcpy(pvec, pjoint, jointprobsize*sizeof(double));
	UNPROTECT(1);

	delete rnet;
	return rvec;
}

SEXP findParentPool(SEXP cnet, SEXP rnode)
{
	int node, i, poolsize;
	SEXP rvec = R_NilValue;
	int *ppool, *pvec;


	PROTECT(cnet);
	RMixNet *rnet = new RMixNet(cnet);
	UNPROTECT(1);
	if(!rnet)
		return rvec;

	PROTECT(rnode = AS_INTEGER(rnode));

	//printf("length(rnode) = %d\n", length(rnode));

	if(length(rnode) == 1) {
		node = INTEGER_VALUE(rnode);
		UNPROTECT(1);
		if(node < 1 || node > rnet->numNodes())
			return rvec;
		node--;
		if(rnet->findParentPool(node) != ERR_CATNET_OK) {
			delete rnet;
			return rvec;
		}
	}
	else {
		pvec = INTEGER_POINTER(rnode);
		int *pnodes = (int*)CATNET_MALLOC(length(rnode)*sizeof(int));
		for(i = 0; i < length(rnode); i++) {
			pnodes[i] = pvec[i] - 1;
			//printf("%d: %d\n", i, pnodes[i]);
		}
		UNPROTECT(1);
		if(rnet->findParentPool(pnodes, length(rnode)) != ERR_CATNET_OK) {
			delete rnet;
			return rvec;
		}
	}

	poolsize = rnet->parPoolSize();
	ppool = rnet->parPool();

	PROTECT(rvec = NEW_INTEGER(poolsize));
	pvec = INTEGER_POINTER(rvec);
	for(i = 0; i < poolsize; i++) {
		pvec[i] = ppool[i] + 1;
	}
	UNPROTECT(1);

	delete rnet;
	return rvec;
}

char *gen_prob_string(int node, SEXP parlist, int paridx, SEXP catlist, SEXP problist, char *str) {
	int j, npar;
	SEXP parprobs, pcats;
	char *newstr, *aux, *aux2, *aux3;

	//char ss[128];

	if(!str) {
		str = (char*)CATNET_MALLOC(1);
		str[0] = 0;
	}

	if(paridx >= length(parlist)) {
		pcats = VECTOR_ELT(catlist, node);
		newstr = (char*)CATNET_MALLOC(((strlen(str)+1+32)*length(pcats))*sizeof(char));

		newstr[0] = 0;
		for(j = 0; j < length(pcats); j++) {
			sprintf(newstr, "%s%s%s %f\n", newstr, str, CHAR(STRING_ELT(pcats, j)), NUMERIC_POINTER(problist)[j]);
		}

		CATNET_FREE(str);
		str = newstr;
		return str;
	}

	npar = INTEGER_POINTER(parlist)[paridx] - 1;
	pcats = VECTOR_ELT(catlist, npar);

	newstr = (char*)CATNET_MALLOC(sizeof(char));
	newstr[0] = 0;
	for(j = 0; j < length(pcats); j++) {
		parprobs = VECTOR_ELT(problist, j);

		aux = (char*)CATNET_MALLOC((strlen(str)+1+8)*sizeof(char));

		//sprintf(aux, "%d, %d, %d,  %s\n", node, npar, length(parprobs), CHAR(STRING_ELT(pcats, j)));
		//printf(aux);

		sprintf(aux, "%s%s", str, CHAR(STRING_ELT(pcats, j)));
		aux2 = gen_prob_string(node, parlist, paridx + 1, catlist, parprobs, aux);

		aux3 = (char*)CATNET_MALLOC((strlen(newstr)+strlen(aux2)+2)*sizeof(char));
		sprintf(aux3, "%s%s", newstr, aux2);
		//sprintf(ss, "5 CATNET_FREE: %p, %d\n", newstr, strlen(newstr)); printf(ss);
		CATNET_FREE(newstr);
		newstr = aux3;

		CATNET_FREE(aux2);
	}
	CATNET_FREE(str);
	str = newstr;

	return str;
}

SEXP prob_string(SEXP rnodes, SEXP rparents, SEXP rcatlist, SEXP rproblist) {

	int node;
	SEXP nodepars, nodeproblist, pstr;
	char *str = NULL, *newstr, *aux;

	PROTECT(rnodes = AS_LIST(rnodes));
	PROTECT(rparents = AS_LIST(rparents));
	PROTECT(rcatlist = AS_LIST(rcatlist));
	PROTECT(rproblist = AS_LIST(rproblist));

	for(node = 0; node < length(rnodes); node++) {
		nodepars = VECTOR_ELT(rparents, node);
		nodeproblist = VECTOR_ELT(rproblist, node);
		newstr = gen_prob_string(node, nodepars, 0, rcatlist, nodeproblist, NULL);
		if(str) {
			aux = (char*)CATNET_MALLOC((strlen(str)+strlen(newstr)+1+16)*sizeof(char));
			sprintf(aux, "%sNode [%d]:\n%s", str, node, newstr);
			CATNET_FREE(str);
			CATNET_FREE(newstr);
			str = aux;
		}
		else
			str = newstr;
	}

	PROTECT(pstr = allocVector(STRSXP, 1));
	SET_STRING_ELT(pstr, 0, mkChar(str));

	UNPROTECT(5);
	return pstr;
}


void gen_prob_vector(int node, SEXP parlist, int paridx, SEXP catlist, SEXP problist, double *&pvec, int &nvec) {
	int j, npar;
	SEXP parprobs, pcats;
	double *newvec;

	if(!pvec) {
		pvec = (double*)CATNET_MALLOC(sizeof(double));
		nvec = 0;
	}

	if(paridx >= length(parlist)) {
		pcats = VECTOR_ELT(catlist, node);
		if (length(problist) != length(pcats)) {
			printf("gen_prob_vector: Format error in problist.");
			return;
		}
		newvec = (double*)CATNET_MALLOC((nvec + length(pcats))*sizeof(double));
		memcpy(newvec, pvec, nvec*sizeof(double));
		for(j = 0; j < length(pcats); j++) {
			newvec[nvec+j] = NUMERIC_POINTER(problist)[j];
		}
		CATNET_FREE(pvec);
		pvec = newvec;
		nvec += length(pcats);
		return;
	}

	npar = INTEGER_POINTER(parlist)[paridx] - 1;
	pcats = VECTOR_ELT(catlist, npar);
	if (length(problist) != length(pcats)) {
		printf("gen_prob_vector: Format error in problist.");
		return;
	}
	for(j = 0; j < length(pcats); j++) {
		parprobs = VECTOR_ELT(problist, j);
		gen_prob_vector(node, parlist, paridx + 1, catlist, parprobs, pvec, nvec);
	}
}

SEXP prob_vector(SEXP rnodes, SEXP rparents, SEXP rcatlist, SEXP rproblist) {

	int node, nvec;
	SEXP nodepars, nodeproblist, pstr;
	SEXP rvec;
	double *pvec, *prvec;

	PROTECT(rnodes = AS_LIST(rnodes));
	PROTECT(rparents = AS_LIST(rparents));
	PROTECT(rcatlist = AS_LIST(rcatlist));
	PROTECT(rproblist = AS_LIST(rproblist));

	PROTECT(pstr = allocVector(VECSXP, length(rnodes)));

	for(node = 0; node < length(rnodes); node++) {
		nodepars = VECTOR_ELT(rparents, node);
		nodeproblist = VECTOR_ELT(rproblist, node);
		pvec = 0;
		nvec = 0;
		gen_prob_vector(node, nodepars, 0, rcatlist, nodeproblist, pvec, nvec);
		PROTECT(rvec = NEW_NUMERIC(nvec));
		prvec = NUMERIC_POINTER(rvec);
		memcpy(prvec, pvec, nvec*sizeof(double));
		CATNET_FREE(pvec);
		SET_VECTOR_ELT(pstr, node, rvec);
		UNPROTECT(1);
	}

	UNPROTECT(5);
	return pstr;
}

SEXP quickSort(SEXP rlist) {
	int nlist;
	double *plist, *psortlist;
	SEXP rvec;

	PROTECT(rlist = AS_NUMERIC(rlist));
	nlist = length(rlist);
	if(nlist < 2)
		return R_NilValue;
	plist = NUMERIC_POINTER(rlist);
	_quick_sort<double>(plist, nlist);
	
	PROTECT(rvec = NEW_NUMERIC(nlist));
	psortlist = NUMERIC_POINTER(rvec);
	memcpy(psortlist, plist, nlist*sizeof(double));
	UNPROTECT(2);
	return rvec;
}

} // extern "C"
