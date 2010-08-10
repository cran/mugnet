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
 * rmixnet.cpp
 *
 *  Created on: Sep 21, 2009
 *      Author: Nikolay Balov
 */

#include "utils.h"
#include "rmixnet.h"
#include "rpoissonmix.h"

extern "C" {

char *gen_prob_string(int node, SEXP parlist, int paridx, SEXP catlist, SEXP problist, char *str);
SEXP prob_string(SEXP rnodes, SEXP rparents, SEXP rcatlist, SEXP rproblist);
void gen_prob_vector(int node, SEXP parlist, int paridx, SEXP catlist, SEXP problist, double *&pvec, int &nvec);
SEXP prob_vector(SEXP rnodes, SEXP rparents, SEXP rcatlist, SEXP rproblist); 

} // extern "C"


RPoissonMix::RPoissonMix(SEXP cnet) {

	double *pvec, *pbeta;
	int nnode, i, nvec, *pn;
	char const *pstr;
	//char str[256];

	SEXP rname, rnodes, rparents, rcatlist, rproblist, pf, 
		nodepars, nodeproblist, pint, rnodeprob, 
		rbetas, rsigmas;

	if(!isS4(cnet))
		return;

	rname = GET_SLOT(cnet, install("objectName"));
	rnodes = GET_SLOT(cnet, install("nodes"));
	rparents = GET_SLOT(cnet, install("parents"));
	rcatlist = GET_SLOT(cnet, install("categories"));
	rproblist = GET_SLOT(cnet, install("probabilities"));

	if(rnodes == R_NilValue || rparents == R_NilValue || rcatlist == R_NilValue || rproblist == R_NilValue) {
		return;
	}

	PROTECT(rname = AS_CHARACTER(rname));
	PROTECT(rnodes = AS_LIST(rnodes));
	PROTECT(rparents = AS_LIST(rparents));
	PROTECT(rcatlist = AS_LIST(rcatlist));
	PROTECT(rproblist = AS_LIST(rproblist));

	pint = GET_SLOT(cnet, install("numnodes"));
	m_numNodes = INTEGER_POINTER(pint)[0];
	pint = GET_SLOT(cnet, install("maxParents"));
	m_maxParents = INTEGER_POINTER(pint)[0];
	pint = GET_SLOT(cnet, install("maxCategories"));
	m_maxCategories = INTEGER_POINTER(pint)[0];

	pint = GET_SLOT(cnet, install("complexity"));
	m_complexity = INTEGER_POINTER(pint)[0];
	pint = GET_SLOT(cnet, install("likelihood"));
	m_loglik = NUMERIC_POINTER(pint)[0];

	if (length(rproblist) != m_numNodes) {
		UNPROTECT(5);
		warning("length(rproblist) != m_numNodes");
		return;
	}
	m_nodeNames = (char**) CATNET_MALLOC(m_numNodes * sizeof(char*));
	memset(m_nodeNames, 0, m_numNodes * sizeof(char*));
	m_numParents = (int*) CATNET_MALLOC(m_numNodes * sizeof(int));
	memset(m_numParents, 0, m_numNodes * sizeof(int));
	m_parents = (int**) CATNET_MALLOC(m_numNodes * sizeof(int*));
	memset(m_parents, 0, m_numNodes * sizeof(int*));
	m_numCategories = (int*) CATNET_MALLOC(m_numNodes * sizeof(int));
	memset(m_numCategories, 0, m_numNodes * sizeof(int));
	m_pProbLists = (PROB_LIST<double>**) CATNET_MALLOC(m_numNodes * sizeof(PROB_LIST<double>*));
	for (i = 0; i < m_numNodes; i++)
		m_pProbLists[i] = 0;
	
	for(nnode = 0; nnode < m_numNodes; nnode++) {
		pf = VECTOR_ELT(rnodes, nnode);
		m_nodeNames[nnode] = 0;
		if(IS_VECTOR(pf)) {
			pstr = CHAR(asChar(pf));
			if(strlen(pstr) >= MAX_NODE_NAME) {
				m_nodeNames[nnode] = (char*) CATNET_MALLOC((strlen(pstr)+1) * sizeof(char));
				strcpy(m_nodeNames[nnode], pstr);
			}
			else {
				m_nodeNames[nnode] = (char*) CATNET_MALLOC(MAX_NODE_NAME * sizeof(char));
				sprintf(m_nodeNames[nnode], "N%d", nnode);
			}
		}

		pf = VECTOR_ELT(rparents, nnode);
		m_numParents[nnode] = 0;
		m_parents[nnode] = 0;
		if (IS_VECTOR(pf)) {
			m_numParents[nnode] = length(pf);
			if(m_maxParents < m_numParents[nnode])
				m_maxParents = m_numParents[nnode];
			pn = INTEGER_POINTER(pf);
			m_parents[nnode] = (int*) CATNET_MALLOC(m_numParents[nnode] * sizeof(int));
			for(i = 0; i < m_numParents[nnode]; i++) {
				m_parents[nnode][i] = pn[i] - 1;
			}
		}

		pf = VECTOR_ELT(rcatlist, nnode);
		m_numCategories[nnode] = length(pf);
	}

	// get probabilities
	if(!strcmp(CHAR(asChar(rname)), "mgNetwork") || 
		!strcmp(CHAR(asChar(rname)), "catNetwork")) {
		for (nnode = 0; nnode < m_numNodes; nnode++) {
			nodepars = VECTOR_ELT(rparents, nnode);
			nodeproblist = VECTOR_ELT(rproblist, nnode);
			pvec = 0;
			nvec = 0;
			gen_prob_vector(nnode, nodepars, 0, rcatlist, nodeproblist, pvec, nvec);
			//sprintf(str, "%d, %d,  %p, %d\n", nnode, m_numParents[nnode], pvec, nvec);
			//cout << str;
			setCondProb(nnode, pvec, nvec);
			CATNET_FREE(pvec);
		}
	}
	else if(!strcmp(CHAR(asChar(rname)), "catNetworkC")) {
		for (nnode = 0; nnode < m_numNodes; nnode++) {
			rnodeprob = VECTOR_ELT(rproblist, nnode);
			setCondProb(nnode, NUMERIC_POINTER(rnodeprob), length(rnodeprob));   
		}
	}

	UNPROTECT(5);

	rbetas = GET_SLOT(cnet, install("betas"));
	rsigmas = GET_SLOT(cnet, install("sigmas"));
	PROTECT(rbetas = AS_LIST(rbetas));
	PROTECT(rsigmas = AS_NUMERIC(rsigmas));

	if(GET_LENGTH(rbetas) == m_numNodes && GET_LENGTH(rsigmas) == m_numNodes) {
		for(nnode = 0; nnode < m_numNodes; nnode++) {
			pf = VECTOR_ELT(rbetas, nnode);
			pbeta = NUMERIC_POINTER(pf);
			setNodeBetas(nnode, pbeta);
		}
		pbeta = NUMERIC_POINTER(rsigmas);
		setSigmas(pbeta, m_numNodes);
	}

	UNPROTECT(2);

	//printf("New RMixNet %d, %d, %d\n", m_numNodes, m_maxParents, m_maxCategories);
}

SEXP RPoissonMix::genRcatnet(const char * objectName = (const char*)"catNetwork") {

	char str[256];
	int node, i, *pslotcats, *pn;
	SEXP plist, ppars, pcats, pnodeprobs, strnames, pint;

	if(strcmp(objectName, "catNetwork")) 
		return R_NilValue;

	SEXP cnet = PROTECT(NEW_OBJECT(MAKE_CLASS(objectName)));

	PROTECT(plist = allocVector(STRSXP, 1));
	SET_STRING_ELT(plist, 0, mkChar(objectName));
	SET_SLOT(cnet, install("objectName"), plist);
	UNPROTECT(1);

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = m_numNodes;
	SET_SLOT(cnet, install("numnodes"), pint);
	UNPROTECT(1);

	PROTECT(strnames = allocVector(STRSXP, m_numNodes));
	for(node = 0; node < m_numNodes; node++) {
		if(m_nodeNames && m_nodeNames[node])
			SET_STRING_ELT(strnames, node, mkChar(m_nodeNames[node]));
		else {
			sprintf(str, "N%d", node+1);
			SET_STRING_ELT(strnames, node, mkChar(str));
		}
	}
	SET_SLOT(cnet, install("nodes"), strnames);
 	UNPROTECT(1);

	PROTECT(plist = allocVector(STRSXP, m_numNodes));
	for(node = 0; node < m_numNodes; node++)
		SET_STRING_ELT(plist, node, mkChar(""));
	SET_SLOT(cnet, install("color"), plist);
	UNPROTECT(1);

	SET_SLOT(cnet, install("meta"), mkString("catNetwork object"));

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = m_maxParents;
	SET_SLOT(cnet, install("maxParents"), pint);
	UNPROTECT(1);

	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	for(node = 0; node < m_numNodes; node++) {
		if(m_numParents[node] <= 0) {
			SET_VECTOR_ELT(plist, node, R_NilValue);
			continue;
		}
		PROTECT(ppars = NEW_INTEGER(m_numParents[node]));
		pn = INTEGER_POINTER(ppars);
		for(i = 0; i < m_numParents[node]; i++)
			// remember to increase the index by 1
			pn[i] = m_parents[node][i] + 1;
		SET_VECTOR_ELT(plist, node, ppars);
		UNPROTECT(1);
	}
	SET_SLOT(cnet, install("parents"), plist);
	UNPROTECT(1);

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = m_maxCategories;
	SET_SLOT(cnet, install("maxCategories"), pint);
	UNPROTECT(1);

	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	for(node = 0; node < m_numNodes; node++) {
		if(m_numCategories[node] > m_maxCategories)
			break;
		PROTECT(pcats = allocVector(STRSXP, m_numCategories[node]));
		for(i = 0; i < m_numCategories[node]; i++) {
			if(m_catIndices && m_catIndices[node])
				sprintf(str, "%d", m_catIndices[node][i]);
			else
				sprintf(str, "C%d", i+1);
			SET_STRING_ELT(pcats, i, mkChar(str));
		}
		SET_VECTOR_ELT(plist, node, pcats);
		UNPROTECT(1);
	}
	SET_SLOT(cnet, install("categories"), plist);
	UNPROTECT(1);

	pslotcats = (int*)CATNET_MALLOC(m_maxParents*sizeof(int));
	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	for(node = 0; node < m_numNodes; node++) {
		memset(pslotcats, 0, m_maxParents*sizeof(int));
		pnodeprobs = genProbList(node, 0, pslotcats);
		SET_VECTOR_ELT(plist, node, pnodeprobs);
		if(pnodeprobs != R_NilValue)
			UNPROTECT(1);
	}
	SET_SLOT(cnet, install("probabilities"), plist);
	UNPROTECT(1);
	CATNET_FREE(pslotcats);

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = complexity();
	SET_SLOT(cnet, install("complexity"), pint);
	UNPROTECT(1);
	
	PROTECT(pint = NEW_NUMERIC(1));
	NUMERIC_POINTER(pint)[0] = getLoglik();
	SET_SLOT(cnet, install("likelihood"), pint);
	UNPROTECT(1);

	UNPROTECT(1); // cnet

	return cnet;
}

SEXP RPoissonMix::genProbList(int node, int paridx, int *pcats) {
	int j, npar;
	SEXP problist;
	double *pslot, *pp;

	if(m_pProbLists == 0 || m_pProbLists[node] == 0 || pcats == 0 || paridx < 0)
		return R_NilValue;

	if(paridx >= m_numParents[node]) {
		pslot = m_pProbLists[node]->find_slot(0, pcats, 0);
		PROTECT(problist = NEW_NUMERIC(m_numCategories[node]));
		pp = NUMERIC_POINTER(problist);
		for(j = 0; j < m_numCategories[node]; j++) 
			pp[j] = pslot[j];
		return problist;
	}

	npar = m_parents[node][paridx];
	PROTECT(problist = allocVector(VECSXP, m_numCategories[npar]));
	for(j = 0; j < m_numCategories[npar]; j++) {
		pcats[paridx] = j;
		SET_VECTOR_ELT(problist, j, genProbList(node, paridx + 1, pcats));
		UNPROTECT(1);
	}

	return problist;
}


SEXP RPoissonMix::genRmixnet(const char * objectName = (const char*)"mgNetwork") {

	char str[256];
	int node, i, *pslotcats, *pn;
	double *pf;
	SEXP plist, ppars, pcats, pnodeprobs, strnames, pint;

	if(strcmp(objectName, "mgNetwork")) 
		return R_NilValue;

	SEXP cnet = PROTECT(NEW_OBJECT(MAKE_CLASS(objectName)));

	PROTECT(plist = allocVector(STRSXP, 1));
	SET_STRING_ELT(plist, 0, mkChar(objectName));
	SET_SLOT(cnet, install("objectName"), plist);
	UNPROTECT(1);

	PROTECT(plist = allocVector(STRSXP, 1));
	SET_STRING_ELT(plist, 0, mkChar("Pois"));
	SET_SLOT(cnet, install("model"), plist);
	UNPROTECT(1);

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = m_numNodes;
	SET_SLOT(cnet, install("numnodes"), pint);
	UNPROTECT(1);

	PROTECT(strnames = allocVector(STRSXP, m_numNodes));
	for(node = 0; node < m_numNodes; node++) {
		if(m_nodeNames && m_nodeNames[node])
			SET_STRING_ELT(strnames, node, mkChar(m_nodeNames[node]));
		else {
			sprintf(str, "N%d", node+1);
			SET_STRING_ELT(strnames, node, mkChar(str));
		}
	}
	SET_SLOT(cnet, install("nodes"), strnames);
 	UNPROTECT(1);

	SET_SLOT(cnet, install("meta"), mkString("mgNetwork object"));

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = m_maxParents;
	SET_SLOT(cnet, install("maxParents"), pint);
	UNPROTECT(1);

	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	for(node = 0; node < m_numNodes; node++) {
		if(m_numParents[node] <= 0) {
			SET_VECTOR_ELT(plist, node, R_NilValue);
			continue;
		}
		PROTECT(ppars = NEW_INTEGER(m_numParents[node]));
		pn = INTEGER_POINTER(ppars);
		for(i = 0; i < m_numParents[node]; i++)
			// remember to increase the index by 1
			pn[i] = m_parents[node][i] + 1;
		SET_VECTOR_ELT(plist, node, ppars);
		UNPROTECT(1);
	}
	SET_SLOT(cnet, install("parents"), plist);
	UNPROTECT(1);

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = m_maxCategories;
	SET_SLOT(cnet, install("maxCategories"), pint);
	UNPROTECT(1);

	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	for(node = 0; node < m_numNodes; node++) {
		if(m_numCategories[node] > m_maxCategories)
			break;
		PROTECT(pcats = allocVector(STRSXP, m_numCategories[node]));
		for(i = 0; i < m_numCategories[node]; i++) {
			if(m_catIndices && m_catIndices[node])
				sprintf(str, "%d", m_catIndices[node][i]);
			else
				sprintf(str, "C%d", i+1);
			SET_STRING_ELT(pcats, i, mkChar(str));
		}
		SET_VECTOR_ELT(plist, node, pcats);
		UNPROTECT(1);
	}
	SET_SLOT(cnet, install("categories"), plist);
	UNPROTECT(1);

	pslotcats = (int*)CATNET_MALLOC(m_maxParents*sizeof(int));
	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	if(!strcmp(objectName, "mgNetworkC")) {
		for(node = 0; node < m_numNodes; node++) {
			if(!m_pProbLists[node])
				continue;
			PROTECT(pnodeprobs = NEW_NUMERIC(m_pProbLists[node]->nProbSize));
			pf = NUMERIC_POINTER(pnodeprobs);
			memcpy(pf, m_pProbLists[node]->pProbs, m_pProbLists[node]->nProbSize * sizeof(double));
			SET_VECTOR_ELT(plist, node, pnodeprobs);
			UNPROTECT(1);
		}
	}
	else { 
		for(node = 0; node < m_numNodes; node++) {
			memset(pslotcats, 0, m_maxParents*sizeof(int));
			pnodeprobs = genProbList(node, 0, pslotcats);
			SET_VECTOR_ELT(plist, node, pnodeprobs);
			if(pnodeprobs != R_NilValue)
				UNPROTECT(1);
		}
	}
	SET_SLOT(cnet, install("probabilities"), plist);
	UNPROTECT(1);
	CATNET_FREE(pslotcats);

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = complexity();
	SET_SLOT(cnet, install("complexity"), pint);
	UNPROTECT(1);
	
	PROTECT(pint = NEW_NUMERIC(1));
	NUMERIC_POINTER(pint)[0] = loglik();
	SET_SLOT(cnet, install("likelihood"), pint);
	UNPROTECT(1);

	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	for(node = 0; node < m_numNodes; node++) {
		if(m_numCategories[node] <= 0) {
			SET_VECTOR_ELT(plist, node, R_NilValue);
			continue;
		}
		PROTECT(ppars = NEW_NUMERIC(m_numCategories[node]));
		pf = NUMERIC_POINTER(ppars);
		for(i = 0; i < m_numCategories[node]; i++)
			pf[i] = m_lambdas[node][i];
		SET_VECTOR_ELT(plist, node, ppars);
		UNPROTECT(1);
	}
	SET_SLOT(cnet, install("betas"), plist);
	UNPROTECT(1);

	PROTECT(pnodeprobs = NEW_NUMERIC(m_numNodes));
	pf = NUMERIC_POINTER(pnodeprobs);
	memset(pf, 0, m_numNodes * sizeof(double));
	SET_SLOT(cnet, install("sigmas"), pnodeprobs);
	UNPROTECT(1);

	UNPROTECT(1); // cnet
	return cnet;
}

