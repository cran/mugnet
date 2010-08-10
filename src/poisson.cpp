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
 * poisson.cpp
 *
 *  Created on: May 6, 2010
 *      Author: Nikolay Balov
 */


#include "utils.h"
#include "poisson.h"

void CPoissonNet::setNodeBetas(int nnode, double *pbetas) {
	if(nnode < 0 || nnode >= m_numNodes)
		return;
	if(!m_lambdas) {
		m_lambdas = (double**)CATNET_MALLOC(m_numNodes*sizeof(double*));
		memset(m_lambdas, 0, m_numNodes*sizeof(double*));
	}
	if(!m_loglambdas) {
		m_loglambdas = (double**)CATNET_MALLOC(m_numNodes*sizeof(double*));
		memset(m_loglambdas, 0, m_numNodes*sizeof(double*));
	}
	if(!m_lambdas[nnode])
		m_lambdas[nnode] = (double*)CATNET_MALLOC(m_numCategories[nnode]*sizeof(double));
	if(!m_loglambdas[nnode])
		m_loglambdas[nnode] = (double*)CATNET_MALLOC(m_numCategories[nnode]*sizeof(double));
	memcpy(m_lambdas[nnode], pbetas, m_numCategories[nnode]*sizeof(double));
	for(int c = 0; c < m_numCategories[nnode]; c++) {
		if(m_lambdas[nnode][c] > 0)
			m_loglambdas[nnode][c] = log((double)m_lambdas[nnode][c]);
		else
			m_loglambdas[nnode][c] = 0;
	}
}

const double **CPoissonNet::setBetas(double **pbetas, int nbetas) {
	int i;
	if(nbetas != m_numNodes)
		return 0;
	if(m_lambdas) { 
		for (i = 0; i < m_numNodes; i++) {
			if(m_lambdas[i])
				CATNET_FREE(m_lambdas[i]);
			m_lambdas[i] = 0;
		}
		CATNET_FREE(m_lambdas);
		m_lambdas = 0;
	}
	m_lambdas = (double**)CATNET_MALLOC(m_numNodes*sizeof(double*));
	if(m_loglambdas) { 
		for (i = 0; i < m_numNodes; i++) {
			if(m_loglambdas[i])
				CATNET_FREE(m_loglambdas[i]);
			m_loglambdas[i] = 0;
		}
		CATNET_FREE(m_loglambdas);
		m_loglambdas = 0;
	}
	m_loglambdas = (double**)CATNET_MALLOC(m_numNodes*sizeof(double*));
	for (i = 0; i < m_numNodes; i++) {
		m_lambdas[i] = 0;
		if(m_numCategories[i] > 0) {
			m_lambdas[i] = (double*)CATNET_MALLOC(m_numCategories[i]*sizeof(double));
			m_loglambdas[i] = (double*)CATNET_MALLOC(m_numCategories[i]*sizeof(double));
		}
		else
			continue;
		if(pbetas[i])
			memcpy(m_lambdas[i], pbetas[i], m_numCategories[i]*sizeof(double));
		else
			memset(m_lambdas[i], 0, m_numCategories[i]*sizeof(double));
		for(int c = 0; c < m_numCategories[i]; c++) {
			if(m_lambdas[i][c] > 0) 
				m_loglambdas[i][c] = log((double)m_lambdas[i][c]);
			else
				m_loglambdas[i][c] = 0;
		}
	}
	return (const double **)m_lambdas;
}

void CPoissonNet::setNodeSigma(int nnode, double fSigma) {
}

const double *CPoissonNet::setSigmas(double *psigmas, int nsigmas) {
	return (const double *)NULL;
}

const double **CPoissonNet::betas() {
	return (const double **)m_lambdas;
}

const double *CPoissonNet::sigmas() {
	return (const double *)NULL;
}

const double *CPoissonNet::pc() {
	return (const double*)m_pc;
}

const double *CPoissonNet::pcC() {
	return (const double*)m_pcC;
}

const double *CPoissonNet::qc() {
	return (const double*)m_qc;
}

const double *CPoissonNet::qcC() {
	return (const double*)m_qcC;
}

int CPoissonNet::pc_size() {
	if(m_curNode >= 0 && m_curNode < m_numNodes && m_numCategories)
		return m_numCategories[m_curNode];
	return 0;
}

int CPoissonNet::pcC_size() {
	if(m_curNode >= 0 && m_curNode < m_numNodes && m_parCatSetSize > 0)
		return m_parCatSetSize;
	return 0;
}


void CPoissonNet::_release() {
	int i;

	if(m_pCatnetSamples)
		CATNET_FREE(m_pCatnetSamples);
	m_pCatnetSamples = 0; 
	m_nCatnetSamples = 0;

	if(m_lambdas) { 
		for (i = 0; i < m_numNodes; i++) {
			if(m_lambdas[i])
				CATNET_FREE(m_lambdas[i]);
			m_lambdas[i] = 0;
		}
		CATNET_FREE(m_lambdas);
		m_lambdas = 0;
	}
	if(m_loglambdas) { 
		for (i = 0; i < m_numNodes; i++) {
			if(m_loglambdas[i])
				CATNET_FREE(m_loglambdas[i]);
			m_loglambdas[i] = 0;
		}
		CATNET_FREE(m_loglambdas);
		m_loglambdas = 0;
	}
	if(m_nodeLoglik)
		CATNET_FREE(m_nodeLoglik);
	m_nodeLoglik = 0;

	release_pqcC();
	// First release CPoissonNet's resources then its parents CATNET
	CATNET<char, MAX_NODE_NAME, double>::_release();
}

void CPoissonNet::release_pqcC() {
	if(m_pc)
		CATNET_FREE(m_pc);
	m_pc = 0;
	if(m_pcC)
		CATNET_FREE(m_pcC);
	m_pcC = 0;
	if(m_qc)
		CATNET_FREE(m_qc);
	m_qc = 0;
	if(m_qcC)
		CATNET_FREE(m_qcC);
	m_qcC = 0;
}

void CPoissonNet::_reset() {
	CATNET<char, MAX_NODE_NAME, double>::_reset();
	m_lambdas = 0;
	m_loglambdas = 0;
	m_curNode =0;
	m_pc = m_pcC = 0;
	m_qc = m_qcC = 0;
	m_parCatSetSize = 0;
	m_nodeLoglik = 0;
	m_pCatnetSamples = 0; 
	m_nCatnetSamples = 0;
}

CPoissonNet& CPoissonNet::operator =(const CPoissonNet &cnet) {

	CATNET<char, MAX_NODE_NAME, double> ::init(
			(int)cnet.m_numNodes, (int)cnet.m_maxParents, (int)cnet.m_maxCategories,
			(const char**)cnet.m_nodeNames, (const int*)cnet.m_numParents, 
			(const int**)cnet.m_parents,
			(const int*)cnet.m_numCategories, (const PROB_LIST<double> **)cnet.m_pProbLists, 
			(int)cnet.m_complexity, (double)cnet.m_loglik);

	if(m_numNodes <= 0)
		return *this;

	int i;
	m_lambdas = (double**)CATNET_MALLOC(m_numNodes*sizeof(double*));
	m_loglambdas = (double**)CATNET_MALLOC(m_numNodes*sizeof(double*));
	m_nodeLoglik = (double*)CATNET_MALLOC(m_numNodes*sizeof(double));
	for (i = 0; i < m_numNodes; i++) {
		m_lambdas[i] = 0;
		if(m_numCategories[i] > 0) {
			m_lambdas[i] = (double*)CATNET_MALLOC(m_numCategories[i]*sizeof(double));
			m_loglambdas[i] = (double*)CATNET_MALLOC(m_numCategories[i]*sizeof(double));
		}
		else
			continue;
		if(cnet.m_lambdas && cnet.m_lambdas[i])
			memcpy(m_lambdas[i], cnet.m_lambdas[i], m_numCategories[i]*sizeof(double));
		else
			memset(m_lambdas[i], 0, m_numCategories[i]*sizeof(double));
		for(int c = 0; c < m_numCategories[i]; c++) {
			if(m_lambdas[i][c] > 0)
				m_loglambdas[i][c] = log((double)m_lambdas[i][c]);
			else
				m_loglambdas[i][c] = 0;
		}
	}
	if(cnet.m_nodeLoglik)
		memcpy(m_nodeLoglik, cnet.m_nodeLoglik, m_numNodes*sizeof(double));
	else
		memset(m_nodeLoglik, 0, m_numNodes*sizeof(double));

	return *this;
}

double CPoissonNet::setNodeLoglik(int nnode, double lik) {
	m_nodeLoglik[nnode] = lik;
	return m_nodeLoglik[nnode];
}

double CPoissonNet::getNodeLoglik(int nnode) {
	return m_nodeLoglik[nnode];
}
	
double CPoissonNet::loglik() {
	int i;
	if(!m_nodeLoglik)
		return 0;
	m_loglik = 0;
	for(i = 0; i < m_numNodes; i++)
		m_loglik += m_nodeLoglik[i];
	return m_loglik;
}

double CPoissonNet::estimateParameters(int nnode, double *psamples, int nsamples, 
					double *pBetas, double *pSigma) {
	double ftemp, fsum;
	int ic, j, c_setSize;
	double *plogbeta;

	if(m_curNode != nnode)
		return -FLT_MAX;

	c_setSize = m_numCategories[nnode];
	plogbeta = (double*)CATNET_MALLOC(c_setSize*sizeof(double));
	if(!plogbeta)
		return -FLT_MAX;

	for(ic = 0; ic < c_setSize; ic++) {
		ftemp = 0;
		fsum = 0;
		for(j = 0; j < nsamples; j++) {
			ftemp += m_qc[j*c_setSize + ic] * psamples[j * m_numNodes + nnode]; 
			fsum += m_qc[j*c_setSize + ic];
		}
		if(fsum > 0) {
			/* update only if possible, don't put 0, remember the last good value instead */
			pBetas[ic] = ftemp / fsum;
		}
	}

	for(ic = 0; ic < c_setSize; ic++) {
		if(pBetas[ic] > 0)
			plogbeta[ic] = log(pBetas[ic]);
		else
			plogbeta[ic] = -FLT_MAX;
	}

	fsum = 0;
	for(j = 0; j < nsamples; j++) {
		for(ic = 0; ic < c_setSize; ic++) {		
			ftemp = psamples[j * m_numNodes + nnode]*plogbeta[ic] - pBetas[ic];
			fsum += m_qc[j*c_setSize + ic] * ftemp;
		}
	}

	CATNET_FREE(plogbeta);

	*pSigma = 0;

	return fsum;
}

int CPoissonNet::findNodeMarginalProb(int nnode, double *psamples, int nsamples) {

	int res, j, ic, numcats;
	double fsum, fmax, *paux;

	if(nnode < 0 || nnode >= m_numNodes || !m_lambdas || !psamples || nsamples < 1)
		return ERR_CATNET_PARAM;

	if(m_pc) 
		CATNET_FREE(m_pc);
	m_pc = 0;
	if(m_qc) 
		CATNET_FREE(m_qc);
	m_qc = 0;

	m_curNode = nnode;

	numcats = m_numCategories[nnode];

	m_pc = (double*)CATNET_MALLOC(numcats*sizeof(double));
	if(!m_pc)
		return ERR_CATNET_MEM;
	memset(m_pc, 0, numcats*sizeof(double));

	if(!m_pCatnetSamples || m_nCatnetSamples < 1) {
		res = marginalProb(nnode);
		if(res != ERR_CATNET_OK)
			return res;
		for(ic = 0; ic < numcats; ic++)
			m_pc[ic] = getCatProb(ic);
	}
	else {
		memset(m_pc, 0, numcats*sizeof(double));
		for(j = 0; j < m_nCatnetSamples; j++) {
			ic = m_pCatnetSamples[j * m_numNodes + nnode];
			m_pc[ic] += 1;
		}
		fsum = 1/(double)m_nCatnetSamples;
		for(ic = 0; ic < numcats; ic++)
			m_pc[ic] *= fsum;
	}

	m_qc = (double*)CATNET_MALLOC(numcats*nsamples*sizeof(double));
	if(!m_qc)
		return ERR_CATNET_MEM;
	memset(m_qc, 0, numcats*nsamples*sizeof(double));

	paux = (double*)CATNET_MALLOC(numcats*sizeof(double));
	for (j = 0; j < nsamples; j++) {
		fmax = -FLT_MAX;
		for(ic = 0; ic < numcats; ic++) {
			paux[ic] = psamples[j * m_numNodes + nnode] * m_loglambdas[nnode][ic] - m_lambdas[nnode][ic];
			if(fmax < paux[ic])
				fmax = paux[ic];
		}
		fsum = 0;
		for(ic = 0; ic < numcats; ic++) {
			paux[ic] = m_pc[ic]*exp(paux[ic] - fmax); 
			fsum += paux[ic];
		}
		if(fsum > 0) 
			fsum = 1/fsum;
		for(ic = 0; ic < numcats; ic++) {
			m_qc[j*numcats + ic] = fsum*paux[ic];
		}
	}
	CATNET_FREE(paux);

	return ERR_CATNET_OK;
}


int CPoissonNet::findNodeJointProb(int nnode, int *parnodes, int numpars, double *psamples, int nsamples) {

	int res, i, j, ic, icc;
	int *pnodes, numnodes, *pcats, numcats;
	double fsum, fdiffsum, fmax, *paux;

	if(nnode < 0 || nnode >= m_numNodes || !m_lambdas || 
		!psamples || nsamples < 1)
		return ERR_CATNET_PARAM;
		
	if(!parnodes || !numpars)
		return ERR_CATNET_OK;
		
	// pnodes has all parent nodes plus nnode attached at the end
	pnodes = (int*)CATNET_MALLOC((numpars+1)*sizeof(int));
	memcpy(pnodes, parnodes, numpars*sizeof(int));
	pnodes[numpars] = nnode;
	numnodes = numpars + 1;

	if(m_pcC) 
		CATNET_FREE(m_pcC);
	m_pcC = 0;
	if(m_qcC) 
		CATNET_FREE(m_qcC);
	m_qcC = 0;

	m_curNode = nnode;

	numcats = m_numCategories[nnode];

	// at this point m_nBlockSizes == numnodes and m_pBlockSizes are filled
	m_parCatSetSize = 1;
	for(i = 0; i < numnodes; i++) 
		m_parCatSetSize *= m_numCategories[pnodes[i]];

	paux = (double*)CATNET_MALLOC(m_parCatSetSize*sizeof(double));
	pcats = (int*)CATNET_MALLOC(m_parCatSetSize*numnodes*sizeof(int));
	m_pcC = (double*)CATNET_MALLOC(m_parCatSetSize*sizeof(double));
	m_qcC = (double*)CATNET_MALLOC(m_parCatSetSize*nsamples*sizeof(double));
	if(!paux || !m_pcC || !m_qcC)
		return ERR_CATNET_MEM;

	memset(m_pcC, 0, m_parCatSetSize*sizeof(double));
	memset(m_qcC, 0, m_parCatSetSize*nsamples*sizeof(double));

	if(!m_pCatnetSamples || m_nCatnetSamples < 1) {
		// find the exact joint probability for small sets
		res = marginalProb(pnodes, numnodes);
		if(res != ERR_CATNET_OK) {
			CATNET_FREE(pnodes);
			return res;
		}

		if(m_margProbSize != m_parCatSetSize) {
			// fatal  error
			error("m_margProbSize != m_parCatSetSize");
			return ERR_CATNET_PROC;
		}

		for(ic = 0; ic < m_parCatSetSize; ic++) {
			// determine the category indices
			fsum = ic;
			for(i = 0; i < numnodes; i++) {
				pcats[ic * numnodes + i] = (int)(fsum/m_pBlockSizes[i]);
				fsum -= pcats[ic * numnodes + i] * m_pBlockSizes[i];
			}
			m_pcC[ic] = getCatProb(pcats + ic * numnodes, numnodes);
		}
	}
	else {
		if(!m_pBlockSizes || m_nBlockSizes < numpars) {
			if(m_pBlockSizes)
				CATNET_FREE(m_pBlockSizes);
			m_nBlockSizes = numpars;
			m_pBlockSizes = (int*) CATNET_MALLOC(m_nBlockSizes * sizeof(int));
		}

		m_pBlockSizes[numnodes - 1] = 1;
		for (i = numnodes - 2; i >= 0; i--)
			m_pBlockSizes[i] = m_pBlockSizes[i + 1] * m_numCategories[pnodes[i + 1]];

		memset(m_pcC, 0, m_parCatSetSize*sizeof(double));

		for(j = 0; j < m_nCatnetSamples; j++) {
			ic = 0;
			for (i = 0; i < numnodes; i++)
				ic += (m_pBlockSizes[i] * m_pCatnetSamples[j * m_numNodes + pnodes[i]]);
			m_pcC[ic] += 1;
		}

		fsum = 1/(double)m_nCatnetSamples;
		for(ic = 0; ic < m_parCatSetSize; ic++)
			m_pcC[ic] *= fsum;

		for(ic = 0; ic < m_parCatSetSize; ic++) {
			// determine the category indices
			fsum = ic;
			for(i = 0; i < numnodes; i++) {
				pcats[ic * numnodes + i] = (int)(fsum/m_pBlockSizes[i]);
				fsum -= pcats[ic * numnodes + i] * m_pBlockSizes[i];
			}
		}
	}

	for (j = 0; j < nsamples; j++) {
		fmax = -FLT_MAX;
		for(ic = 0; ic < m_parCatSetSize; ic++) {
			fdiffsum = 0;
			for(i = 0; i < numnodes; i++) {
				fdiffsum += (psamples[j * m_numNodes + pnodes[i]] * m_loglambdas[pnodes[i]][pcats[ic * numnodes + i]] - m_lambdas[pnodes[i]][pcats[ic * numnodes + i]]);
			}
			paux[ic] = fdiffsum;
			if(fmax < paux[ic])
				fmax = paux[ic];	
		}

		fsum = 0;
		for(ic = 0; ic < m_parCatSetSize; ic++) {
			paux[ic] = m_pcC[ic]*exp(paux[ic] - fmax);
			fsum += paux[ic];
		}
		if(fsum > 0) {
			fsum = 1/fsum;
			for(ic = 0; ic < m_parCatSetSize; ic++) 
				m_qcC[j*m_parCatSetSize + ic] = fsum*paux[ic];
		}
	}

	CATNET_FREE(pnodes);
	CATNET_FREE(pcats);
	CATNET_FREE(paux);

	if(m_pc) 
		CATNET_FREE(m_pc);
	if(m_qc) 
		CATNET_FREE(m_qc);		
	m_pc = (double*)CATNET_MALLOC(numcats*sizeof(double));
	m_qc = (double*)CATNET_MALLOC(numcats*nsamples*sizeof(double));

	for(icc = 0; icc < numcats; icc++) {
		fsum = 0;
		for(ic = 0; ic < m_parCatSetSize; ic += numcats) 
			fsum += m_pcC[ic + icc];
		m_pc[icc] = fsum;
		for (j = 0; j < nsamples; j++) {
			fsum = 0;
			for(ic = 0; ic < m_parCatSetSize; ic += numcats) 
				fsum += m_qcC[j*m_parCatSetSize + ic + icc];
			m_qc[j*numcats + icc] = fsum;
		}
	}

	return ERR_CATNET_OK;
}

int CPoissonNet::findParentsJointProb(int *parnodes, int numpars, double *psamples, int nsamples) {

	int res, i, j, ic, k;
	int *pcats;
	double fdiffsum, fsum, *paux;

	if(!m_lambdas || !psamples || nsamples < 1)
		return ERR_CATNET_PARAM;
	if(!parnodes || !numpars)
		return ERR_CATNET_PARAM;

	if(m_pcC) CATNET_FREE(m_pcC);
	m_pcC = 0;
	if(m_qcC) CATNET_FREE(m_qcC);
	m_qcC = 0;

	// at this point m_nBlockSizes == numnodes and m_pBlockSizes are filled
	m_parCatSetSize = 1;
	for(i = 0; i < numpars; i++)
		m_parCatSetSize *= m_numCategories[parnodes[i]];

	m_pcC = (double*)CATNET_MALLOC(m_parCatSetSize*sizeof(double));
	pcats = (int*)CATNET_MALLOC(m_parCatSetSize*numpars*sizeof(int));
	if(!m_pcC || !pcats)
		return ERR_CATNET_MEM;
	memset(m_pcC, 0, m_parCatSetSize*sizeof(double));

	if(!m_pCatnetSamples || m_nCatnetSamples < 1) {

		res = marginalProb(parnodes, numpars);
		if(res != ERR_CATNET_OK) {
			CATNET_FREE(pcats);
			return res;
		}

		if(m_margProbSize != m_parCatSetSize) {
			CATNET_FREE(pcats);
			return ERR_CATNET_PROC;
		}

		for(ic = 0; ic < m_parCatSetSize; ic++) {
			// determine the category indices
			fsum = ic;
			for(i = 0; i < numpars; i++) {
				pcats[ic * numpars + i] = (int)(fsum/m_pBlockSizes[i]);
				fsum -= pcats[ic * numpars + i] * m_pBlockSizes[i];
			}
			m_pcC[ic] = getCatProb(pcats + ic * numpars, numpars);
		}
	}
	else {
		if(!m_pBlockSizes || m_nBlockSizes < numpars) {
			if(m_pBlockSizes)
				CATNET_FREE(m_pBlockSizes);
			m_nBlockSizes = numpars;
			m_pBlockSizes = (int*) CATNET_MALLOC(m_nBlockSizes * sizeof(int));
		}

		m_pBlockSizes[numpars - 1] = 1;
		for (i = numpars - 2; i >= 0; i--) {
			m_pBlockSizes[i] = m_pBlockSizes[i + 1] * m_numCategories[parnodes[i + 1]];
			
		}

		for(j = 0; j < m_nCatnetSamples; j++) {
			ic = 0;
			for (i = 0; i < numpars; i++)
				ic += (m_pBlockSizes[i] * m_pCatnetSamples[j * m_numNodes + parnodes[i]]);
			m_pcC[ic] += 1;
		}

		fsum = 1/(double)m_nCatnetSamples;
		for(ic = 0; ic < m_parCatSetSize; ic++)
			m_pcC[ic] *= fsum;

		for(ic = 0; ic < m_parCatSetSize; ic++) {
			// determine the category indices
			fsum = ic;
			for(i = 0; i < numpars; i++) {
				pcats[ic * numpars + i] = (int)(fsum/m_pBlockSizes[i]);
				fsum -= pcats[ic * numpars + i] * m_pBlockSizes[i];
			}
		}
	}

	paux = (double*)CATNET_MALLOC(m_parCatSetSize*sizeof(double));
	m_qcC = (double*)CATNET_MALLOC(m_parCatSetSize*nsamples*sizeof(double));
	if(!m_qcC || !paux)
		return ERR_CATNET_MEM;

	for (j = 0; j < nsamples; j++) {

		fsum = 0;
		for(ic = 0; ic < m_parCatSetSize; ic++) {
			// determine the category indices
			k = ic;
			for(i = 0; i < numpars; i++) {
				pcats[i] = (int)(k/m_pBlockSizes[i]);
				k -= pcats[i] * m_pBlockSizes[i];
			}
			fdiffsum = 0;
			for(i = 0; i < numpars; i++) {
				fdiffsum += (psamples[j * m_numNodes + parnodes[i]] * m_loglambdas[parnodes[i]][pcats[i]] - m_lambdas[parnodes[i]][pcats[i]]);
			}
			paux[ic] = m_pcC[ic]*exp(fdiffsum);
			fsum += paux[ic];
		}
		if(fsum > 0) fsum = 1/fsum;
		for(ic = 0; ic < m_parCatSetSize; ic++)
			m_qcC[j*m_parCatSetSize + ic] = fsum*paux[ic];
	}

	CATNET_FREE(pcats);
	CATNET_FREE(paux);

	return ERR_CATNET_OK;
}

int* CPoissonNet::catnetSample(int nsamples) {

	int *porder;
	int i, j, k, nnode, *pnodepars, *pnodesample;
	double u, v, *pnodeprob;
	PROB_LIST<double>* pProbList;

	if(nsamples < 1 || m_numNodes < 1)
		return 0;
	porder = getOrder();
	if(!porder) 
		return 0;

	if(m_pCatnetSamples)
		CATNET_FREE(m_pCatnetSamples);
	m_nCatnetSamples = nsamples;
	m_pCatnetSamples = (int*)CATNET_MALLOC(m_nCatnetSamples * m_numNodes * sizeof(int)); 
	
	pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));

	for(k = 0; k < m_numNodes; k++) {
		nnode = porder[k];
		pnodepars = m_parents[nnode];
		pProbList = (PROB_LIST<double>*)getNodeProb(nnode);

		for (j = 0; j < m_nCatnetSamples; j++) {
			for (i = 0; i < m_numParents[nnode]; i++) {
				if (pnodepars[i] < 0 || pnodepars[i] >= m_numNodes)
					break;
				pnodesample[i] = (int)m_pCatnetSamples[j * m_numNodes + pnodepars[i]];
			}
			pnodeprob = pProbList->find_slot(0, pnodesample, 0);

			u = (double)rand() / (double)RAND_MAX;
			v = 0;
			for(i = 0; i < m_numCategories[nnode]; i++) {
				v += pnodeprob[i];
				if(u <= v)
					break;
			}
			m_pCatnetSamples[j * m_numNodes + nnode] = (double)i;
		}
	}

	CATNET_FREE(pnodesample);
	CATNET_FREE(porder);

	return m_pCatnetSamples;
}

int CPoissonNet::sample(double *psamples, int nsamples) {

	int *porder;
	int i, j, k, nnode, *pnodepars, *pnodesample;
	double u, v, *pnodeprob;
	PROB_LIST<double>* pProbList; 

	if(!m_lambdas)
		return ERR_CATNET_INIT;
	
	if(!psamples || nsamples < 1)
		return ERR_CATNET_PARAM;	

	porder = getOrder();
	if(!porder) {
		return ERR_CATNET_PROC;
	}

	pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));

	for(k = 0; k < m_numNodes; k++) {
		nnode = porder[k];
		pnodepars = m_parents[nnode];
		pProbList = (PROB_LIST<double>*)getNodeProb(nnode);

		for (j = 0; j < nsamples; j++) {
			for (i = 0; i < m_numParents[nnode]; i++) {
				if (pnodepars[i] < 0 || pnodepars[i] >= m_numNodes)
					break;
				pnodesample[i] = (int)psamples[j * m_numNodes + pnodepars[i]];
			}
			pnodeprob = pProbList->find_slot(0, pnodesample, 0);

			u = (double)rand() / (double)RAND_MAX;
			v = 0;
			for(i = 0; i < m_numCategories[nnode]; i++) {
				v += pnodeprob[i];
				if(u <= v)
					break;
			}
			psamples[j * m_numNodes + nnode] = (double)i;
		}
	}

	CATNET_FREE(pnodesample);

	for(k = 0; k < m_numNodes; k++) {
		nnode = porder[k];
		for (j = 0; j < nsamples; j++) {
			i = (int)psamples[j * m_numNodes + nnode];
			v = exp(-m_lambdas[nnode][i]);
			u = 1;
			i = 0;
			while(u > v) {
				i++;
				u *= ((double)rand() / (double)RAND_MAX);
			}
			psamples[j * m_numNodes + nnode] = (double)(i - 1);
		}
	}

	CATNET_FREE(porder);

	return ERR_CATNET_OK;
}

int CPoissonNet::predict(double *psamples, int nsamples) {

	int *porder;
	int j, k, ic, icC, nnode, imax, cC_setSize;
	double fsum, faux, *pnodeprob;
	PROB_LIST<double>* pProbList;

	if(!m_lambdas)
		return ERR_CATNET_INIT;
	
	if(!psamples || nsamples < 1)
		return ERR_CATNET_PARAM;	

	porder = getOrder();
	if(!porder)
		return ERR_CATNET_PROC;

	for(k = 0; k < m_numNodes; k++) {
		nnode = porder[k];
		pProbList = (PROB_LIST<double>*)getNodeProb(nnode);
		pnodeprob = pProbList->pProbs;

		if(m_numParents[nnode] <= 0) {
			for (j = 0; j < nsamples; j++) {
				// predict only NA ones
				if(!isnan(psamples[j*m_numNodes + nnode]))
					continue;
				fsum = 0;
				for(ic = 0; ic < m_numCategories[nnode]; ic++) {
					fsum += (m_lambdas[nnode][ic] * pnodeprob[ic]);
				}
				psamples[j*m_numNodes + nnode] = fsum;
			}
			continue;
		}
 
		if(findParentsJointProb(m_parents[nnode], m_numParents[nnode], psamples, nsamples) != ERR_CATNET_OK)
			continue;
		cC_setSize = pcC_size();
		// now m_qcC contains P(x_j, j\in J_i | y_j^s, j\in J_i), j = 1,...,numsamples
		for (j = 0; j < nsamples; j++) {
			// predict only NA ones
			if(!isnan(psamples[j*m_numNodes + nnode]))
				continue;
if(1) {
			fsum = 0;
			for(ic = 0; ic < m_numCategories[nnode]; ic++) {
				faux = 0;
				for(icC = 0; icC < cC_setSize; icC++) {
					faux += pnodeprob[icC*m_numCategories[nnode] + ic] * 
						m_qcC[j*cC_setSize + icC];
				}
				fsum += m_lambdas[nnode][ic] * faux;
			}
			psamples[j*m_numNodes + nnode] = fsum;
}
else {
			imax = 0;
			fsum = -FLT_MAX;
			for(ic = 0; ic < m_numCategories[nnode]; ic++) {
				faux = 0;
				for(icC = 0; icC < cC_setSize; icC++) {
					faux += pnodeprob[ic*m_numCategories[nnode] + ic] * 
						m_qcC[j*cC_setSize + icC];
				}
				if(fsum < faux) {
					fsum = faux;
					imax = ic;
				}
			}
			psamples[j*m_numNodes + nnode] = m_lambdas[nnode][imax];
}
		}
	}

	CATNET_FREE(porder);

	return ERR_CATNET_OK;
}

double CPoissonNet::findLogNodeLikelihood(int nnode, double *psamples, int nsamples) {
	
	double faux, fsum, fsumsum, factoriel, fLogLik;
	int j, cC_setSize, ic, icC;
	double *pnodeprob;
	PROB_LIST<double>* pProbList;

	if(!m_lambdas || nnode < 0 || nnode >= m_numNodes || m_lambdas[nnode] <= 0)
		return -FLT_MAX;
	
	if(!psamples || nsamples < 1)
		return -FLT_MAX;

	pProbList = (PROB_LIST<double>*)getNodeProb(nnode);
	pnodeprob = pProbList->pProbs;

	cC_setSize = 0;
	if(m_numParents[nnode] > 0) {
		if(findParentsJointProb(m_parents[nnode], m_numParents[nnode], psamples, nsamples) != ERR_CATNET_OK)
			return -FLT_MAX;
		// now m_qcC contains P(x_j, j\in J_i | y_j^s, j\in J_i), j = 1,...,numsamples
		cC_setSize = pcC_size();
	}

	fLogLik = 0;
	for(j = 0; j < nsamples; j++) {

		faux = (int)psamples[j * m_numNodes + nnode];
		factoriel = 0;
		while(faux > 0)
			factoriel += log(faux--);

		fsumsum = 0;
		for(ic = 0; ic < m_numCategories[nnode]; ic++) {
			faux = psamples[j * m_numNodes + nnode]*m_loglambdas[nnode][ic] - m_lambdas[nnode][ic];
			if(m_numParents[nnode] > 0) {
				fsum = 0;
				for(icC = 0; icC < cC_setSize; icC++) {
					fsum += pnodeprob[icC*m_numCategories[nnode] + ic] * 
						m_qcC[j*cC_setSize + icC];
					}
			}
			else
				fsum = pnodeprob[ic];
			fsumsum += exp(faux - factoriel) * fsum;
		}
		fLogLik += log(fsumsum);
	}
	return fLogLik;
}

double CPoissonNet::findNodeLoglik(int nnode, double *psamples, int nsamples) {
	
	double faux, fLogLik;
	int j, cC_setSize, ic, icC, iC;
	double *pnodeprob;
	PROB_LIST<double>* pProbList;
	
	if(!m_lambdas || nnode < 0 || nnode >= m_numNodes || m_lambdas[nnode] <= 0)
		return -FLT_MAX;
	
	if(!psamples || nsamples < 1)
		return -FLT_MAX;

	pProbList = (PROB_LIST<double>*)getNodeProb(nnode);
	pnodeprob = pProbList->pProbs;

	cC_setSize = 0;
	if(m_numParents[nnode] > 0) {
		if(findNodeJointProb(nnode, m_parents[nnode], m_numParents[nnode], psamples, nsamples) != ERR_CATNET_OK)
			return -FLT_MAX;
		// now m_qcC contains P(x, i, x_j, j\in J_i | y_i, y_j^s, j\in J_i), j = 1,...,numsamples
		cC_setSize = pcC_size();
	}
	else {
		findNodeMarginalProb(nnode, psamples, nsamples);
	}
	
	fLogLik = 0;
	if(m_numParents[nnode] > 0) {
		for(j = 0; j < nsamples; j++) {
			iC = 0;
			for(icC = 0; icC < cC_setSize; icC += m_numCategories[nnode]) {
				for(ic = 0; ic < m_numCategories[nnode]; ic++) {
					faux = psamples[j * m_numNodes + nnode]*m_loglambdas[nnode][ic] - m_lambdas[nnode][ic];
					if(pnodeprob[iC*m_numCategories[nnode] + ic] > 0)
						fLogLik += m_qcC[j*cC_setSize + icC + ic] * (faux + log(pnodeprob[iC*m_numCategories[nnode] + ic]));
					else
						fLogLik += m_qcC[j*cC_setSize + icC + ic] * faux;
				}
				iC++;
			}
		}
	}
	else {
		for(j = 0; j < nsamples; j++) {
			for(ic = 0; ic < m_numCategories[nnode]; ic++) {
				faux = psamples[j * m_numNodes + nnode]*m_loglambdas[nnode][ic] - m_lambdas[nnode][ic];
				if(pnodeprob[ic] > 0)
					fLogLik += m_qc[j*m_numCategories[nnode] + ic] * (faux + log(pnodeprob[ic]));
				else
					fLogLik += m_qc[j*m_numCategories[nnode] + ic] * faux;
			}
		}
	}
	return fLogLik;
}
