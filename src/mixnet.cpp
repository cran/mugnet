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
 * mixnet.cpp
 *
 *  Created on: Nov 11, 2010
 *      Author: Nikolay Balov
 */


#include "utils.h"
#include "mixnet.h"

void CMixNet::setNodeBetas(int nnode, double *pbetas) {
	if(nnode < 0 || nnode >= m_numNodes)
		return;
	if(!m_betas) {
		m_betas = (double**)CATNET_MALLOC(m_numNodes*sizeof(double*));
		memset(m_betas, 0, m_numNodes*sizeof(double*));
	}
	if(!m_betas[nnode])
		m_betas[nnode] = (double*)CATNET_MALLOC(m_numCategories[nnode]*sizeof(double));
	memcpy(m_betas[nnode], pbetas, m_numCategories[nnode]*sizeof(double));
}

const double **CMixNet::setBetas(double **pbetas, int nbetas) {
	int i;
	if(nbetas != m_numNodes)
		return 0;
	if(m_betas) { 
		for (i = 0; i < m_numNodes; i++) {
			if(m_betas[i])
				CATNET_FREE(m_betas[i]);
			m_betas[i] = 0;
		}
		CATNET_FREE(m_betas);
		m_betas = 0;
	}
	m_betas = (double**)CATNET_MALLOC(m_numNodes*sizeof(double*));
	for (i = 0; i < m_numNodes; i++) {
		m_betas[i] = 0;
		if(m_numCategories[i] > 0)
			m_betas[i] = (double*)CATNET_MALLOC(m_numCategories[i]*sizeof(double));
		else
			continue;
		if(pbetas[i])
			memcpy(m_betas[i], pbetas[i], m_numCategories[i]*sizeof(double));
		else
			memset(m_betas[i], 0, m_numCategories[i]*sizeof(double));
	}
	return (const double **)m_betas;
}

void CMixNet::setNodeSigma(int nnode, double fSigma) {
	if(nnode < 0 || nnode >= m_numNodes)
		return;
	if(!m_sigmas) 
		m_sigmas = (double*)CATNET_MALLOC(m_numNodes*sizeof(double));
	m_sigmas[nnode] = fSigma;
}

const double *CMixNet::setSigmas(double *psigmas, int nsigmas) {
	if(nsigmas != m_numNodes)
		return 0;
	if(m_sigmas) {
		CATNET_FREE(m_sigmas);
		m_sigmas = 0; 
	}
	m_sigmas = (double*)CATNET_MALLOC(m_numNodes*sizeof(double));
	if(psigmas)
		memcpy(m_sigmas, psigmas, m_numNodes*sizeof(double));
	else
		memset(m_sigmas, 0, m_numNodes*sizeof(double));
	return (const double *)m_sigmas;
}

const double **CMixNet::betas() {
	return (const double **)m_betas;
}

const double *CMixNet::sigmas() {
	return (const double *)m_sigmas;
}

const double *CMixNet::pc() {
	return (const double*)m_pc;
}

const double *CMixNet::pcC() {
	return (const double*)m_pcC;
}

const double *CMixNet::qc() {
	return (const double*)m_qc;
}

const double *CMixNet::qcC() {
	return (const double*)m_qcC;
}

int CMixNet::pc_size() {
	if(m_curNode >= 0 && m_curNode < m_numNodes && m_numCategories)
		return m_numCategories[m_curNode];
	return 0;
}

int CMixNet::pcC_size() {
	if(m_curNode >= 0 && m_curNode < m_numNodes && m_parCatSetSize > 0)
		return m_parCatSetSize;
	return 0;
}


void CMixNet::_release() {
	int i;

	if(m_pCatnetSamples)
		CATNET_FREE(m_pCatnetSamples);
	m_pCatnetSamples = 0; 
	m_nCatnetSamples = 0;

	if(m_betas) { 
		for (i = 0; i < m_numNodes; i++) {
			if(m_betas[i])
				CATNET_FREE(m_betas[i]);
			m_betas[i] = 0;
		}
		CATNET_FREE(m_betas);
		m_betas = 0;
	}
	if(m_sigmas) {
		CATNET_FREE(m_sigmas);
		m_sigmas = 0; 
	}

	if(m_nodeLoglik)
		CATNET_FREE(m_nodeLoglik);
	m_nodeLoglik = 0;

	release_pqcC();
	// First release CMixNet's resources then its parents CATNET
	CATNET<char, double>::_release();
}

void CMixNet::release_pqcC() {
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

void CMixNet::_reset() {
	CATNET<char, double>::_reset();
	m_betas = 0;
	m_sigmas = 0;
	m_curNode =0;
	m_pc = m_pcC = 0;
	m_qc = m_qcC = 0;
	m_parCatSetSize = 0;
	m_nodeLoglik = 0;
	m_pCatnetSamples = 0; 
	m_nCatnetSamples = 0;
}

CMixNet& CMixNet::operator =(const CMixNet &cnet) {

	CATNET<char, double> ::init(
			(int)cnet.m_numNodes, (int)cnet.m_maxParents, (int)cnet.m_maxCategories,
			(const char**)cnet.m_nodeNames, (const int*)cnet.m_numParents, 
			(const int**)cnet.m_parents,
			(const int*)cnet.m_numCategories, (const PROB_LIST<double> **)cnet.m_pProbLists, 
			(int)cnet.m_complexity, (double)cnet.m_loglik);

	if(m_numNodes <= 0)
		return *this;

	int i;
	m_betas = (double**)CATNET_MALLOC(m_numNodes*sizeof(double*));
	m_sigmas = (double*)CATNET_MALLOC(m_numNodes*sizeof(double));
	m_nodeLoglik = (double*)CATNET_MALLOC(m_numNodes*sizeof(double));
	for (i = 0; i < m_numNodes; i++) {
		m_betas[i] = 0;
		if(m_numCategories[i] > 0)
			m_betas[i] = (double*)CATNET_MALLOC(m_numCategories[i]*sizeof(double));
		else
			continue;
		if(cnet.m_betas && cnet.m_betas[i])
			memcpy(m_betas[i], cnet.m_betas[i], m_numCategories[i]*sizeof(double));
		else
			memset(m_betas[i], 0, m_numCategories[i]*sizeof(double));
	}
	if(cnet.m_sigmas)
		memcpy(m_sigmas, cnet.m_sigmas, m_numNodes*sizeof(double));
	else
		memset(m_sigmas, 0, m_numNodes*sizeof(double));
	if(cnet.m_nodeLoglik)
		memcpy(m_nodeLoglik, cnet.m_nodeLoglik, m_numNodes*sizeof(double));
	else
		memset(m_nodeLoglik, 0, m_numNodes*sizeof(double));

	return *this;
}

double CMixNet::setNodeLoglik(int nnode, double lik) {
	m_nodeLoglik[nnode] = lik;
	return m_nodeLoglik[nnode];
}

double CMixNet::getNodeLoglik(int nnode) {
	return m_nodeLoglik[nnode];
}
	
double CMixNet::loglik() {
	int i;
	if(!m_nodeLoglik)
		return 0;
	m_loglik = 0;
	for(i = 0; i < m_numNodes; i++)
		m_loglik += m_nodeLoglik[i];
	return m_loglik;
}

double CMixNet::estimateParameters(int nnode, double *psamples, int nsamples, 
					double *pBetas, double *pSigma) {
	double ftemp, fsum, fBetaSigmaLoglik;
	int ic, j, c_setSize;

	if(m_curNode != nnode)
		return -FLT_MAX;

	c_setSize = m_numCategories[nnode];
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
		//Rprintf("fsum = %f, pBetasNext[%d][%d] = %f\n", fsum, nnode, ic, pBetasNext[nnode][ic]);
	}

	fsum = 0;
	for(ic = 0; ic < c_setSize; ic++) {
		for(j = 0; j < nsamples; j++) {
			ftemp = psamples[j * m_numNodes + nnode] - pBetas[ic];
			ftemp = ftemp * ftemp;
			fsum += m_qc[j*c_setSize + ic] * ftemp;
		}
	}
	*pSigma = 0;
	if(fsum > 0) 
		*pSigma = fsum / nsamples;
	//else Rprintf("fsum == 0\n");
	if(*pSigma > 0)
		fBetaSigmaLoglik = -0.5*(log(PI2**pSigma) + 1)*nsamples;
	else 
		fBetaSigmaLoglik = FLT_MAX;
	//Rprintf("node = %d, fsum = %f, fBetaSigmaLoglik = %f, sigma = %f\n", nnode, fsum, fBetaSigmaLoglik, *pSigma);
	return fBetaSigmaLoglik;
}

int CMixNet::findNodeMarginalProb(int nnode, double *psamples, int nsamples) {

	int res, j, ic, numcats;
	double sigma2, fdiff, fsum, *paux;

	if(nnode < 0 || nnode >= m_numNodes || !m_sigmas || !m_betas || 
		!psamples || nsamples < 1)
		return ERR_CATNET_PARAM;

	if(m_pc) 
		CATNET_FREE(m_pc);
	m_pc = 0;
	if(m_qc) 
		CATNET_FREE(m_qc);
	m_qc = 0;

	m_curNode = nnode;

	sigma2 = FLT_MAX;
	if(m_sigmas[nnode] > 0)
		sigma2 =  1 / (2*m_sigmas[nnode]);
	numcats = m_numCategories[nnode];

	m_pc = (double*)CATNET_MALLOC(numcats*sizeof(double));

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
	paux = (double*)CATNET_MALLOC(numcats*sizeof(double));
	for (j = 0; j < nsamples; j++) {
		fsum = 0;
		for(ic = 0; ic < numcats; ic++) {
			fdiff = psamples[j * m_numNodes + nnode] - m_betas[nnode][ic];
			paux[ic] = m_pc[ic]*exp(-sigma2*fdiff*fdiff); 
			fsum += paux[ic];
		}
		if(fsum > 0) 
			fsum = 1/fsum;
		else
			fsum = FLT_MAX;
		for(ic = 0; ic < numcats; ic++) {
			m_qc[j*numcats + ic] = fsum*paux[ic];
		}
	}
	CATNET_FREE(paux);
	paux = 0;

	return ERR_CATNET_OK;
}

// m_pcC returns marginal P(X_i,X_{Pa_i}) 
// m_qcC returns P(X_i,X_{Pa_i}|Y_i,Y_{Pa_i}) for each sample instance y
int CMixNet::findNodeJointProb(int nnode, int *parnodes, int numpars, double *psamples, int nsamples) {

	int res, i, j, ic, icc;
	int *pnodes, numnodes, *pcats, numcats;
	double *psigma2, fdiff, fdiffsum, fsum, *paux;

//Rprintf("CMixNet::parentSetProb %d, %d, %p, %d\n", nnode, numpars, psamples, nsamples);

	if(nnode < 0 || nnode >= m_numNodes || !m_sigmas || !m_betas || 
		!psamples || nsamples < 1)
		return ERR_CATNET_PARAM;
		
	if(!parnodes || !numpars)
		return ERR_CATNET_OK;
		
	// pnodes has all parent nodes plus nnode attached at the end
	pnodes = (int*)CATNET_MALLOC((numpars+1)*sizeof(int));
	memcpy(pnodes, parnodes, numpars*sizeof(int));
	pnodes[numpars] = nnode;
	numnodes = numpars + 1;

//Rprintf("pnodes = ");
//for(i = 0; i < numnodes; i++)
//Rprintf("%d ", pnodes[i]);
//Rprintf("\n");

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
	for(i = 0; i < numnodes; i++) {
		m_parCatSetSize *= m_numCategories[pnodes[i]];
	}

	psigma2 = (double*)CATNET_MALLOC(numnodes*sizeof(double));
	memset(psigma2, 0, numnodes*sizeof(double));
	for(i = 0; i < numnodes; i++) {
		if(m_sigmas[pnodes[i]] > 0)
			psigma2[i] = 1 / (2*m_sigmas[pnodes[i]]);
		else
			psigma2[i] = FLT_MAX;
	}

	paux = (double*)CATNET_MALLOC(m_parCatSetSize*sizeof(double));
	pcats = (int*)CATNET_MALLOC(m_parCatSetSize*numnodes*sizeof(int));
	m_pcC = (double*)CATNET_MALLOC(m_parCatSetSize*sizeof(double));
	m_qcC = (double*)CATNET_MALLOC(m_parCatSetSize*nsamples*sizeof(double));

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
		for (i = numnodes - 2; i >= 0; i--) {
			m_pBlockSizes[i] = m_pBlockSizes[i + 1] * m_numCategories[pnodes[i + 1]];
			
		}

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
		fsum = 0;
		for(ic = 0; ic < m_parCatSetSize; ic++) {
			// a blunder: only one beta per nnode
			// would be nice to have m_betas[nnode][pcats[ic]]
			//fdiff = psamples[j * m_numNodes + nnode] - m_betas[nnode][pcats[numnodes-1]];
			// The way to GO
			fdiffsum = 0;
			for(i = 0; i < numnodes; i++) {
				fdiff = psamples[j * m_numNodes + pnodes[i]] - m_betas[pnodes[i]][pcats[ic * numnodes + i]];
				fdiffsum += psigma2[i]*fdiff*fdiff;
			}
			paux[ic] = m_pcC[ic]*exp(-fdiffsum);
			fsum += paux[ic];
		}
		if(fsum > 0) 
			fsum = 1/fsum;
		else
			fsum = FLT_MAX;
		for(ic = 0; ic < m_parCatSetSize; ic++)
			m_qcC[j*m_parCatSetSize + ic] = fsum*paux[ic];
	}
	CATNET_FREE(pnodes);
	CATNET_FREE(psigma2);
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

// m_pcC returns marginal P(X_{Pa_i}) 
// returns P(X_{Pa_i}|Y_{Pa_i}) for each sample instance y
int CMixNet::findParentsJointProb(int *parnodes, int numpars, double *psamples, int nsamples) {

	int res, i, j, ic;
	int *pcats;
	double *psigma2, fdiff, fdiffsum, fsum, *paux;

	if(!m_sigmas || !m_betas || !psamples || nsamples < 1)
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

	psigma2 = (double*)CATNET_MALLOC(numpars*sizeof(double));
	memset(psigma2, 0, numpars*sizeof(double));
	for(i = 0; i < numpars; i++) {
		if(m_sigmas[parnodes[i]] > 0)
			psigma2[i] = 1 / (2*m_sigmas[parnodes[i]]);
		else
			psigma2[i] = FLT_MAX;
	}

	paux = (double*)CATNET_MALLOC(m_parCatSetSize*sizeof(double));
	m_qcC = (double*)CATNET_MALLOC(m_parCatSetSize*nsamples*sizeof(double));
	if(!m_qcC || !paux)
		return ERR_CATNET_MEM;

	for (j = 0; j < nsamples; j++) {
		fsum = 0;
		for(ic = 0; ic < m_parCatSetSize; ic++) {
			fdiffsum = 0;
			for(i = 0; i < numpars; i++) {
				fdiff = psamples[j * m_numNodes + parnodes[i]] - m_betas[parnodes[i]][pcats[ic * numpars + i]];
				fdiffsum += psigma2[i]*fdiff*fdiff;
			}
			paux[ic] = m_pcC[ic]*exp(-fdiffsum);
			fsum += paux[ic];
		}
		if(fsum > 0) 
			fsum = 1/fsum;
		else
			fsum = FLT_MAX;
		for(ic = 0; ic < m_parCatSetSize; ic++)
			m_qcC[j*m_parCatSetSize + ic] = fsum*paux[ic];
	}

	CATNET_FREE(psigma2);
	CATNET_FREE(pcats);
	CATNET_FREE(paux);

	return ERR_CATNET_OK;
}

int* CMixNet::catnetSample(int nsamples) {

	int *porder;
	int i, j, k, nnode, *pnodepars, *pnodesample;
	double u, v, *pnodeprob;
	PROB_LIST<double>* pProbList;

	if(nsamples < 1 || m_numNodes < 1)
		return 0;
	porder = getOrder();
	if(!porder) {
		warning("porder is NULL");
		return 0;
	}

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
			if(!pnodeprob) {
				if(m_pCatnetSamples)
					CATNET_FREE(m_pCatnetSamples);
				m_pCatnetSamples = 0;
				m_nCatnetSamples = 0;
				CATNET_FREE(pnodesample);
				CATNET_FREE(porder);
				return 0;
			}
			u = (double)rand() / (double)RAND_MAX;
			v = 0;
			for(i = 0; i < m_numCategories[nnode]; i++) {
				v += pnodeprob[i];
				if(u <= v)
					break;
			}
			m_pCatnetSamples[j * m_numNodes + nnode] = (int)i;
		}
	}

	CATNET_FREE(pnodesample);
	CATNET_FREE(porder);

	return m_pCatnetSamples;
}

int CMixNet::sample(double *psamples, int nsamples) {

	int *porder;
	int i, j, k, nnode, *pnodepars, *pnodesample;
	double u, v, *pnodeprob;
	PROB_LIST<double>* pProbList; 

	if(!m_sigmas || !m_betas)
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
		u = sqrt(m_sigmas[nnode]);
		for (j = 0; j < nsamples; j++) {
			v = _gen_std_normal_var<double>();
			i = (int)psamples[j * m_numNodes + nnode];
			psamples[j * m_numNodes + nnode] = m_betas[nnode][i] + v*u;
		}
	}
	
	CATNET_FREE(porder);

	return ERR_CATNET_OK;
}

int CMixNet::predict(double *psamples, int nsamples) {

	int *porder;
	int j, k, ic, icC, nnode, cC_setSize;
	double fsum, faux, *pnodeprob;
	PROB_LIST<double>* pProbList;

	if(!m_sigmas || !m_betas)
		return ERR_CATNET_INIT;
	
	if(!psamples || nsamples < 1)
		return ERR_CATNET_PARAM;	

	porder = getOrder();
	if(!porder)
		return ERR_CATNET_PROC;

	for(k = 0; k < m_numNodes; k++) {
		/* work with the sample distribution of pCurNet if necessary */
		set_sample_cache();

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
					fsum += (m_betas[nnode][ic] * pnodeprob[ic]);
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
			fsum = 0;
			for(ic = 0; ic < m_numCategories[nnode]; ic++) {
				faux = 0;
				for(icC = 0; icC < cC_setSize; icC++) {
					faux += pnodeprob[icC*m_numCategories[nnode] + ic] * 
						m_qcC[j*cC_setSize + icC];
				}
				fsum += m_betas[nnode][ic] * faux;
			}
			psamples[j*m_numNodes + nnode] = fsum;
		}

		if(m_pCatnetSamples)
			CATNET_FREE(m_pCatnetSamples);
		m_nCatnetSamples = 0;
		m_pCatnetSamples = 0; 
	}

	CATNET_FREE(porder);

	return ERR_CATNET_OK;
}

double CMixNet::findLogNodeLikelihood(int nnode, double *psamples, int nsamples) {
	
	double faux, fsum, fsumsum, fLogLik, fconst1, fconst2;
	int j, cC_setSize, ic, icC;
	double *pnodeprob;
	PROB_LIST<double>* pProbList;
	
	if(!m_sigmas || !m_betas || nnode < 0 || nnode >= m_numNodes || m_sigmas[nnode] <= 0)
		return -FLT_MAX;
	
	if(!psamples || nsamples < 1)
		return -FLT_MAX;

	pProbList = (PROB_LIST<double>*)getNodeProb(nnode);
	pnodeprob = pProbList->pProbs;

	//Rprintf("prob %d: ", nnode);
	//for(ic = 0; ic < pProbList->nProbSize; ic++) {
		//Rprintf(" %f ", pnodeprob[ic]);
	//}
	//Rprintf("\n");
			
	cC_setSize = 0;
	if(m_numParents[nnode] > 0) {
		if(findParentsJointProb(m_parents[nnode], m_numParents[nnode], psamples, nsamples) != ERR_CATNET_OK)
			return -FLT_MAX;
		// now m_qcC contains P(x_j, j\in J_i | y_j^s, j\in J_i), j = 1,...,numsamples
		cC_setSize = pcC_size();
	}

	fconst1 = -0.5 * log(PI2*m_sigmas[nnode]);
	fconst2 = -0.5 / m_sigmas[nnode];

	fLogLik = 0;
	for(j = 0; j < nsamples; j++) {
		fsumsum = 0;
		for(ic = 0; ic < m_numCategories[nnode]; ic++) {
			faux = psamples[j * m_numNodes + nnode] - m_betas[nnode][ic];
			faux = fconst2 * faux * faux;
			if(m_numParents[nnode] > 0) {
				fsum = 0;
				for(icC = 0; icC < cC_setSize; icC++) {
					fsum += pnodeprob[icC*m_numCategories[nnode] + ic] * 
						m_qcC[j*cC_setSize + icC];
					//Rprintf("      %f, %f\n", pnodeprob[icC*m_numCategories[nnode] + ic], m_qcC[j*cC_setSize + icC]);
				}
			}
			else
				fsum = pnodeprob[ic];
			//Rprintf("faux = %f, fsum = %f\n", faux, fsum);
			fsumsum += exp(faux) * fsum;
		}
		fLogLik += log(fsumsum);
	}
	fLogLik += nsamples*fconst1;
	//Rprintf("nnode = %d:  fLogLik = %f\n", nnode, fLogLik);
	return fLogLik;
}


double CMixNet::findNodeLoglik(int nnode, double *psamples, int nsamples) {
	
	double faux, fLogLik, fconst1, fconst2;
	int res, j, cC_setSize, ic, icC, iC;
	double *pnodeprob;
	PROB_LIST<double>* pProbList;
	
	if(!m_sigmas || !m_betas || nnode < 0 || nnode >= m_numNodes || m_sigmas[nnode] <= 0)
		return -FLT_MAX;
	
	if(!psamples || nsamples < 1)
		return -FLT_MAX;

	pProbList = (PROB_LIST<double>*)getNodeProb(nnode);
	pnodeprob = pProbList->pProbs;
	
	fconst1 = -0.5 * log(PI2*m_sigmas[nnode]);
	fconst2 = -0.5 / m_sigmas[nnode];

	cC_setSize = 0;
	if(m_numParents[nnode] > 0) {
		if(findNodeJointProb(nnode, m_parents[nnode], m_numParents[nnode], psamples, nsamples) != ERR_CATNET_OK)
			return -FLT_MAX;
		// now m_qcC contains P(x_i, x_j, j\in J_i | y_i, y_j^s, j\in J_i), j = 1,...,numsamples
		cC_setSize = pcC_size();
	}
	else {
		res = findNodeMarginalProb(nnode, psamples, nsamples);
		if(res == ERR_CATNET_MEM) {
			int nsize = (int)exp(log((double)(maxCategories()))*(1+maxParentSet()))*1000;
			catnetSample(nsize);
			res = findNodeMarginalProb(nnode, psamples, nsamples);
		}
		if(res != ERR_CATNET_OK)
			return -FLT_MAX;
	}

	fLogLik = 0;
	if(m_numParents[nnode] > 0) {
		for(j = 0; j < nsamples; j++) {
			iC = 0;
			for(icC = 0; icC < cC_setSize; icC += m_numCategories[nnode]) {
				for(ic = 0; ic < m_numCategories[nnode]; ic++) {
					faux = psamples[j * m_numNodes + nnode] - m_betas[nnode][ic];
					faux = fconst2 * faux * faux;
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
				faux = psamples[j * m_numNodes + nnode] - m_betas[nnode][ic];
				faux = fconst2 * faux * faux;
				if(pnodeprob[ic] > 0)
					fLogLik += m_qc[j*m_numCategories[nnode] + ic] * (faux + log(pnodeprob[ic]));
				else
					fLogLik += m_qc[j*m_numCategories[nnode] + ic] * faux;
			}
		}
	}
	
	fLogLik += nsamples*fconst1;
	//Rprintf("nnode = %d:  fLogLik = %f\n", nnode, fLogLik);
	return fLogLik;
}

int CMixNet::estimateNodeProb(int nnode, int *parnodes, int numpars, double *psamples, int nsamples) {

	int res, i, j, ic, iC;
	int *pcats, numcats;
	double *psigma2, fsigma2, fdiff, fdiffsum, fsum, *paux;
	double *pnodeprob;
	PROB_LIST<double>* pProbList;

	if(nnode < 0 || nnode >= m_numNodes || !m_sigmas || !m_betas || 
		!psamples || nsamples < 1)
		return ERR_CATNET_PARAM;

	pProbList = (PROB_LIST<double>*)getNodeProb(nnode);
	pnodeprob = pProbList->pProbs;

	numcats = m_numCategories[nnode];
	fsigma2 = 0.5 / m_sigmas[nnode];

	if(numpars <= 0 || !parnodes) {
		fsum = 0;
		for(ic = 0; ic < numcats; ic++) {
			fdiffsum = 0;
			for(j = 0; j < nsamples; j++) {
				fdiff = psamples[j * m_numNodes + nnode] - m_betas[nnode][ic];
				fdiff = fsigma2 * fdiff * fdiff;
				// P(Y_i|X_i=c)
				fdiffsum += exp(-fdiff);
			}
			pnodeprob[ic] = fdiffsum;
			fsum += fdiffsum;
		}
		if(fsum > 0)
			fsum = 1/fsum;
		for(ic = 0; ic < numcats; ic++)
			pnodeprob[ic] *= fsum;
		return ERR_CATNET_OK;
	}

	// at this point m_nBlockSizes == numnodes and m_pBlockSizes are filled
	m_parCatSetSize = 1;
	for(i = 0; i < numpars; i++) {
		m_parCatSetSize *= m_numCategories[parnodes[i]];
	}

	if(m_pcC) 
		CATNET_FREE(m_pcC);
	m_pcC = 0;
	if(m_qcC) 
		CATNET_FREE(m_qcC);
	m_qcC = 0;

	// m_pcC returns marginal P(X_{Pa_i}) 
	m_pcC = (double*)CATNET_MALLOC(m_parCatSetSize*sizeof(double));
	pcats = (int*)CATNET_MALLOC(m_parCatSetSize*numpars*sizeof(int));
	if(!m_pcC || !pcats)
		return ERR_CATNET_MEM;
	memset(m_pcC, 0, m_parCatSetSize*sizeof(double));

	if(!m_pCatnetSamples || m_nCatnetSamples < 1) {
		// find the exact joint probability for small sets
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

	psigma2 = (double*)CATNET_MALLOC(numpars*sizeof(double));
	memset(psigma2, 0, numpars*sizeof(double));
	for(i = 0; i < numpars; i++) {
		if(m_sigmas[parnodes[i]] > 0)
			psigma2[i] = 0.5 / m_sigmas[parnodes[i]];
		else
			psigma2[i] = FLT_MAX;
	}

	paux = (double*)CATNET_MALLOC(m_parCatSetSize*sizeof(double));
	m_qcC = (double*)CATNET_MALLOC(m_parCatSetSize*nsamples*sizeof(double));

	for (j = 0; j < nsamples; j++) {
		fsum = 0;
		for(ic = 0; ic < m_parCatSetSize; ic++) {
			fdiffsum = 0;
			for(i = 0; i < numpars; i++) {
				fdiff = psamples[j * m_numNodes + parnodes[i]] - m_betas[parnodes[i]][pcats[ic * numpars + i]];
				fdiffsum += psigma2[i]*fdiff*fdiff;
			}
			paux[ic] = m_pcC[ic]*exp(-fdiffsum);
			fsum += paux[ic];
		}
		if(fsum > 0) 
			fsum = 1/fsum;
		// set P(X_{Pa_i}|Y_{Pa_i}) in m_qcC
		for(ic = 0; ic < m_parCatSetSize; ic++)
			m_qcC[j*m_parCatSetSize + ic] = fsum*paux[ic];
	}

	for(iC = 0; iC < m_parCatSetSize; iC++) {
		fsum = 0;
		for(ic = 0; ic < numcats; ic++) {
			fdiffsum = 0;
			for(j = 0; j < nsamples; j++) {
				fdiff = psamples[j * m_numNodes + nnode] - m_betas[nnode][ic];
				fdiff = fsigma2 * fdiff * fdiff;
				// P(Y_i|X_i=c)P(Pa_i=C|Y_{Pa_i})
				fdiffsum += exp(-fdiff)*m_qcC[j*m_parCatSetSize + iC];
			}
			pnodeprob[iC*numcats + ic] = fdiffsum;
			fsum += fdiffsum;
		}
		if(fsum > 0) {
			fsum = 1/fsum;
			for(ic = 0; ic < numcats; ic++)
				pnodeprob[iC*numcats + ic] *= fsum;
		}
		else {
			// set uniform prob
			fsum = 1/numcats;
			for(ic = 0; ic < numcats; ic++)
				pnodeprob[iC*numcats + ic] = fsum;
			fsum = 0;
			for(ic = 0; ic < numcats-1; ic++)
				fsum += pnodeprob[iC*numcats + ic];
			pnodeprob[iC*numcats + numcats-1] = 1 - fsum;
		}
	}

	CATNET_FREE(psigma2);
	CATNET_FREE(pcats);
	CATNET_FREE(paux);

	// we don't need them any more
	if(m_pcC) 
		CATNET_FREE(m_pcC);
	m_pcC = 0;
	if(m_qcC) 
		CATNET_FREE(m_qcC);
	m_qcC = 0;

	// update the betas and sigma 
	res = findNodeMarginalProb(nnode, psamples, nsamples);
	if(res != ERR_CATNET_OK)
		return res;
	estimateParameters(nnode, psamples, nsamples, m_betas[nnode], &m_sigmas[nnode]);

	if(m_pc) 
		CATNET_FREE(m_pc);
	m_pc = 0;
	if(m_qc) 
		CATNET_FREE(m_qc);
	m_qc = 0;

	return ERR_CATNET_OK;
}

int CMixNet::setProbability(double *psamples, int nsamples) {

	int *porder;
	int k, nnode;

	if(!m_sigmas || !m_betas)
		return ERR_CATNET_INIT;
	
	if(!psamples || nsamples < 1)
		return ERR_CATNET_PARAM;	

	porder = getOrder();
	if(!porder)
		return ERR_CATNET_PROC;

	for(k = 0; k < m_numNodes; k++) {
		
		/* work with the sample distribution of pCurNet if necessary */
		set_sample_cache();
/// check this out !
		nnode = porder[k];

		if(estimateNodeProb(nnode, m_parents[nnode], m_numParents[nnode], psamples, nsamples) != ERR_CATNET_OK)
			break;

		if(m_pCatnetSamples)
			CATNET_FREE(m_pCatnetSamples);
		m_nCatnetSamples = 0;
		m_pCatnetSamples = 0; 
	}

	CATNET_FREE(porder);

	return ERR_CATNET_OK;
}

void CMixNet::set_sample_cache(int becho) {
	if(	maxParentSet() > 3 
		|| maxCategories()*maxParentSet() > 9 
		|| m_numNodes*maxCategories()*maxParentSet() > 256
	) {
		int nsize = (int)exp(log((double)(maxCategories()))*(1+maxParentSet()))*1000;
		if(becho)
			Rprintf("simulate sample cache %d\n", nsize);
		catnetSample(nsize);
	}
}

