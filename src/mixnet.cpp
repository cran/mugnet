/*
 * mixnet.cpp
 *
 *  Created on: Nov 16, 2009
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
	CATNET<char, MAX_NODE_NAME, double>::_release();
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
	if(m_pCatnetSamples)
		CATNET_FREE(m_pCatnetSamples);
	m_pCatnetSamples = 0; 
	m_nCatnetSamples = 0;
}

void CMixNet::_reset() {
	CATNET<char, MAX_NODE_NAME, double>::_reset();
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

	CATNET<char, MAX_NODE_NAME, double> ::init(
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

	sigma2 = 0;
	if(m_sigmas[nnode] > 0)
		sigma2 =  - 1 / (2*m_sigmas[nnode]);
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
			paux[ic] = m_pc[ic]*exp(sigma2*fdiff*fdiff); 
			fsum += paux[ic];
		}
		if(fsum > 0) 
			fsum = 1/fsum;
		for(ic = 0; ic < numcats; ic++) {
			m_qc[j*numcats + ic] = fsum*paux[ic];
		}
	}
	CATNET_FREE(paux);
	paux = 0;

	return ERR_CATNET_OK;
}


int CMixNet::findNodeJointProb(int nnode, int *parnodes, int numpars, double *psamples, int nsamples) {

	int res, i, j, ic, icc, k;
	int *pnodes, numnodes, *pcats, numcats;
	double *psigma2, fdiff, fdiffsum, fsum, *paux;

//printf("CMixNet::parentSetProb %d, %d, %p, %d\n", nnode, numpars, psamples, nsamples);

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

//printf("pnodes = ");
//for(i = 0; i < numnodes; i++)
//printf("%d ", pnodes[i]);
//printf("\n");

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
	}

	paux = (double*)CATNET_MALLOC(m_parCatSetSize*sizeof(double));
	pcats = (int*)CATNET_MALLOC(numnodes*sizeof(int));
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
				pcats[i] = (int)(fsum/m_pBlockSizes[i]);
				fsum -= pcats[i] * m_pBlockSizes[i];
			}
			m_pcC[ic] = getCatProb(pcats, numnodes);
		}
	}
	else {
		m_pBlockSizes[numnodes - 1] = 1;
		for (i = numnodes - 2; i >= 0; i--) {
			m_pBlockSizes[i] = m_pBlockSizes[i + 1] * m_numCategories[pnodes[i + 1]];
			
		}

		memset(m_pcC, 0, m_nCatnetSamples*sizeof(double));

		for(j = 0; j < m_nCatnetSamples; j++) {
			ic = 0;
			for (i = 0; i < numnodes; i++)
				ic += (m_pBlockSizes[i] * m_pCatnetSamples[j * m_numNodes + pnodes[i]]);
			m_pcC[ic] += 1;
		}

		fsum = 1/(double)m_nCatnetSamples;
		for(ic = 0; ic < m_parCatSetSize; ic++)
			m_pcC[ic] *= fsum;
	}

	for (j = 0; j < nsamples; j++) {
		fsum = 0;
		for(ic = 0; ic < m_parCatSetSize; ic++) {
			// determine the category indices
			k = ic;
			for(i = 0; i < numnodes; i++) {
				pcats[i] = (int)(k/m_pBlockSizes[i]);
				k -= pcats[i] * m_pBlockSizes[i];
			}

			// a blunder: only one beta per nnode
			// would be nice to have m_betas[nnode][pcats[ic]]
			//fdiff = psamples[j * m_numNodes + nnode] - m_betas[nnode][pcats[numnodes-1]];
			// The way to GO
			fdiffsum = 0;
			for(i = 0; i < numnodes; i++) {
				fdiff = psamples[j * m_numNodes + pnodes[i]] - m_betas[pnodes[i]][pcats[i]];
				fdiffsum += psigma2[i]*fdiff*fdiff;
			}
			paux[ic] = m_pcC[ic]*exp(-fdiffsum);
			fsum += paux[ic];
		}
		if(fsum > 0) fsum = 1/fsum;
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

int CMixNet::findParentsJointProb(int *parnodes, int numpars, double *psamples, int nsamples) {

	int res, i, j, ic, k;
	int *pcats;
	double *psigma2, fdiff, fdiffsum, fsum, *paux;

	if(!m_sigmas || !m_betas || !psamples || nsamples < 1)
		return ERR_CATNET_PARAM;
	if(!parnodes || !numpars)
		return ERR_CATNET_PARAM;

	res = marginalProb(parnodes, numpars);
	if(res != ERR_CATNET_OK) {
		return res;
	}

	if(m_pcC) CATNET_FREE(m_pcC);
	m_pcC = 0;
	if(m_qcC) CATNET_FREE(m_qcC);
	m_qcC = 0;

	// at this point m_nBlockSizes == numnodes and m_pBlockSizes are filled
	m_parCatSetSize = 1;
	for(i = 0; i < numpars; i++)
		m_parCatSetSize *= m_numCategories[parnodes[i]];
	//m_parCatSet = (int*)CATNET_MALLOC(m_parCatSetSize*sizeof(int));

	if(m_margProbSize != m_parCatSetSize) {
		// fatal  error
		return ERR_CATNET_PROC;
	}

	psigma2 = (double*)CATNET_MALLOC(numpars*sizeof(double));
	memset(psigma2, 0, numpars*sizeof(double));
	for(i = 0; i < numpars; i++) {
		if(m_sigmas[parnodes[i]] > 0)
			psigma2[i] = 1 / (2*m_sigmas[parnodes[i]]);
	}

	paux = (double*)CATNET_MALLOC(m_parCatSetSize*sizeof(double));
	pcats = (int*)CATNET_MALLOC(numpars*sizeof(int));

	m_pcC = (double*)CATNET_MALLOC(m_parCatSetSize*sizeof(double));
	for(ic = 0; ic < m_parCatSetSize; ic++) {
		// determine the category indices
		fsum = ic;
		for(i = 0; i < numpars; i++) {
			pcats[i] = (int)(fsum/m_pBlockSizes[i]);
			fsum -= pcats[i] * m_pBlockSizes[i];
		}
		m_pcC[ic] = getCatProb(pcats, numpars);
	}

	m_qcC = (double*)CATNET_MALLOC(m_parCatSetSize*nsamples*sizeof(double));
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
				fdiff = psamples[j * m_numNodes + parnodes[i]] - m_betas[parnodes[i]][pcats[i]];
				fdiffsum += psigma2[i]*fdiff*fdiff;
			}
			paux[ic] = m_pcC[ic]*exp(-fdiffsum);
			fsum += paux[ic];
		}
		if(fsum > 0) fsum = 1/fsum;
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
	int j, k, ic, icC, nnode, imax, cC_setSize;
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
if(1) {
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
			psamples[j*m_numNodes + nnode] = m_betas[nnode][imax];
}
		}
	}

	CATNET_FREE(porder);

	return ERR_CATNET_OK;
}

double CMixNet::findNodeLoglik(int nnode, double *psamples, int nsamples) {
	
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

	cC_setSize = 0;
	if(m_numParents[nnode] > 0) {
		if(findParentsJointProb(m_parents[nnode], m_numParents[nnode], psamples, nsamples) != ERR_CATNET_OK)
			return -FLT_MAX;
		// now m_qcC contains P(x_j, j\in J_i | y_j^s, j\in J_i), j = 1,...,numsamples
		cC_setSize = pcC_size();
	}

	fconst1 = log(1 / sqrt(PI2*m_sigmas[nnode]));
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
				}
			}
			else
				fsum = pnodeprob[ic];
			fsumsum += exp(faux) * fsum;
		}
		fLogLik += log(fsumsum);
	}
	fLogLik += nsamples*fconst1;
	return fLogLik;
}

