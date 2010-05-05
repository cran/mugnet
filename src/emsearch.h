/*
 * emsearch.h
 *
 *  Created on: Nov 16, 2009
 *      Author: Nikolay Balov
 */

#include <float.h>
#include <stdio.h>

#include "utils.h"
#include "inetparams.h"

#ifndef EMSEARCH_H
#define EMSEARCH_H

template<class t_prob>
class EMSEARCH {
protected:
	int m_nNets;
	I_NETPARAMS<t_prob> **m_pCurNets, **m_pNextNets;
	t_prob *m_emLoglik;

	/* returns increasing subsets of 'parset' of size 'parsize' */
	void combinationSets(int **&plist, int &nlist, int *curset, int *parset, int nparset, 
				int parid, int parsize) {
		int i, ancestor;
 
		if(parid < 0 || parid >= parsize)
			return;
		ancestor = -1;
		if(parid > 0)
			ancestor = curset[parid-1];

		if(parid == parsize - 1) {
			for(i = 0; i < nparset; i++) {
				if(parset[i] <= ancestor)
					continue;
				int **pnewlist = (int**)CATNET_MALLOC((nlist+1)*sizeof(int*));
				if(nlist > 0)
					memcpy(pnewlist, plist, nlist*sizeof(int*));
				pnewlist[nlist] = (int*)CATNET_MALLOC(parsize*sizeof(int));
				if(curset) {
					memcpy(pnewlist[nlist], curset, (parsize-1)*sizeof(int));
				}
				pnewlist[nlist][parsize-1] = parset[i];

				CATNET_FREE(plist);
				plist = pnewlist;
				nlist++;
			}
			if(curset) {
				CATNET_FREE(curset);
				curset = 0;
			}
			return;
		}
		for(i = 0; i < nparset; i++) {
			if(parset[i] <= ancestor)
				continue;
			int *pnewset = (int*)CATNET_MALLOC((parid+1)*sizeof(int));
			if(curset && parid > 0)
				memcpy(pnewset, curset, parid*sizeof(int));
			pnewset[parid] = parset[i];
			combinationSets(plist, nlist, pnewset, parset, nparset, parid+1, parsize);
		}
		if(curset) {
			CATNET_FREE(curset);
			curset = 0;
		}
	}

public:
	EMSEARCH() {
		m_nNets = 0;
		m_pCurNets = 0;
		m_pNextNets = 0;
		m_emLoglik = 0;
	}

	~EMSEARCH() {
		int i;
		if(m_pCurNets) {
			for(i = 0; i < m_nNets; i++)
				if(m_pCurNets[i]) {
					m_pCurNets[i]->releaseRef();
					m_pCurNets[i] = 0;
				}
			CATNET_FREE(m_pCurNets);
		}
		m_pCurNets = 0;
		if(m_pNextNets)
			CATNET_FREE(m_pNextNets);
		m_pNextNets = 0;
		m_nNets = 0;
		if(m_emLoglik)
			CATNET_FREE(m_emLoglik);
		m_emLoglik = 0;
	}

	/* psamples and perturbations are sample=columns and node=rows. */
	/* Each parentsPool[i] is numnodes long ! */
int estimate(
	I_NETPARAMS<t_prob> *pBaseNet,
	int numSamples, t_prob *pSamples, int *perturbations,
	int maxParentSet, int maxComplexity,
	int **parentsPool, int **fixedParentsPool, 
	int emIterations, double stopDelta, int netSelection = 0, int becho = FALSE) {

	int i, j, k, d, ncomb, ncombMaxLogLik, nnode, numNodes;
	int maxCategories, complx, emiter, curSelComplexity;
	const int *pNodeNumCats;
	int *parset, parsetsize, *fixparset, fixparsetsize, *parcats;
	int *paux, **pcomblist, ncomblist, maxpars, ballow, bfixallow;
		
	int nCurNet, nLoopNet;
	I_NETPARAMS<t_prob> *pCurNet, *pNewNet, **pCurCatnetList;

	t_prob ftemp, fsum, fsumsum, curSelLoglik, logdiff;
	int c_setSize, cC_setSize, ic, icC;
	const t_prob *pc, *qc, *pcC, *qcC;

	t_prob **pBetasNext, *pSigmasNext;
	t_prob *pNodeBetas, fNodeSigma;

	t_prob *psamples, *psubsamples;
	int numsamples;

	t_prob fLogLik, fMaxLogLik, fBetaSigmaLoglik, *pprobs;
	PROB_LIST<t_prob> *pProbMaxNode;

	if(!pBaseNet || pBaseNet->numNodes() <= 0 || pBaseNet->complexity() < pBaseNet->numNodes())
		return ERR_CATNET_PARAM;

	if(numSamples < 1 || !pSamples)
		return ERR_CATNET_PARAM;

	if(emIterations < 2)
		emIterations = 2;

	numNodes = pBaseNet->numNodes();
	//maxParentSet < pBaseNet->maxParentSet();
	maxCategories = pBaseNet->maxCategories();
	pNodeNumCats = pBaseNet->numCategories();

	if(maxComplexity < pBaseNet->complexity())
		maxComplexity = pBaseNet->complexity();

	if(m_pCurNets) {
		for(i = 0; i < m_nNets; i++)
			if(m_pCurNets[i]) {
				m_pCurNets[i]->releaseRef();
				m_pCurNets[i] = 0;
			}
		CATNET_FREE(m_pCurNets);
	}
	m_pCurNets = 0;

	m_nNets = maxComplexity + 1;
	m_pCurNets = (I_NETPARAMS<t_prob>**)CATNET_MALLOC(m_nNets*sizeof(I_NETPARAMS<t_prob>*));
	memset(m_pCurNets, 0, m_nNets*sizeof(I_NETPARAMS<t_prob>*));

	// allocate a new net identical to pBaseNet
	complx = pBaseNet->complexity();
	m_pCurNets[complx] = pBaseNet->clone();
	// set the lowest possible likelihood
	if(m_pCurNets[complx])
		m_pCurNets[complx]->setLoglik(-FLT_MAX);

	if(m_pNextNets)
		CATNET_FREE(m_pNextNets);
	m_pNextNets = (I_NETPARAMS<t_prob>**)CATNET_MALLOC(m_nNets*sizeof(I_NETPARAMS<t_prob>*));
	memset(m_pNextNets, 0, m_nNets*sizeof(I_NETPARAMS<t_prob>*));

	pCurCatnetList = (I_NETPARAMS<t_prob>**)CATNET_MALLOC(m_nNets*sizeof(I_NETPARAMS<t_prob>*));

	parset = (int*)CATNET_MALLOC(numNodes*sizeof(int));
	fixparset = (int*)CATNET_MALLOC(numNodes*sizeof(int));
	parcats = (int*)CATNET_MALLOC(maxParentSet*sizeof(int));

	if(m_emLoglik)
		CATNET_FREE(m_emLoglik);
	m_emLoglik = (t_prob*)CATNET_MALLOC(m_nNets*emIterations*sizeof(t_prob)); 

	pBetasNext = (t_prob**)CATNET_MALLOC(numNodes*sizeof(t_prob*));
	pNodeBetas = (t_prob*)CATNET_MALLOC(pBaseNet->maxCategories()*sizeof(t_prob));
	for(i = 0; i < numNodes; i++) {
		pBetasNext[i] = (t_prob*)CATNET_MALLOC(pNodeNumCats[i]*sizeof(t_prob));
	}
	pSigmasNext = (t_prob*)CATNET_MALLOC(numNodes*sizeof(t_prob));

	psubsamples = 0;
	if(perturbations) {
		psubsamples = (t_prob*)CATNET_MALLOC(numNodes*numSamples*sizeof(t_prob));
	}

	/* Main EM Loop */
	emiter = 0;
	nLoopNet = 0;
	while(emiter < emIterations) {

	switch(netSelection) {
	case 1: // AIC
		pCurNet = 0;
		fLogLik = -FLT_MAX;
		for(nCurNet = 0; nCurNet < m_nNets; nCurNet++) {
			if(!m_pCurNets[nCurNet])
				continue;
			// net complexity = nCurNet + numNodes*(1 + maxCategories)
			complx = m_pCurNets[nCurNet]->getLoglik() - m_pCurNets[nCurNet]->complexity();
			if(complx > fLogLik) {
				fLogLik = complx;
				pCurNet = m_pCurNets[nCurNet];
			}
		}
		break;
	case 2:  // BIC
		pCurNet = 0;
		fLogLik = -FLT_MAX;
		ftemp = 0.5*log(numSamples);
		for(nCurNet = 0; nCurNet < m_nNets; nCurNet++) {
			if(!m_pCurNets[nCurNet])
				continue;
			complx = m_pCurNets[nCurNet]->getLoglik() - ftemp * m_pCurNets[nCurNet]->complexity();
			if(complx > fLogLik) {
				fLogLik = complx;
				pCurNet = m_pCurNets[nCurNet];
			}
		}
		break;
	default:
		// optimize w.r.t. max complexity network
		pCurNet = 0;
		for(nCurNet = 0; nCurNet < m_nNets; nCurNet++)
			if(m_pCurNets[nCurNet]) {
				pCurNet = m_pCurNets[nCurNet];
			}
	}
	if(!pCurNet)
		break;
	memset(m_pNextNets, 0, m_nNets*sizeof(I_NETPARAMS<t_prob>*));	

	curSelComplexity = pCurNet->complexity();
	curSelLoglik = pCurNet->getLoglik();

	// work with the sample distribution of pCurNet if necessary
	if(pCurNet->maxParentSet() > 2 || pCurNet->maxCategories()*pCurNet->maxParentSet() > 5 || 
		numNodes*pCurNet->maxCategories()*pCurNet->maxParentSet() > 80) {
		k = (int)exp(log(pCurNet->maxCategories())*(1+pCurNet->maxParentSet())) * 50;
		if(becho)
			printf("simulate probability with sample of size %d\n", k);
		pCurNet->catnetSample(k);
	}

	//printf("pCurNet: complx = %d, logLik = %f\n", curSelComplexity, curSelLoglik);

	/* update betas and sigmas */

	/* copy the current beta set !!! */
	for(i = 0; i < numNodes; i++) 
		memcpy(pBetasNext[i], pCurNet->betas()[i], pNodeNumCats[i]*sizeof(t_prob));
	memcpy(pSigmasNext, pCurNet->sigmas(), numNodes*sizeof(t_prob));

	/* and create the no-parents-network */
	pNewNet = pCurNet->clone();
	if(!pNewNet)
		break;

	for(nnode = 0; nnode < numNodes; nnode++) {

		pNewNet->setParents(nnode, 0, 0);

		if(perturbations) {
			numsamples = 0;
			for(j = 0; j < numSamples; j++) {
				if(!perturbations[j * numNodes + nnode]) {
					memcpy(psubsamples + numsamples*numNodes, pSamples + j*numNodes, 
						numNodes*sizeof(t_prob));
					numsamples++;
				}
			}
			psamples = psubsamples;
			//printf("perturbations: nnode = %d, numsamples = %d\n", nnode, numsamples);
		}
		else {
			numsamples = numSamples;
			psamples = pSamples;
		}

		if(numsamples <= 0) {
			pCurNet->setNodeLoglik(nnode, 0);
			pNewNet->setNodeLoglik(nnode, 0);
			//printf("numsamples = 0,   pNewNet: nnode = %d, fNodeLogLik = %f\n", nnode, pNewNet->getNodeLoglik(nnode));
			continue;
		}

		pCurNet->findNodeMarginalProb(nnode, psamples, numsamples);
		c_setSize = pCurNet->pc_size();
		pc = pCurNet->pc();
		qc = pCurNet->qc();

		for(ic = 0; ic < pNodeNumCats[nnode]; ic++) {
			ftemp = 0;
			fsum = 0;
			for(j = 0; j < numsamples; j++) {
				ftemp += qc[j*c_setSize + ic] * psamples[j * numNodes + nnode]; 
				fsum += qc[j*c_setSize + ic];
			}
			if(fsum > 0) {
				/* update only if possible, don't put 0, remember the last good value instead */
				fsum = 1/fsum;
				pBetasNext[nnode][ic] = ftemp * fsum;
			}
			//printf("fsum = %f, pBetasNext[%d][%d] = %f\n", fsum, nnode, ic, pBetasNext[nnode][ic]);
		}

		fsum = 0;
		for(ic = 0; ic < pNodeNumCats[nnode]; ic++) {
			for(j = 0; j < numsamples; j++) {
				ftemp = psamples[j * numNodes + nnode] - pBetasNext[nnode][ic];
				ftemp = ftemp * ftemp;
				fsum += qc[j*c_setSize + ic] * ftemp;
			}
		}
		if(fsum > 0) 
			pSigmasNext[nnode] = fsum / numsamples;
		if(pSigmasNext[nnode] > 0)
			fBetaSigmaLoglik = -0.5*(log(PI2*pSigmasNext[nnode]) + 1)*numsamples;
		else
			fBetaSigmaLoglik = -FLT_MAX;

//printf("pSigmasNext[%d] = %f\n", nnode, pSigmasNext[nnode]);
//printf("fsum = %f, numsamples = %d, fBetaSigmaLoglik = %f\n", fsum, numsamples, fBetaSigmaLoglik);
//printf("node = %d, fBetaSigmaLoglik = %f\n", nnode+1, fBetaSigmaLoglik);

		/* refresh probabilities  */
		pprobs = (t_prob*)CATNET_MALLOC(c_setSize*sizeof(t_prob));

		fsumsum = 0;
		ftemp = 0;
		for(ic = 0; ic < c_setSize/*pNodeNumCats[nnode]*/; ic++) {
			fsum = 0;
			for(j = 0; j < numsamples; j++)
				fsum += qc[j*c_setSize + ic];
			ftemp += fsum * log(fsum);
			fsumsum += fsum;
			pprobs[ic] = fsum;
		}
		if(fsumsum > 0) {
			fLogLik = (ftemp - fsumsum*log(fsumsum));
			fsumsum = 1/ fsumsum;
		}
		else
			fLogLik = -FLT_MAX;
		for(ic = 0; ic < c_setSize; ic++)
			pprobs[ic] = pprobs[ic] * fsumsum;
		pProbMaxNode = new PROB_LIST<t_prob>(pNodeNumCats[nnode], 
			maxCategories, 0, 0, (double*)pprobs, c_setSize);

		CATNET_FREE(pprobs);
		if(!pProbMaxNode)
			CATNET_MEM_ERR();

		/* should return 0 */
		pNewNet->setNodeProb(nnode, pProbMaxNode);
		if(pProbMaxNode) {
			delete pProbMaxNode;
			pProbMaxNode = 0;
		}

		pCurNet->release_pqcC();

		// you'll need these
		pCurNet->setNodeLoglik(nnode, fBetaSigmaLoglik);
		pNewNet->setNodeLoglik(nnode, fBetaSigmaLoglik + fLogLik);

		//printf("pNewNet: nnode = %d, fNodeLogLik = %f\n", nnode, pNewNet->getNodeLoglik(nnode));
	} /* nnode */

	// refresh curNet's betas and sigmas
	pCurNet->setBetas(pBetasNext, numNodes);
	pCurNet->setSigmas(pSigmasNext, numNodes);

	pNewNet->setBetas(pBetasNext, numNodes);
	pNewNet->setSigmas(pSigmasNext, numNodes);

	// setup the seed-net
	complx = pNewNet->complexity();
	if(complx < m_nNets)
		m_pNextNets[complx] = pNewNet;
	else {
		pNewNet->releaseRef();
		break;
	}

	//printf("pNewNet %p: complx = %d, fNodeLogLik = %f\n", pNewNet, complx, pNewNet->loglik());

	/* update G_i */
	for(nnode = 0; nnode < numNodes; nnode++) {

		if(perturbations) {
			numsamples = 0;
			for(j = 0; j < numSamples; j++) {
				if(!perturbations[j * numNodes + nnode]) {
					memcpy(psubsamples + numsamples*numNodes, pSamples + j*numNodes, 
						numNodes*sizeof(t_prob));
					numsamples++;
				}
			}
			psamples = psubsamples;
		}
		else {
			numsamples = numSamples;
			psamples = pSamples;
		}

		fixparsetsize = 0;
		parsetsize = 0;
		for(j = 0; j < nnode; j++) {
			ballow = 1;
			if(parentsPool && parentsPool[nnode]) {
				ballow = 0;
				for(k = 0; k < numNodes; k++) {
					if(j == parentsPool[nnode][k])
						ballow = 1;
				}
			}
			if(parentsPool && !parentsPool[nnode])
				ballow = 0;
			bfixallow = 0;
			if(fixedParentsPool && fixedParentsPool[nnode]) {
				for(k = 0; k < numNodes; k++) {
					if(j == fixedParentsPool[nnode][k])
						bfixallow = 1;
				}
			}
			if(!ballow)
				continue;
			if(bfixallow) {
			  fixparset[fixparsetsize] = j;
			  fixparsetsize++;
			}
			else {
			  parset[parsetsize] = j;
			  parsetsize++;
			}
		}

		maxpars = maxParentSet;
		if(maxpars > parsetsize + fixparsetsize)
			maxpars = parsetsize + fixparsetsize;

		if(numsamples <= 0) {
			// prevent parent search
			maxpars = 0;
		}

		memset(pCurCatnetList, 0, m_nNets*sizeof(I_NETPARAMS<t_prob>*));

		for(d = fixparsetsize + 1; d <= maxpars; d++) {

			pcomblist = 0;
			ncomblist = 0;
			combinationSets(pcomblist, ncomblist, 0, parset, parsetsize, 
				0, d - fixparsetsize);
			//printf("nnode = %d, d = %d, ncomblist = %d\n", nnode,  d, ncomblist);
			if(fixparsetsize > 0) {
		        	if(!pcomblist || ncomblist < 1) {
		        	    	pcomblist = (int**)CATNET_MALLOC(1*sizeof(int*));
		            		pcomblist[0] = 0;
		            		ncomblist = 1;
		          	}
		        	for(k = 0; k < ncomblist; k++) {
		            		paux = (int*)CATNET_MALLOC(d*sizeof(int));
					for(j = 0; j < fixparsetsize; j++)
		            			paux[j] = fixparset[j];
			        	if(pcomblist[k] && d > fixparsetsize) {
		            			memcpy(paux + fixparsetsize, pcomblist[k], 
							(d-fixparsetsize)*sizeof(int));
					}
		            		if(pcomblist[k])
		            			CATNET_FREE(pcomblist[k]); 
		           		pcomblist[k] = paux;
				}
			}

			pProbMaxNode = 0;
			fMaxLogLik = -FLT_MAX;
			ncombMaxLogLik = -1;
			fNodeSigma = 0;

			for(ncomb = 0; ncomb < ncomblist; ncomb++) {
				if(pCurNet->findNodeJointProb(nnode, pcomblist[ncomb], d, 
					psamples, numsamples) != ERR_CATNET_OK)
					continue; 
				// find betas and sigmas
				c_setSize = pCurNet->pc_size();
				pc = pCurNet->pc();
				qc = pCurNet->qc();
				for(ic = 0; ic < pNodeNumCats[nnode]; ic++) {
					ftemp = 0;
					fsum = 0;
					for(j = 0; j < numsamples; j++) {
						ftemp += qc[j*c_setSize + ic] * psamples[j * numNodes + nnode]; 
						fsum += qc[j*c_setSize + ic];
					}
					if(fsum > 0) {
						fsum = 1/fsum;
						pBetasNext[nnode][ic] = ftemp * fsum;
					}
				}
				fsum = 0;
				for(ic = 0; ic < pNodeNumCats[nnode]; ic++) {
					for(j = 0; j < numsamples; j++) {
						ftemp = psamples[j * numNodes + nnode] - pBetasNext[nnode][ic];
						ftemp = ftemp * ftemp;
						fsum += qc[j*c_setSize + ic] * ftemp;
					}
				}
				if(fsum > 0) 
					pSigmasNext[nnode] = fsum / numsamples;
				if(pSigmasNext[nnode] > 0)
					fBetaSigmaLoglik = -0.5*(log(PI2*pSigmasNext[nnode]) + 1)*numsamples;
				else
					fBetaSigmaLoglik = -FLT_MAX;

				cC_setSize = pCurNet->pcC_size();
				pcC = pCurNet->pcC();
				qcC = pCurNet->qcC();
				pprobs = (t_prob*)CATNET_MALLOC(cC_setSize*sizeof(t_prob));

				// nnode cats are in the last slots of pcC and qcC
				fLogLik = 0;
				for(icC = 0; icC < cC_setSize; icC += pNodeNumCats[nnode]) {
					fsumsum = 0;
					ftemp = 0;
					for(ic = 0; ic < pNodeNumCats[nnode]; ic++) {
						fsum = 0;
						for(j = 0; j < numsamples; j++)
							fsum += qcC[j*cC_setSize + icC + ic];
						ftemp += fsum * log(fsum);
						fsumsum += fsum;
						pprobs[icC + ic] = fsum;
					}
					fLogLik += (ftemp - fsumsum*log(fsumsum));
					if(fsumsum > 0) fsumsum = 1/fsumsum;
					for(ic = 0; ic < pNodeNumCats[nnode]; ic++)
						pprobs[icC + ic] = pprobs[icC + ic] * fsumsum;
				}

				// sum of two conditionals
				fLogLik += fBetaSigmaLoglik;

				if(fMaxLogLik < fLogLik) {
					memcpy(pNodeBetas, pBetasNext[nnode], pNodeNumCats[nnode]*sizeof(t_prob));
					fNodeSigma = pSigmasNext[nnode];

					if(pProbMaxNode)
						delete pProbMaxNode;
					for(i = 0; i < d; i++)
						parcats[i] = 
							pNodeNumCats[pcomblist[ncomb][i]];
					pProbMaxNode = 
						new PROB_LIST<t_prob>(pNodeNumCats[nnode],
						maxCategories, d, parcats,
						pprobs, cC_setSize);
					if(!pProbMaxNode)
						CATNET_MEM_ERR();
					fMaxLogLik = fLogLik;
					ncombMaxLogLik = ncomb;
				}
			
				CATNET_FREE(pprobs);
				pCurNet->release_pqcC();

			} /* for ncomb */				

			if(ncombMaxLogLik >= 0 && pProbMaxNode) {
				for(k = 0; k < m_nNets; k++) {
					if(!m_pNextNets[k]) 
						continue;
					pNewNet = m_pNextNets[k]->clone();
					//fLogLik = pNewNet->getNodeLoglik(nnode);

					pNewNet->setNodeBetas(nnode, pNodeBetas);
					pNewNet->setNodeSigma(nnode, fNodeSigma);

					pNewNet->setNodeLoglik(nnode, fMaxLogLik);
					pNewNet->setParents(nnode, pcomblist[ncombMaxLogLik], d);
					complx = pNewNet->complexity();
					if(complx > maxComplexity/*==m_nNets+1*/) {
						pNewNet->releaseRef();
						continue;
					}
					/* should return 0 */
					pNewNet->setNodeProb(nnode, pProbMaxNode);
					fLogLik = pNewNet->loglik();
					if(pCurCatnetList[complx]) {
						if(pCurCatnetList[complx]->loglik() < fLogLik) {
							pCurCatnetList[complx]->releaseRef();
							pCurCatnetList[complx] = pNewNet;
						}
						else {
							pNewNet->releaseRef();
						}
					}
					else {
						pCurCatnetList[complx] = pNewNet;
					}
				} /* for k */
			}

			if(pProbMaxNode) {
				delete pProbMaxNode;
				pProbMaxNode = 0;
			}
			/* release the combination set */
        		for(ncomb = 0; ncomb < ncomblist; ncomb++) {
          			if(pcomblist[ncomb])
            				CATNET_FREE(pcomblist[ncomb]);
          			pcomblist[ncomb] = NULL;
			} /* for ncomb */
        		CATNET_FREE(pcomblist);
        		pcomblist = 0;
			ncomblist = 0;
		} /* for d */

		/* merge m_pNextNets and pCurCatnetList */		
		for(j = 0; j < m_nNets; j++) {
			if(m_pNextNets[j]) {
				if(pCurCatnetList[j]) {
					if(m_pNextNets[j]->loglik() < pCurCatnetList[j]->loglik()) {
						m_pNextNets[j]->releaseRef();
						m_pNextNets[j] = pCurCatnetList[j];
						pCurCatnetList[j] = 0;
					}
					else {
						pCurCatnetList[j]->releaseRef();
						pCurCatnetList[j] = 0;
					}
				}
			}
			else {
				m_pNextNets[j] = pCurCatnetList[j];
				pCurCatnetList[j] = 0;
			}
		}

	} /* nnode */

	logdiff = FLT_MAX;
	if(m_pNextNets[curSelComplexity])
		logdiff = m_pNextNets[curSelComplexity]->loglik() - curSelLoglik;

	i = 0;
	for(j = 0; j < m_nNets; j++)  {

		if(m_pCurNets[j])
			m_pCurNets[j]->releaseRef();
		m_pCurNets[j] = 0;

		m_emLoglik[j + m_nNets*emiter] = 0;
		if(m_pNextNets[j]) {
			i++;
			m_emLoglik[j + m_nNets*emiter] = m_pNextNets[j]->loglik();
			//printf("pNextNets[%d]  = %p\n", j, m_pNextNets[j]);
		}
		m_pCurNets[j] = m_pNextNets[j];
		m_pNextNets[j] = 0;
	} 

	emiter++;

	if(becho)
		printf("EM Iteration %d: #nets = %d\n", emiter, i);

	/* a local maximum is achieved */
	if(logdiff < stopDelta) {
		emIterations = emiter;
		break;
	}

	} /* while(emiter < emIterations) */

	for(i = 0; i < numNodes; i++)
		CATNET_FREE(pBetasNext[i]);
	CATNET_FREE(pBetasNext);
	CATNET_FREE(pNodeBetas);
	CATNET_FREE(pSigmasNext);

	CATNET_FREE(pCurCatnetList);
	CATNET_FREE(parset);
	CATNET_FREE(fixparset);
	CATNET_FREE(parcats);

	if(psubsamples)
		CATNET_FREE(psubsamples);

	switch(netSelection) {
	case 1: // AIC
		j = 0;
		fLogLik = -FLT_MAX;
		for(nCurNet = 0; nCurNet < m_nNets; nCurNet++) {
			if(!m_pCurNets[nCurNet])
				continue;
			// net complexity = nCurNet + numNodes*(1 + maxCategories)
			complx = m_pCurNets[nCurNet]->getLoglik() - m_pCurNets[nCurNet]->complexity();
			if(complx > fLogLik) {
				fLogLik = complx;
				j = nCurNet;
			}
		}
		break;
	case 2:  // BIC
		j = 0;
		fLogLik = -FLT_MAX;
		ftemp = 0.5*log(numSamples);
		for(nCurNet = 0; nCurNet < m_nNets; nCurNet++) {
			if(!m_pCurNets[nCurNet])
				continue;
			complx = m_pCurNets[nCurNet]->getLoglik() - ftemp*m_pCurNets[nCurNet]->complexity();
			if(complx > fLogLik) {
				fLogLik = complx;
				j = nCurNet;
			}
		}
		break;
	default:
		// optimize w.r.t. most complex network
		j = 0;
		for(nCurNet = 0; nCurNet < m_nNets; nCurNet++)
			if(m_pCurNets[nCurNet]) {
				j = nCurNet;
			}
	}

	if(j < m_nNets) {
		pBetasNext = (double**)m_pCurNets[j]->betas();
		pSigmasNext = (double*)m_pCurNets[j]->sigmas();
		if(becho) {
			printf("Loglik sequence for network with complexity %d: ", j);
			for(emiter = 0; emiter < emIterations; emiter++)
				printf("  %.4f", m_emLoglik[j + m_nNets*emiter]);
			printf("\n");
			for(i = 0; i < numNodes; i++) {
				printf("  beta[%d] = ", i+1);
				for(ic = 0; ic < pNodeNumCats[i]; ic++)
					printf("%.4f  ", pBetasNext[i][ic]);
				printf("\n");
				printf("    sigma[%d] = %.4f\n", i+1, pSigmasNext[i]);
			}
			printf("\n");
		}
	}

	return ERR_CATNET_OK;
}

};

#endif /* EMSEARCH_H */
