/*
 * catnet.h
 *
 *  Created on: Nov 16, 2009
 *      Author: Nikolay Balov
 */

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <float.h>

#include "problist.h"

#ifndef CATNET_CLASS_H_
#define CATNET_CLASS_H_

//#define DEBUG_ON	1 

#ifdef DEBUG_ON 
#include <iostream.h>
#endif

#define ERR_CATNET_OK		0
#define ERR_CATNET_PARAM	-1
#define ERR_CATNET_MEM		-2
#define ERR_CATNET_INIT		-3
#define ERR_CATNET_PROC		-4

template<class t_node, int t_node_size, class t_prob>
class CATNET {
protected:
	/* nodes are assumed ordered */
	int m_numNodes;
	t_node **m_nodeNames;
	int m_maxParents;
	int *m_numParents;
	int **m_parents;
	int m_maxCategories;
	int *m_numCategories;
	int **m_catIndices;
	int m_complexity;
	t_prob m_loglik;
	PROB_LIST<t_prob> **m_pProbLists;

	int m_nBlockSizes;
	int *m_pBlockSizes;
	t_prob *m_jointProb;
	int m_jointProbSize;
	int m_jointPoolSize;
	int *m_jointPool;
	t_prob *m_margProb;
	int m_margProbSize;
	int *m_margIndex;
	int m_margIndexSize;

public:
	CATNET() {
		_reset();
	}

	CATNET(int nnodes, int maxpars, int maxcats = 2, const t_node **nodes = 0,
			const int * pnumpars = 0, const int **ppars = 0, const int *pcats = 0) {
		_reset();
		init(nnodes, maxpars, maxcats, nodes, pnumpars, ppars, pcats);
	}

	~CATNET() {
		_release();
	}

	CATNET<t_node, t_node_size, t_prob>& operator =(const CATNET<t_node,
			t_node_size, t_prob> &cnet) {
		init(cnet.m_numNodes, cnet.m_maxParents, cnet.m_maxCategories,
				(const t_node**)cnet.m_nodeNames, (const int*)cnet.m_numParents, 
				(const int**)cnet.m_parents,
				(const int*)cnet.m_numCategories, (const PROB_LIST<t_prob> **)cnet.m_pProbLists, 
				cnet.m_complexity, cnet.m_loglik);
		return *this;
	}

protected:
	virtual void _release() {
		int i;
		for (i = 0; i < m_numNodes; i++) {
			if (m_pProbLists && m_pProbLists[i]) {
				delete m_pProbLists[i];
				m_pProbLists[i] = 0;
			}
			if (m_parents && m_parents[i]) {
				CATNET_FREE(m_parents[i]);
				m_parents[i] = 0;
			}
			if (m_nodeNames && m_nodeNames[i]) {
				CATNET_FREE(m_nodeNames[i]);
				m_nodeNames[i] = 0;
			}
			if(m_catIndices && m_catIndices[i]) {
				CATNET_FREE(m_catIndices[i]);
				m_catIndices[i] = 0;
			}
		}
		if (m_numParents)
			CATNET_FREE(m_numParents);
		if (m_parents)
			CATNET_FREE(m_parents);
		if (m_numCategories)
			CATNET_FREE(m_numCategories);
		if (m_nodeNames)
			CATNET_FREE(m_nodeNames);
		if(m_catIndices)
			CATNET_FREE(m_catIndices);
		if (m_pProbLists)
			CATNET_FREE(m_pProbLists);
//printf("CATNET::_release: m_jointProb=%p, m_margProb=%p\n", m_jointProb, m_margProb);
		if(m_jointPool)
			CATNET_FREE(m_jointPool);
		if(m_pBlockSizes)
			CATNET_FREE(m_pBlockSizes);
		if(m_jointProb)
			CATNET_FREE(m_jointProb);
		if(m_margProb)
			CATNET_FREE(m_margProb);
		if(m_margIndex)
			CATNET_FREE(m_margIndex);

		_reset();
	}

	/* do not make this virtual ever */
	void _reset() {
		m_numNodes = 0;
		m_maxParents = 0;
		m_maxCategories = 0;
		m_nodeNames = 0;
		m_numParents = 0;
		m_parents = 0;
		m_numCategories = 0;
		m_catIndices = 0;
		m_pProbLists = 0;
		m_complexity = 0;
		m_loglik = 0;

		m_jointPoolSize = 0;
		m_jointPool = 0;
		m_nBlockSizes = 0;
		m_pBlockSizes = 0;
		m_jointProbSize = 0;
		m_jointProb = 0;
		m_margProbSize = 0;
		m_margProb = 0;
		m_margIndexSize = 0;
		m_margIndex = 0;
	}

public:
	void init(int nnodes, int maxpars, int maxcats = 2, const t_node **nodes = 0,
			const int *pnumpars = 0, const int **ppars = 0, const int *pcats = 0, const PROB_LIST<t_prob> **pprobs = 0, int complexity = 0, t_prob loglik = 0) {

		_release();

		int i, j, *nodeparcats;

		if (nnodes < 1 || maxpars < 0)
			return;
		m_numNodes = nnodes;
		m_maxParents = maxpars;
		m_maxCategories = maxcats;

		m_numParents = (int*) CATNET_MALLOC(m_numNodes * sizeof(int));
		m_parents = (int**) CATNET_MALLOC(m_numNodes * sizeof(int*));
		m_numCategories = (int*) CATNET_MALLOC(m_numNodes * sizeof(int));
		m_catIndices = (int**) CATNET_MALLOC(m_numNodes * sizeof(int*));
		for (i = 0; i < m_numNodes; i++) {
			m_numParents[i] = 0;
			m_parents[i] = 0;
			m_numCategories[i] = m_maxCategories;
			m_catIndices[i] = 0;
		}
		m_nodeNames = (t_node**) CATNET_MALLOC(m_numNodes * sizeof(t_node*));
		if (nodes) {	
			for (i = 0; i < m_numNodes; i++) {
				m_nodeNames[i] = 0;
				if (!nodes[i])
					continue;
				m_nodeNames[i] = (t_node*) CATNET_MALLOC(t_node_size * sizeof(t_node));
				if(strlen(nodes[i]) < t_node_size)
					strcpy(m_nodeNames[i], nodes[i]);
				else {
					memcpy(m_nodeNames[i], nodes[i], (t_node_size-1) * sizeof(t_node));
					m_nodeNames[i][t_node_size-1] = 0;
				}
			}
		}
		else {
			memset(m_nodeNames, 0, m_numNodes * sizeof(t_node*));
		}

		if (pnumpars > 0 && ppars != 0) {
			memcpy(m_numParents, pnumpars, m_numNodes * sizeof(int));
			for (i = 0; i < m_numNodes; i++) {
				m_parents[i] = (int*) CATNET_MALLOC(m_numParents[i] * sizeof(int));
				memset(m_parents[i], 0, m_numParents[i] * sizeof(int));
				if(ppars[i])
					memcpy(m_parents[i], ppars[i], m_numParents[i] * sizeof(int));
			}
		}
		if (pcats != 0)
			memcpy(m_numCategories, pcats, m_numNodes * sizeof(int));
		
		nodeparcats = (int*)CATNET_MALLOC((m_maxParents+1)*sizeof(int));

		m_pProbLists = (PROB_LIST<t_prob>**) CATNET_MALLOC(m_numNodes * sizeof(PROB_LIST<t_prob>*));
		memset(m_pProbLists, 0, m_numNodes * sizeof(PROB_LIST<t_prob>*));
		for (i = 0; i < m_numNodes; i++) {
			if (pprobs && pprobs[i]) {
				m_pProbLists[i] = new PROB_LIST<t_prob>;
				*m_pProbLists[i] = *pprobs[i];
			}
			else {
				for(j = 0; j < m_numParents[i]; j++)
					nodeparcats[j] = m_numCategories[m_parents[i][j]]; 
				m_pProbLists[i] = new PROB_LIST<t_prob>(
					m_numCategories[i], m_maxCategories, m_numParents[i], nodeparcats);
			}
		}

		CATNET_FREE(nodeparcats);

		m_complexity = complexity;
		m_loglik = loglik;
	}

	void setNodesOrder(const int *porder) {
		int i;
		char str[256];

		// reset node names only
		if (!porder) 
			return;
		if(!m_nodeNames)
			m_nodeNames = (t_node**) CATNET_MALLOC(m_numNodes * sizeof(t_node*));	
		for (i = 0; i < m_numNodes; i++) {
			m_nodeNames[i] = 0;
			if (porder[i] < 1 || porder[i] > m_numNodes)
				break;
			sprintf(str, "N%d", (int)porder[i]);
			if(t_node_size <= strlen(str))
				m_nodeNames[i] = (t_node*) CATNET_MALLOC((strlen(str)+1) * sizeof(t_node));
			else
				m_nodeNames[i] = (t_node*) CATNET_MALLOC(t_node_size * sizeof(t_node));
			strcpy((char*)m_nodeNames[i], str);
		}
	}

	void setCategoryIndices(int *pNumCats, int **pCatIndices) {
		int i;
		// reset node names only
		if (!pNumCats || !pCatIndices || !m_numCategories) 
			return;
		//if(m_catIndices)
		//	CATNET_FREE(m_catIndices);
		if(!m_catIndices)
			m_catIndices = (int**) CATNET_MALLOC(m_numNodes * sizeof(int*));
		for (i = 0; i < m_numNodes; i++) {
			if(pNumCats[i] != m_numCategories[i])
				break;
			if(!m_catIndices[i])
				m_catIndices[i] = (int*) CATNET_MALLOC(m_numCategories[i] * sizeof(int));
			memcpy(m_catIndices[i], pCatIndices[i], m_numCategories[i] * sizeof(int));
		}
	}

	void setCondProb(int node, t_prob *pcondprob, int condprobsize) {
		if (m_numNodes < 1 || node < 0 || node >= m_numNodes)
			return;
		if (!m_pProbLists)
			m_pProbLists = (PROB_LIST<t_prob>**) CATNET_MALLOC(m_numNodes
					* sizeof(PROB_LIST<t_prob>*));
		if (m_pProbLists && m_pProbLists[node])
			delete m_pProbLists[node];
		m_pProbLists[node] = 0;
		if (m_numParents[node] < 0 || m_numParents[node] > m_maxParents)
			return;
		int *parcats = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		for (int i = 0; i < m_numParents[node]; i++) {
			parcats[i] = m_numCategories[m_parents[node][i]];
			//cout << "parcats[i] = " << parcats[i] << "\n";
		}
		m_pProbLists[node] = new PROB_LIST<t_prob> (m_numCategories[node],
				m_maxCategories, m_numParents[node], parcats, pcondprob,
				condprobsize);
		CATNET_FREE(parcats);
	}

	int numNodes() {
		return m_numNodes;
	}

	const t_node** nodeNames() {
		return (const t_node**)m_nodeNames;
	}

	int maxCategories() {
		return m_maxCategories;
	}

	const int* numCategories() {
		return (const int*)m_numCategories;
	}

	int numCategories(int node) {
		if (node < 0 || node >= m_numNodes)
			return 0;
		return m_numCategories[node];
	}

	int maxParents() {
		return m_maxParents;
	}

	const int** parents() {
		return (const int**)m_parents;
	}

	const int* numParents() {
		return (const int*)m_numParents;
	}

	int numParents(int node) {
		if (node < 0 || node >= m_numNodes)
			return 0;
		return m_numParents[node];
	}

	const int* getParents(int node) {
		if (node < 0 || node >= m_numNodes)
			return 0;
		return (const int*)m_parents[node];
	}

	int setParents(int node, int* parents, int numparents) {
		if (node < 0 || node >= m_numNodes)
			return 0;
		if (!m_numParents || !m_parents)
			return 0;
		if(m_numParents[node] != numparents) {
			m_numParents[node] = numparents;
			if(m_parents[node])
				CATNET_FREE(m_parents[node]);
			m_parents[node] = (int*) CATNET_MALLOC(m_numParents[node] * sizeof(int));
		}
		memcpy(m_parents[node], parents, m_numParents[node] * sizeof(int));

		if (!m_pProbLists)
			m_pProbLists = (PROB_LIST<t_prob>**) CATNET_MALLOC(m_numNodes * sizeof(PROB_LIST<t_prob>*));
		else {
			if(m_pProbLists[node])
				delete m_pProbLists[node];
		}
		int *parcats = (int*)CATNET_MALLOC(m_numParents[node]*sizeof(int));
		for(int i = 0; i < m_numParents[node]; i++) 
			parcats[i] = m_numCategories[parents[i]];	
		m_pProbLists[node] = new PROB_LIST<t_prob>(m_numCategories[node], m_maxCategories, m_numParents[node], parcats);
		CATNET_FREE(parcats);
		return m_numParents[node];
	}

	int complexity() {
		int i, j, c;
		m_complexity = 0;
		for (i = 0; i < m_numNodes; i++) {
			if (!m_parents || !m_parents[i]) {
				m_complexity += (m_numCategories[i]-1);
				continue;
			}
			c = 1;
			for (j = 0; j < m_numParents[i]; j++)
				c *= m_numCategories[m_parents[i][j]];
			m_complexity += c*(m_numCategories[i]-1);
		}
		return m_complexity;
	}

	t_prob getLoglik() {
		return m_loglik;
	}

	t_prob setLoglik(t_prob lik) {
		m_loglik = lik;
		return m_loglik;
	}

	t_prob loglik() {
		int i;
		m_loglik = 0;
		for (i = 0; i < m_numNodes; i++) {
			if(m_pProbLists[i])
				m_loglik += m_pProbLists[i]->loglik;
		}
		return m_loglik;
	}

	const PROB_LIST<t_prob>** probLists() {
		return (const PROB_LIST<t_prob>**)m_pProbLists;
	}

	const PROB_LIST<t_prob>* getNodeProb(int nnode) {
		if (!m_pProbLists || nnode < 0 || nnode >= m_numNodes)
			return 0;
		return m_pProbLists[nnode];
	}

	int setNodeProb(int nnode, PROB_LIST<t_prob> *pprob) {
		if(!m_pProbLists || nnode < 0 || nnode >= m_numNodes)
			return 0;
		if(!m_pProbLists[nnode])
			m_pProbLists[nnode] = new PROB_LIST<t_prob>;
		if (m_pProbLists[nnode] && pprob) {
			*m_pProbLists[nnode] = *pprob;
			// normalize
			//m_pProbLists[nnode]->normalize();
		}
		return 0;
	}

	t_prob *nodeSampleLoglik(int nnode, int *pnodepars, int nodepars,
			int *psamples, int nsamples) {
		int j, ipar;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, floglik = 0;
		int *pnodesample, samp;

		if(!psamples || nsamples < 1)
			return 0;

		pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		pnodepars = m_parents[nnode];
		for (j = 0; j < nsamples; j++) {
			for (ipar = 0; ipar < nodepars; ipar++) {
				pnodesample[ipar] = psamples[j * m_numNodes + pnodepars[ipar]];
			}
			pnodeprob = m_pProbLists[nnode]->find_slot(0, pnodesample, 0);
			samp = psamples[j * m_numNodes + nnode];
			if (samp >= 0 && samp < m_numCategories[nnode] && pnodeprob[samp] > 0)
				floglik += log(pnodeprob[samp]);
		}
		CATNET_FREE(pnodesample);
		return floglik;
	}

	t_prob *sampleLoglik(int *psamples, int nsamples) {
		int i, j, ipar;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob;
		int nodepars;
		int *pnodepars, *pnodesample, samp;

		if(!psamples || nsamples < 1)
			return 0;

		pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		m_loglik = 0;
		for (i = 0; i < m_numNodes; i++) {
			if(!m_pProbLists[i])
				continue;
			pnodepars = m_parents[i];
			nodepars = m_numParents[i];
			for (j = 0; j < nsamples; j++) {
				for (ipar = 0; ipar < nodepars; ipar++) {
					pnodesample[ipar] = psamples[j * m_numNodes + pnodepars[ipar]];
				}
				pnodeprob = m_pProbLists[i]->find_slot(0, pnodesample, 0);
				samp = psamples[j * m_numNodes + i];
				if (pnodeprob && samp >= 0 && samp < m_numCategories[i] && pnodeprob[samp] > 0)
					m_loglik += log(pnodeprob[samp]);
			}
		}
		CATNET_FREE(pnodesample);
		return m_loglik;
	}

	// sets sample conditional probability and returns its log-likelihood
	t_prob setNodeSampleProb(int nnode, int *psamples, int nsamples, int bNormalize = 0) {
		int i, j;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, floglik;
		int *pnodesample, *pnodepars, samp;
		pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));

		if (!m_pProbLists || !psamples || nsamples < 1)
			return 0;

		if(!m_pProbLists[nnode]) {
			// error
			return 0;
		}
		else
			m_pProbLists[nnode]->set_zero();
		pnodepars = m_parents[nnode];

		for (j = 0; j < nsamples; j++) {
			for (i = 0; i < m_numParents[nnode]; i++) {
				if (pnodepars[i] < 0 || pnodepars[i] >= m_numNodes)
					break;
				pnodesample[i] = psamples[j * m_numNodes + pnodepars[i]];
			}
			pnodeprob = m_pProbLists[nnode]->find_slot(0, pnodesample, 0);
			samp = psamples[j * m_numNodes + nnode];
			if (pnodeprob && samp >= 0 && samp < m_numCategories[nnode])
				pnodeprob[samp]++;
		}

		CATNET_FREE(pnodesample);

		// at this point m_pProbLists[nnode] has counts not probabilities 
		floglik = m_pProbLists[nnode]->loglikelihood();

		if(bNormalize)
			m_pProbLists[nnode] -> normalize();

		return(floglik);
	}

	int normalizeProbabilities() {
		int i;
		if(!m_pProbLists)
			return -1;
		for (i = 0; i < m_numNodes; i++) {
			if(m_pProbLists[i])
				m_pProbLists[i] -> normalize();
		}
		return 0;
	}

	int *getOrder() {
		int *porder = 0, *pfound = 0, bhaspar, i, j, k;
		if(m_numNodes < 1 || !m_numParents || !m_parents)
			return 0;

		porder = (int*)CATNET_MALLOC(m_numNodes * sizeof(int));
		pfound = (int*)CATNET_MALLOC(m_numNodes * sizeof(int));
		memset(pfound, 0, m_numNodes * sizeof(int));
		for(i = 0; i < m_numNodes; i++) {
			for(j = 0; j < m_numNodes; j++) {
				if(pfound[j])
					continue;
				if(m_numParents[j] <= 0)
					break;
				bhaspar = 0;
				for(k = 0; k < m_numParents[j]; k++) {
					if(pfound[m_parents[j][k]])
						continue;
					bhaspar = 1;
					break;
				}
				if(!bhaspar)
					break;
			}
			if(j >= m_numNodes || pfound[j]) {
				CATNET_FREE(pfound);
				CATNET_FREE(porder);
				return 0;
			}
			porder[i] = j;
			pfound[j] = 1;
			//printf("porder[%d] = %d\n", i+1, porder[i]+1);
		}
		CATNET_FREE(pfound);
		return porder;
	}

	/*****************************************
	  Joint and Marginal Probability Routines
	******************************************/

	int *parPool() {
		return m_jointPool;
	}

	int parPoolSize() {
		return m_jointPoolSize;
	}

	t_prob *jointProb() {
		return m_jointProb;
	}

	int jointProbSize() {
		return m_jointProbSize;
	}

	t_prob *margProb() {
		return m_margProb;
	}

	int margProbSize() {
		return m_margProbSize;
	}

private:
	int *_findParentPool(int nnode, int &poolsize) {

		int ipar, par, i, j, bfound, parpoolsize, *parpool, *ppool, *paux;

		poolsize = 0;
		if (nnode < 0 || nnode >= m_numNodes || !m_parents || !m_parents[nnode] || !m_numParents[nnode])
			return 0;

		//cout << "nnode = " << nnode+1 << ", m_numParents[nnode] = " << m_numParents[nnode] << "\n";

		ppool = 0;
		for (ipar = 0; ipar < m_numParents[nnode]; ipar++) {
			par = m_parents[nnode][ipar];
			parpool = _findParentPool(par, parpoolsize);

			i = 0;
			while (i < parpoolsize) {
				bfound = 0;
				for (j = 0; j < poolsize; j++) {
					if (parpool[i] == ppool[j]) {
						bfound = 1;
						break;
					}
				}
				if (!bfound) {
					i++;
					continue;
				}
				for (j = i + 1; j < parpoolsize; j++)
					parpool[j - 1] = parpool[j];
				parpoolsize--;
			}

			//cout << "par = " << par+1 << "poolsize = " << poolsize << ", parpoolsize = " << parpoolsize << "\n";
			paux = (int*) CATNET_MALLOC((poolsize + parpoolsize + 1) * sizeof(int));
			if (poolsize > 0 && ppool)
				memcpy(paux, ppool, poolsize * sizeof(int));
			if (parpoolsize > 0 && parpool)
				memcpy(paux + poolsize, parpool, parpoolsize * sizeof(int));
			// release that
			if(parpool)
				CATNET_FREE(parpool);

			// check par
			bfound = 0;
			for (j = 0; j < poolsize + parpoolsize; j++) {
				if (paux[j] == par) {
					bfound = 1;
					break;
				}
			}
			if (bfound) {
				poolsize += parpoolsize;
			} else {
				paux[poolsize + parpoolsize] = par;
				poolsize += (parpoolsize + 1);
			}

			if (ppool)
				CATNET_FREE(ppool);
			ppool = paux;
		}
		return ppool;
	}

	int *_findParentPool(int *pnodes, int numnodes, int &poolsize) {
		int inode, nnode, ipar, par, i, j, bfound, parpoolsize, *parpool, *ppool, *paux;

		ppool = 0;
		poolsize = 0;
		if (numnodes < 0 || !m_parents)
			return 0;

		for(inode = 0; inode < numnodes; inode++) { 

		nnode = pnodes[inode]; 

		if(nnode >= m_numNodes || !m_parents[nnode] || !m_numParents[nnode])
			continue;

		for (ipar = 0; ipar < m_numParents[nnode]; ipar++) {
			par = m_parents[nnode][ipar];
			parpool = _findParentPool(par, parpoolsize);

//printf("\n\nppool: ");
//for (j = 0; j < poolsize; j++) 
//	printf("%d, ", ppool[j]+1);
//printf("\n");
//printf("parpool: ");
//for (j = 0; j < parpoolsize; j++) 
//	printf("%d, ", parpool[j]+1);
//printf("\n");
			i = 0;
			while (i < parpoolsize) {
				bfound = 0;
				for (j = 0; j < poolsize; j++) {
					if (parpool[i] == ppool[j]) {
						bfound = 1;
						break;
					}
				}
				if (!bfound) {
					i++;
					continue;
				}
				for (j = i + 1; j < parpoolsize; j++)
					parpool[j - 1] = parpool[j];
				parpoolsize--;
			}

			//cout << "par = " << par+1 << "poolsize = " << poolsize << ", parpoolsize = " << parpoolsize << "\n";
			paux = (int*) CATNET_MALLOC((poolsize + parpoolsize + 1) * sizeof(int));
			if (poolsize > 0 && ppool)
				memcpy(paux, ppool, poolsize * sizeof(int));
			if (parpoolsize > 0 && parpool)
				memcpy(paux + poolsize, parpool, parpoolsize * sizeof(int));
			// release that
			if(parpool)
				CATNET_FREE(parpool);

			// check par
			bfound = 0;
			for (j = 0; j < poolsize + parpoolsize; j++) {
				if (paux[j] == par) {
					bfound = 1;
					break;
				}
			}
			if (bfound) {
				poolsize += parpoolsize;
			} else {
				paux[poolsize + parpoolsize] = par;
				poolsize += (parpoolsize + 1);
			}

			if (ppool)
				CATNET_FREE(ppool);
			ppool = paux;
		}

		} // inode

		return ppool;
	}


public:
	int findParentPool(int nnode) {
		m_jointPoolSize = 0;
		if(m_jointPool)
			CATNET_FREE(m_jointPool);
		m_jointPool = _findParentPool(nnode, m_jointPoolSize);
		return (m_jointPoolSize == 0)?ERR_CATNET_MEM:ERR_CATNET_OK;
	}

	int findParentPool(int *pnodes, int numnodes) {
		m_jointPoolSize = 0;
		if(m_jointPool)
			CATNET_FREE(m_jointPool);
		m_jointPool = _findParentPool(pnodes, numnodes, m_jointPoolSize);
		return (m_jointPoolSize == 0)?ERR_CATNET_MEM:ERR_CATNET_OK;
	}

	int findJointProb(int nnode) {
		t_prob *parprob;
		int i, ii, i0, ic, ipar, j, par, ipool, poolnode;
		int *paux, *paridx, *pcats;
		PROB_LIST<t_prob> *probnode = 0;

		if (nnode < 0 || nnode >= m_numNodes)
			return ERR_CATNET_PARAM;
		if (!m_pProbLists)
			return ERR_CATNET_INIT;

		m_jointPoolSize = 0;	
		if(m_jointPool)
			CATNET_FREE(m_jointPool);
		m_jointPool = _findParentPool(nnode, m_jointPoolSize);

		// add nnode to the parents list
		paux = (int*) CATNET_MALLOC((m_jointPoolSize + 1) * sizeof(int));
		if (m_jointPool && m_jointPoolSize > 0)
			memcpy(paux, m_jointPool, m_jointPoolSize * sizeof(int));
		if(m_jointPool);
			CATNET_FREE(m_jointPool);
		paux[m_jointPoolSize] = nnode;
		m_jointPoolSize++;
		m_jointPool = paux;

		pcats = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		paridx = (int*) CATNET_MALLOC(m_jointPoolSize * sizeof(int));

		if(m_pBlockSizes)
			CATNET_FREE(m_pBlockSizes);
		m_nBlockSizes = m_jointPoolSize;
		m_pBlockSizes = (int*) CATNET_MALLOC(m_jointPoolSize * sizeof(int));
		m_pBlockSizes[m_jointPoolSize - 1] = 1;
		for (i = m_jointPoolSize - 2; i >= 0; i--) {
			m_pBlockSizes[i] = m_pBlockSizes[i + 1] * m_numCategories[m_jointPool[i + 1]];
		}

		m_jointProbSize = 1;
		for (i = 0; i < m_jointPoolSize; i++) {
			m_jointProbSize *= m_numCategories[m_jointPool[i]];
		}

		if(m_jointProb)
			CATNET_FREE(m_jointProb);
		m_jointProb = (t_prob*) CATNET_MALLOC(m_jointProbSize * sizeof(t_prob));
		for (i = 0; i < m_jointProbSize; i++)
			m_jointProb[i] = 1.0;

		for (ipool = 0; ipool < m_jointPoolSize; ipool++) {
			poolnode = m_jointPool[ipool];
			probnode = m_pProbLists[poolnode];

			//cout << ipool << ", poolnode = " << poolnode + 1 << "\n";

			if (m_numParents[poolnode] == 0) {
				for (ii = 0; ii < m_jointProbSize; ii += (m_pBlockSizes[ipool]
						* m_numCategories[poolnode])) {
					for (ic = 0; ic < m_numCategories[poolnode]; ic++) {
						i0 = ic * m_pBlockSizes[ipool];
						for (i = 0; i < m_pBlockSizes[ipool]; i++) {
							m_jointProb[ii + i0 + i] *= probnode->pProbs[ic];
						}
					}
				}
				continue;
			}

			//cout << "Parents: ";
			memset(paridx, 0, m_jointPoolSize * sizeof(int));
			for (ipar = 0; ipar < m_numParents[poolnode]; ipar++) {
				par = m_parents[poolnode][ipar];
				for (i = 0; i < ipool; i++)
					if (par == m_jointPool[i])
						break;
				if (i >= ipool) {
					// bad
					break;
				}
				paridx[i] = ipar + 1;
			}

			memset(pcats, 0, m_maxParents * sizeof(int));
			for (ii = 0; ii < m_jointProbSize; ii += (m_pBlockSizes[ipool] * m_numCategories[poolnode])) {
				for (j = 0; j < ipool; j++) {
					if (paridx[j] > 0) {
						ic = (int) (ii / m_pBlockSizes[j]);
						ic -= m_numCategories[poolnode] * (int) (ic / m_numCategories[poolnode]);
						pcats[paridx[j] - 1] = ic;
					}
				}
				parprob = probnode->find_slot(0, pcats, 0);
				if (!m_jointPool)
					continue;
				for (ic = 0; ic < m_numCategories[poolnode]; ic++) {
					i0 = ic * m_pBlockSizes[ipool];
					for (i = 0; i < m_pBlockSizes[ipool]; i++) {
						m_jointProb[ii + i0 + i] *= parprob[ic];
					}
				}
			}
		}

		CATNET_FREE(pcats);
		CATNET_FREE(paridx);

		return ERR_CATNET_OK;
	}

	int marginalProb(int nnode) {
		int nodecats, ic, k;
		nodecats = m_numCategories[nnode];

		if (findJointProb(nnode))
			return ERR_CATNET_MEM;

		if(m_pBlockSizes)
			CATNET_FREE(m_pBlockSizes);
		m_nBlockSizes = 1;
		m_pBlockSizes = (int*) CATNET_MALLOC(m_nBlockSizes * sizeof(int));
		m_pBlockSizes[0] = 1;

		m_margProbSize = nodecats;
		if(m_margProb)
			CATNET_FREE(m_margProb);
		m_margProb = (t_prob*) CATNET_MALLOC(m_margProbSize * sizeof(t_prob));
		for (ic = 0; ic < nodecats; ic++) {
			k = ic;
			m_margProb[ic] = 0;
			while (k < m_jointProbSize) {
				m_margProb[ic] += m_jointProb[k];
				k += nodecats;
			}
		}

		return ERR_CATNET_OK;
	}

	int findJointProb(int *pnodes, int numnodes) {

		t_prob *parprob;
		int inode, nnode;
		int i, ii, i0, ic, ipar, j, par, ipool, poolnode;
		int *paux, *paridx, *pcats;
		PROB_LIST<t_prob> *probnode = 0;

		if (!m_pProbLists)
			return ERR_CATNET_INIT;

		m_jointPoolSize = 0;	
		if(m_jointPool)
			CATNET_FREE(m_jointPool);
		m_jointPool = _findParentPool(pnodes, numnodes, m_jointPoolSize);

		//if(exp(log(m_jointPoolSize+1) * m_maxCategories) > MAX_MEM_PROB) {
		//	// MC simulation of the joint probability
		//}
		//printf("m_jointPoolSize = %d\n", m_jointPoolSize);

		// add nnode to the parents list
		paux = (int*) CATNET_MALLOC((m_jointPoolSize + numnodes) * sizeof(int));
		if (m_jointPool && m_jointPoolSize > 0) 
			memcpy(paux, m_jointPool, m_jointPoolSize * sizeof(int));
		if(m_jointPool)
			CATNET_FREE(m_jointPool);
		m_jointPool = paux;
		for(inode = 0; inode < numnodes; inode++) { 
			nnode = pnodes[inode];
			for (i = 0; i < m_jointPoolSize; i++)
				if(nnode == m_jointPool[i])
					break;
			if(i < m_jointPoolSize)
				continue;
			m_jointPool[m_jointPoolSize] = nnode;
			m_jointPoolSize++;
		}

		pcats = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		paridx = (int*) CATNET_MALLOC(m_jointPoolSize * sizeof(int));

		if(m_pBlockSizes)
			CATNET_FREE(m_pBlockSizes);
		m_nBlockSizes = m_jointPoolSize;
		m_pBlockSizes = (int*) CATNET_MALLOC(m_jointPoolSize * sizeof(int));
		m_pBlockSizes[m_jointPoolSize - 1] = 1;
		for (i = m_jointPoolSize - 2; i >= 0; i--) {
			m_pBlockSizes[i] = m_pBlockSizes[i + 1] * m_numCategories[m_jointPool[i + 1]];
			
		}

		m_jointProbSize = 1;
		for (i = 0; i < m_jointPoolSize; i++) {
			m_jointProbSize *= m_numCategories[m_jointPool[i]];
		}

		if(m_jointProb)
			CATNET_FREE(m_jointProb);
		m_jointProb = (t_prob*) CATNET_MALLOC(m_jointProbSize * sizeof(t_prob));
		for (i = 0; i < m_jointProbSize; i++)
			m_jointProb[i] = 1.0;

		for (ipool = 0; ipool < m_jointPoolSize; ipool++) {
			poolnode = m_jointPool[ipool];
			probnode = m_pProbLists[poolnode];

			//cout << ipool << ", poolnode = " << poolnode + 1 << "\n";

			if (m_numParents[poolnode] == 0) {
				for (ii = 0; ii < m_jointProbSize; ii += (m_pBlockSizes[ipool]
						* m_numCategories[poolnode])) {
					for (ic = 0; ic < m_numCategories[poolnode]; ic++) {
						i0 = ic * m_pBlockSizes[ipool];
						for (i = 0; i < m_pBlockSizes[ipool]; i++) {
							m_jointProb[ii + i0 + i] *= probnode->pProbs[ic];
						}
					}
				}
				continue;
			}

			//cout << "Parents: ";
			memset(paridx, 0, m_jointPoolSize * sizeof(int));
			for (ipar = 0; ipar < m_numParents[poolnode]; ipar++) {
				par = m_parents[poolnode][ipar];
				for (i = 0; i < ipool; i++)
					if (par == m_jointPool[i])
						break;
				if (i >= ipool) {
					// bad
					break;
				}
				paridx[i] = ipar + 1;
				//cout << ipar + 1 << ", ";
			}
			//cout << "\n";

			memset(pcats, 0, m_maxParents * sizeof(int));
			for (ii = 0; ii < m_jointProbSize; ii += (m_pBlockSizes[ipool] * m_numCategories[poolnode])) {
				for (j = 0; j < ipool; j++) {
					if (paridx[j] > 0) {
						ic = (int) (ii / m_pBlockSizes[j]);
						ic -= m_numCategories[poolnode] * (int) (ic / m_numCategories[poolnode]);
						//cout << "paridx[i] = " << paridx[j] << ", ic = " << ic << "\n";
						pcats[paridx[j] - 1] = ic;
					}
				}
				parprob = probnode->find_slot(0, pcats, 0);
				if (!parprob)
					continue;
				//cout << ii << ", " << m_pBlockSizes[ipool] << ", " << m_jointProbSize << "\n";
				for (ic = 0; ic < m_numCategories[poolnode]; ic++) {
					i0 = ic * m_pBlockSizes[ipool];
					for (i = 0; i < m_pBlockSizes[ipool]; i++) {
						m_jointProb[ii + i0 + i] *= parprob[ic];
					}
				}
			}
		}

		CATNET_FREE(pcats);
		CATNET_FREE(paridx);

		return ERR_CATNET_OK;
	}

	int marginalProb(int *pnodes, int numnodes) {
		t_prob *paux;
		int ic, k, kk, i, ii, m, imarg, nnode, newsize;
		
		/* ASSUME THAT pnodes are consistent the topological order of the network, i.e. 
		   no node with higher index is a parent of a node with lower one  */
//printf("FindJointProb ");
		if (findJointProb(pnodes, numnodes) != ERR_CATNET_OK)
			return ERR_CATNET_MEM;

		if(m_margProb)
			CATNET_FREE(m_margProb);
		m_margProbSize = m_jointProbSize;
		m_margProb = (t_prob*) CATNET_MALLOC(m_jointProbSize * sizeof(t_prob));
		memcpy(m_margProb, m_jointProb, m_jointProbSize*sizeof(t_prob));

//printf("numnodes = %d, jointProbSize = %d\n", numnodes, m_jointProbSize);

		paux = (t_prob*) CATNET_MALLOC(m_jointProbSize * sizeof(t_prob));

		if(m_margIndex)
			CATNET_FREE(m_margIndex);
		m_margIndexSize = numnodes;
		m_margIndex = (int*)CATNET_MALLOC(m_margIndexSize*sizeof(int));

		/* contraction */
		/* pnodes has the same relative order in m_jointPool */
		/* work with the m_pBlockSizes for the joint prob */
		imarg = 0;
		for(i = 0; i < m_jointPoolSize; i++) {
			nnode = m_jointPool[i];
			for(k = 0; k < numnodes; k++)
				if(nnode == pnodes[k])
					break;
			if(k < numnodes) {
				m_margIndex[k] = imarg++;
				continue;
			}

			newsize = (int)(m_margProbSize/m_numCategories[nnode]);

//printf("imarg = %d, i = %d, nnode=%d, margProbSize = %d, newsize = %d, m_pBlockSizes[i] = %d \n", 
//imarg, i, nnode+1, m_margProbSize, newsize, m_pBlockSizes[i]);

			ii = 0;
			k = 0;
			while(ii < newsize) {
				for(m = 0; m < m_pBlockSizes[i]; m++) {
					kk = k + m;
					paux[ii] = 0;
					for(ic = 0; ic < m_numCategories[nnode]; ic++) {
//printf("ii = %d, k = %d, kk=%d, ic=%d\n", ii, k, kk, ic);
						paux[ii] += m_margProb[kk];
						kk += m_pBlockSizes[i];
					}
					ii++;
				}
				k += m_pBlockSizes[i]*m_numCategories[nnode];
			}

			m_margProbSize = newsize;
			memcpy(m_margProb, paux, newsize * sizeof(t_prob));
		}

		CATNET_FREE(paux);

		if(m_pBlockSizes)
			CATNET_FREE(m_pBlockSizes);
		m_nBlockSizes = numnodes;
		m_pBlockSizes = (int*) CATNET_MALLOC(m_nBlockSizes * sizeof(int));
		m_pBlockSizes[numnodes - 1] = 1;
		for (i = numnodes - 2; i >= 0; i--) {
			nnode = pnodes[i];
			m_pBlockSizes[i] = m_pBlockSizes[i + 1] * m_numCategories[nnode];
		}

		m_margProbSize = m_pBlockSizes[0] * m_numCategories[pnodes[numnodes-1]];
		paux = (t_prob*) CATNET_MALLOC(m_margProbSize * sizeof(t_prob));
		memcpy(paux, m_margProb, m_margProbSize * sizeof(t_prob));
		CATNET_FREE(m_margProb);
		m_margProb = paux;

//printf("m_jointProb=%p, m_margProb=%p\n", m_jointProb, m_margProb);

//printf("margIndex: ");
//for(k=0;k<m_margIndexSize;k++)
//	printf("%d, ", m_margIndex[k]);
//printf("\n");

		return ERR_CATNET_OK;
	}

	t_prob getCatProb(int ncat) {
		int ipos;
		if(!m_margProb)
			return 0;
		ipos = m_pBlockSizes[0]*ncat;
		if(ipos >= m_margProbSize)
			return 0;
		return m_margProb[ipos];
	}

	t_prob getCatProb(int *pcats, int numnodes) {
		int ipos, i, j;
//printf("numnodes = %d, nBlockSizes = %d, margIndexSize = %d\n", numnodes, m_nBlockSizes, m_margIndexSize);
		if(!m_margProb || numnodes != m_nBlockSizes || numnodes != m_margIndexSize)
			return 0;
		ipos = 0;
		for(i = 0; i < numnodes; i++) {
			j = m_margIndex[i];
			ipos += m_pBlockSizes[j]*pcats[i];
		}
		if(ipos >= m_margProbSize)
			return 0;
		return m_margProb[ipos];
	}

};

#endif /* CATNET_CLASS_H_ */
