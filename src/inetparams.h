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
 * netparams.h
 *
 *  Created on: Nov 16, 2009
 *      Author: Nikolay Balov
 */

#ifndef I_NETPARAMS_H_
#define I_NETPARAMS_H_

#define ERR_CATNET_OK		0
#define ERR_CATNET_PARAM	-1
#define ERR_CATNET_MEM		-2
#define ERR_CATNET_INIT		-3
#define ERR_CATNET_PROC		-4

#include "problist.h" 

template<class t_prob>
class I_NETPARAMS {
public:
	virtual int addRef() = 0;
	virtual int releaseRef() = 0;

	virtual I_NETPARAMS<t_prob>* assign(const I_NETPARAMS<t_prob> *pinet) = 0;
	virtual I_NETPARAMS<t_prob>* clone() = 0;

	virtual int numNodes() = 0;
	virtual int maxParentSet() = 0;
	virtual int maxCategories() = 0;
	virtual const int* numCategories() = 0; 
	virtual int setParents(int nnode, int* parents, int numparents) = 0;
	virtual int getNumParents(int nnode) = 0;
	virtual const PROB_LIST<t_prob>* getNodeProb(int nnode) = 0;
	virtual int setNodeProb(int nnode, PROB_LIST<t_prob> *pprob) = 0;
	virtual t_prob setNodeSampleProb(int nnode, int *psamples, int nsamples) = 0;
	virtual int complexity() = 0; 
	virtual t_prob setLoglik(t_prob lik) = 0;
	virtual t_prob getLoglik() = 0;
	virtual t_prob loglik() = 0;
	virtual t_prob setNodeLoglik(int nnode, t_prob lik) = 0;
	virtual t_prob getNodeLoglik(int nnode) = 0;
	virtual t_prob findNodeLoglik(int nnode, t_prob *psamples, int nsamples) = 0;
	virtual void setNodesOrder(const int *porder) = 0;
	virtual int * getOrder() = 0;
	virtual int* catnetSample(int nsamples) = 0;
	virtual void set_sample_cache(int becho = 0) = 0;

	virtual void setNodeBetas(int nnode, t_prob *pbetas) = 0;
	virtual void setNodeSigma(int nnode, t_prob fSigma) = 0;
	virtual const t_prob **setBetas(t_prob **pbetas, int nbetas) = 0;
	virtual const t_prob *setSigmas(t_prob *psigmas, int nsigmas) = 0;
	virtual const t_prob **betas() = 0;
	virtual const t_prob *sigmas() = 0;

	virtual int pc_size() = 0;
	virtual int pcC_size() = 0;
	virtual const t_prob *pc() = 0;
	virtual const t_prob *pcC() = 0;
	virtual const t_prob *qc() = 0;
	virtual const t_prob *qcC() = 0;
	virtual void release_pqcC() = 0;

	virtual int findNodeMarginalProb(int nnode, t_prob *psamples, int nsamples) = 0;
	virtual int findNodeJointProb(int nnode, int *pnodes, int numnodes, t_prob *psamples, int nsamples) = 0;
	virtual double estimateParameters(int nnode, t_prob *psamples, int nsamples, t_prob *pBetas, t_prob *pSigma) = 0;

	virtual int sample(t_prob *psamples, int nsamples) = 0;
	virtual int predict(t_prob *psamples, int nsamples) = 0;
	virtual int setProbability(t_prob *psamples, int nsamples) = 0;
};


#endif /* I_NETPARAMS_H_ */
