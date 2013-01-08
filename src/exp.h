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
 * exp.h
 *
 *  Created on: May 6, 2010
 *      Author: Nikolay Balov
 */

#ifndef EXPNET_CLASS_H
#define EXPNET_CLASS_H

#include "catnet.h"
#include "inetparams.h"

#define MAX_NODE_NAME	16

extern int g_netcounter;

class CExpNet : public CATNET<char, double>, public I_NETPARAMS<double> {
private:
	int m_ref;
	int *m_pCatnetSamples, m_nCatnetSamples;

protected:
	double **m_lambdas, **m_loglambdas;

	int m_curNode;
	int m_parCatSetSize;
	double *m_pc, *m_pcC;
	double *m_qc, *m_qcC;

	double *m_nodeLoglik;

	void _release();
	void _reset();

public:
	CExpNet() {
		_reset();
		m_ref = 0;
		addRef();
	}
	
	CExpNet(int nnodes, int maxpars, int maxcats = 2, const char **nodes = 0,
			const int * pnumpars = 0, const int **ppars = 0, const int *pcats = 0) {
		_reset();
		init(nnodes, maxpars, maxcats, nodes, pnumpars, ppars, pcats);
		m_ref = 0;
		addRef();
	}

	~CExpNet() {
		_release();
	}

	int addRef() {
		m_ref++;
		g_netcounter++;
		return m_ref;
	}

	int releaseRef() {
		g_netcounter--;
		m_ref--;
		if(m_ref <= 0) {
			delete this;
			return 0;
		}
		return m_ref;
	}

	CExpNet& operator =(const CExpNet &cnet);

	/* I_NETPARAMS methods */

	I_NETPARAMS<double> *assign(const I_NETPARAMS<double> *pinet) {
		*this = *((CExpNet*)pinet);
		return (I_NETPARAMS<double>*)this;
	}

	I_NETPARAMS<double> *clone() {
		CExpNet *pinet = new CExpNet; 
		*pinet = *this;
		return (I_NETPARAMS<double>*)pinet;
	}

	int numNodes() {
		return CATNET<char, double>::numNodes();
	}

	int maxParentSet() {
		return CATNET<char, double>::maxParents();
	}

	int maxCategories() {
		return CATNET<char, double>::maxCategories();
	}

	const int* numCategories() {
		return CATNET<char, double>::numCategories();
	} 

	int setParents(int nnode, int* parents, int numparents) {
		return CATNET<char, double>::setParents(nnode, parents, numparents);
	}

	int getNumParents(int nnode) {
		return CATNET<char, double>::numParents(nnode);
	}

	int setNodeProb(int nnode, PROB_LIST<double> *pprob) {
		return CATNET<char, double>::setNodeProb(nnode, pprob);
	}

	const PROB_LIST<double>* getNodeProb(int nnode) {
		return CATNET<char, double>::getNodeProb(nnode);
	}

	double setNodeSampleProb(int nnode, int *psamples, int nsamples) {
		return CATNET<char, double>::setNodeSampleProb(
							nnode, psamples, nsamples, 1/* do normalize */);
	}

	int complexity() {
		int n;
		int numnodes = numNodes();
		const int* pNumCats = numCategories();
		int complx = 0;
		for(n = 0; n < numnodes; n++)
			complx += pNumCats[n];
		return(complx + CATNET<char, double>::complexity());
	}

	double setLoglik(double lik) {
		return CATNET<char, double>::setLoglik(lik);
	}

	double getLoglik() {
		return CATNET<char, double>::getLoglik();
	}

	double loglik();

	double setNodeLoglik(int nnode, double lik);
	double getNodeLoglik(int nnode);
	double findNodeLoglik(int nnode, double *psamples, int nsamples);
	double findLogNodeLikelihood(int nnode, double *psamples, int nsamples);
	
	void setNodesOrder(const int *porder) {
		return CATNET<char, double>::setNodesOrder(porder);
	}

	int *getOrder() {
		return CATNET<char, double>::getOrder();
	}

	int* catnetSample(int nsamples);

	void set_sample_cache(int becho = 0);

	void setNodeBetas(int nnode, double *pbetas);
	void setNodeSigma(int nnode, double fSigma);
	const double **setBetas(double **pbetas, int nbetas);
	const double *setSigmas(double *psigmas, int nsigmas);
	const double **betas();
	const double *sigmas();

	int pc_size();
	int pcC_size();
	const double *pc();
	const double *pcC();
	const double *qc();
	const double *qcC();
	void release_pqcC();

	int findNodeMarginalProb(int nnode, double *psamples, int nsamples);
	int findNodeJointProb(int nnode, int *pnodes, int numnodes, double *psamples, int nsamples);
	int findParentsJointProb(int *pnodes, int numnodes, double *psamples, int nsamples);
	double estimateParameters(int nnode, double *psamples, int nsamples, 
					double *pBetas, double *pSigma);

	int sample(double *psamples, int nsamples);
	int predict(double *psamples, int nsamples);
	
	int estimateNodeProb(int nnode, int *parnodes, int numpars, double *psamples, int nsamples);
	int setProbability(double *psamples, int nsamples);
};

#endif /* EXPNET_CLASS_H */
