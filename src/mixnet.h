/*
 * mixnet.h
 *
 *  Created on: Nov 16, 2009
 *      Author: Nikolay Balov
 */

#ifndef MIXNET_CLASS_H
#define MIXNET_CLASS_H

#include "catnet.h"
#include "inetparams.h"

#define MAX_NODE_NAME	16

extern int g_netcounter;

class CMixNet : public CATNET<char, MAX_NODE_NAME, double>, public I_NETPARAMS<double> {
private:
	int m_ref;
	int *m_pCatnetSamples, m_nCatnetSamples;

protected:
	double **m_betas, *m_sigmas;

	int m_curNode;
	int m_parCatSetSize;
	double *m_pc, *m_pcC;
	double *m_qc, *m_qcC;

	double *m_nodeLoglik;

	void _release();
	void _reset();

public:
	CMixNet() {
		_reset();
		m_ref = 0;
		addRef();
	}
	
	CMixNet(int nnodes, int maxpars, int maxcats = 2, const char **nodes = 0,
			const int * pnumpars = 0, const int **ppars = 0, const int *pcats = 0) {
		_reset();
		init(nnodes, maxpars, maxcats, nodes, pnumpars, ppars, pcats);
		m_ref = 0;
		addRef();
	}

	~CMixNet() {
		_release();
	}

	int addRef() {
		m_ref++;
		g_netcounter++;
		//printf("+netref %d\n", g_netcounter);
		return m_ref;
	}

	int releaseRef() {
		g_netcounter--;
		//printf("-netref %d\n", g_netcounter);
		m_ref--;
		if(m_ref <= 0) {
			delete this;
			return 0;
		}
		return m_ref;
	}

	CMixNet& operator =(const CMixNet &cnet);

	/* I_NETPARAMS methods */

	I_NETPARAMS<double> *assign(const I_NETPARAMS<double> *pinet) {
		*this = *((CMixNet*)pinet);
		return (I_NETPARAMS<double>*)this;
	}

	I_NETPARAMS<double> *clone() {
		CMixNet *pinet = new CMixNet; 
		*pinet = *this;
		return (I_NETPARAMS<double>*)pinet;
	}

	int numNodes() {
		return CATNET<char, MAX_NODE_NAME, double>::numNodes();
	}

	int maxParentSet() {
		return CATNET<char, MAX_NODE_NAME, double>::maxParents();
	}

	int maxCategories() {
		return CATNET<char, MAX_NODE_NAME, double>::maxCategories();
	}

	const int* numCategories() {
		return CATNET<char, MAX_NODE_NAME, double>::numCategories();
	} 

	int setParents(int nnode, int* parents, int numparents) {
		return CATNET<char, MAX_NODE_NAME, double>::setParents(nnode, parents, numparents);
	}

	int getNumParents(int nnode) {
		return CATNET<char, MAX_NODE_NAME, double>::numParents(nnode);
	}

	int setNodeProb(int nnode, PROB_LIST<double> *pprob) {
		return CATNET<char, MAX_NODE_NAME, double>::setNodeProb(nnode, pprob);
	}

	const PROB_LIST<double>* getNodeProb(int nnode) {
		return CATNET<char, MAX_NODE_NAME, double>::getNodeProb(nnode);
	}

	double setNodeSampleProb(int nnode, int *psamples, int nsamples) {
		return CATNET<char, MAX_NODE_NAME, double>::setNodeSampleProb(
							nnode, psamples, nsamples, 1/* do normalize */);
	}

	int complexity() {
		int n;
		int numnodes = numNodes();
		const int* pNumCats = numCategories();
		int complx = 0;
		for(n = 0; n < numnodes; n++)
			complx += (1 + pNumCats[n]);
		return(complx + CATNET<char, MAX_NODE_NAME, double>::complexity());
	}

	double setLoglik(double lik) {
		return CATNET<char, MAX_NODE_NAME, double>::setLoglik(lik);
	}

	double getLoglik() {
		return CATNET<char, MAX_NODE_NAME, double>::getLoglik();
	}

	double loglik();

	double setNodeLoglik(int nnode, double lik);
	double getNodeLoglik(int nnode);
	double findNodeLoglik(int nnode, double *psamples, int nsamples);

	void setNodesOrder(const int *porder) {
		return CATNET<char, MAX_NODE_NAME, double>::setNodesOrder(porder);
	}

	int *getOrder() {
		return CATNET<char, MAX_NODE_NAME, double>::getOrder();
	}

	int* catnetSample(int nsamples);

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

	int sample(double *psamples, int nsamples);
	int predict(double *psamples, int nsamples);
};

#endif /* MIXNET_CLASS_H */
