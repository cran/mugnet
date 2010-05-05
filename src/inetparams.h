/*
 * netparams.h
 *
 *  Created on: Nov 16, 2009
 *      Author: Nikolay Balov
 */

#ifndef I_NETPARAMS_H_
#define I_NETPARAMS_H_

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
	virtual t_prob findNodeLoglik(int nnode, double *psamples, int nsamples) = 0;
	virtual void setNodesOrder(const int *porder) = 0;
	virtual int * getOrder() = 0;
	virtual int* catnetSample(int nsamples) = 0;

	virtual void setNodeBetas(int nnode, double *pbetas) = 0;
	virtual void setNodeSigma(int nnode, double fSigma) = 0;
	virtual const double **setBetas(double **pbetas, int nbetas) = 0;
	virtual const double *setSigmas(double *psigmas, int nsigmas) = 0;
	virtual const double **betas() = 0;
	virtual const double *sigmas() = 0;

	virtual int pc_size() = 0;
	virtual int pcC_size() = 0;
	virtual const t_prob *pc() = 0;
	virtual const t_prob *pcC() = 0;
	virtual const t_prob *qc() = 0;
	virtual const t_prob *qcC() = 0;
	virtual void release_pqcC() = 0;

	virtual int findNodeMarginalProb(int nnode, t_prob *psamples, int nsamples) = 0;
	virtual int findNodeJointProb(int nnode, int *pnodes, int numnodes, t_prob *psamples, int nsamples) = 0;
};


#endif /* I_NETPARAMS_H_ */
