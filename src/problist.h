/*
 * problist.h
 *
 *  Created on: Nov 16, 2009
 *      Author: Nikolay Balov
 */

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <float.h>

#ifndef PROBLIST_H_
#define PROBLIST_H_

#ifdef DEBUG_ON 
#include <iostream.h>
#endif

template<class t_prob>
struct PROB_LIST {
	t_prob *pProbs;
	int nProbSize;
	int numCats;
	int numPars;
	int *numParCats;
	int *pBlockSize;
	t_prob loglik;

	void reset() {
		if (numParCats)
			CATNET_FREE(numParCats);
		if (pBlockSize)
			CATNET_FREE(pBlockSize);
		if (pProbs)
			CATNET_FREE(pProbs);
		numPars = 0;
		numCats = 0;
		numParCats = 0;
		pBlockSize = 0;
		pProbs = 0;
		nProbSize = 0;
		loglik = 0;
	}


	PROB_LIST() {
		numPars = 0;
		numCats = 0;
		numParCats = 0;
		pBlockSize = 0;
		pProbs = 0;
		nProbSize = 0;
		loglik = 0;
	}

	PROB_LIST<t_prob>& operator =(const PROB_LIST<t_prob> &plist) {
		numPars = plist.numPars;
		numCats = plist.numCats;
		if (numParCats)
			CATNET_FREE(numParCats);
		numParCats = (int*) CATNET_MALLOC(numPars * sizeof(int));
		memset(numParCats, 0, numPars * sizeof(int));
		if (plist.numParCats)
			memcpy(numParCats, plist.numParCats, numPars * sizeof(int));
		if (pBlockSize)
			CATNET_FREE(pBlockSize);
		pBlockSize = (int*) CATNET_MALLOC(numPars * sizeof(int));
		memset(pBlockSize, 0, numPars * sizeof(int));
		if (plist.pBlockSize)
			memcpy(pBlockSize, plist.pBlockSize, numPars * sizeof(int));
		nProbSize = plist.nProbSize;
		if (pProbs)
			CATNET_FREE(pProbs);
		pProbs = 0;
		if(nProbSize > 0) {
			pProbs = (t_prob*) CATNET_MALLOC(nProbSize * sizeof(t_prob));
			memset(pProbs, 0, nProbSize * sizeof(t_prob));
			if (plist.pProbs) {
				// doesn't work
				//memcmp(pProbs, plist.pProbs, nProbSize * sizeof(t_prob));
				//printf("pProbs = ");
				for(int i = 0; i < nProbSize; i++) {
					pProbs[i] = plist.pProbs[i];
				//printf("%f ", pProbs[i]);
				}
				//printf("\n");
			}
		}
		loglik = plist.loglik;
		return *this;
	}

	PROB_LIST(int ncats, int nmaxcats = 2, int npars = 0, int *parcats = 0,
			t_prob *pprobs = 0, int probsize = 0) {
		int i;
		if (ncats < 1 || nmaxcats < 1 || npars < 0 || (npars && !parcats))
			return;
		numPars = npars;
		numCats = ncats;
		numParCats = 0;
		pBlockSize = 0;
		pProbs = 0;
		nProbSize = 0;
		loglik = 0;
		//cout << "ncats = " << ncats << ", nmaxcats = " << nmaxcats << ", npars = " << npars << ", 
		//probsize = " << probsize << "\n";
		if (numPars > 0) {
			numParCats = (int*) CATNET_MALLOC(numPars * sizeof(int));
			if (parcats)
				memcpy(numParCats, parcats, numPars * sizeof(int));
			else
				for (i = 0; i < numPars; i++)
					numParCats[i] = nmaxcats;

			pBlockSize = (int*) CATNET_MALLOC(numPars * sizeof(int));
			pBlockSize[numPars - 1] = ncats;
			//cout << "pBlockSize[numPars-1] = " << pBlockSize[numPars-1] << "\n";
			for (i = numPars - 1; i > 0; i--) {
				if (parcats[i] < 1 || parcats[i] > nmaxcats) {
					CATNET_FREE(pBlockSize);
					pBlockSize = 0;
					numPars = 0;
					return;
				}
				pBlockSize[i - 1] = parcats[i] * pBlockSize[i];
				//cout << i << " pBlockSize[i-1] = " << pBlockSize[i-1] << "\n";
			}
			nProbSize = pBlockSize[0] * numParCats[0];
		} else
			nProbSize = ncats;
		
		//cout << "nProbSize = " << nProbSize << "\n";
		pProbs = (t_prob*) CATNET_MALLOC(nProbSize * sizeof(t_prob));
		memset(pProbs, 0, nProbSize * sizeof(t_prob));
		if (pProbs && pprobs) {
			if (probsize != nProbSize) {
				//cout << "probsize != nProbSize\n";
				return;
			}
			memcpy(pProbs, pprobs, nProbSize * sizeof(t_prob));
		}
		else if(numPars == 0) {
			// set the default probability
			t_prob fsum = 0;
			for(i = 0; i < numCats-1; i++) {
				pProbs[i] = (t_prob)1/(t_prob)numCats;
				fsum += pProbs[i];
			}
			pProbs[numCats-1] = 1 - fsum;
			//printf("PROB: ");
			//for(i = 0; i < numCats; i++) 
			//	printf("%f, ", pProbs[i]);
			//printf("\n");
			
		}
	}

	~PROB_LIST() {
//if(ggg > 0) ggg--;
//printf("~PROB_LIST:: %p\t\t%d\n", this, ggg);
		reset();
	}

	void set_zero() {
		if (pProbs)
			memset(pProbs, 0, nProbSize * sizeof(t_prob));
		loglik = 0;
	}

	void normalize() {
		int i, k;
		t_prob pp;
		if (!pProbs)
			return;
		k = 0;
		while (k < nProbSize) {
			pp = 0;
			for (i = 0; i < numCats; i++)
				pp += pProbs[k + i];
			if (pp > 0)
				pp = 1 / pp;
			for (i = 0; i < numCats; i++) {
				pProbs[k + i] *= pp;
			}
			k += numCats;
		}
	}

	t_prob loglikelihood() {
		int i, k;
		t_prob pp;
		loglik = 0;
		k = 0;
//cout << "nProbSize = " << nProbSize << "\n";
//cout << "npars = " << numPars << ", probsize = " << nProbSize << ",loglik: ";
//for(i = 0; i < nProbSize; i++)
//	cout << pProbs[i] << ", ";
//cout << "\n";
		while (k < nProbSize) {
			pp = 0;
			for (i = 0; i < numCats; i++)
				pp += pProbs[k + i];
			if (pp <= 0) {
				/////////////////////////////////////////////////////////////////
				// The next condition determines whether parent sets with  
				// non-sample-populated slots should be abandoned
				//loglik = -FLT_MAX;
				//break;
				/////////////////////////////////////////////////////////////////
			}
			else {
				pp = 1 / pp;
//cout << pp << "\n";
				for (i = 0; i < numCats; i++) {
					if(pProbs[k + i] > 0) {
						loglik += pProbs[k + i] * log(pProbs[k + i] * pp);
					}
				}
			}
			k += numCats;
		}
//cout << "PROB_LIST::loglik = " << loglik << "\n";
		return loglik;  
	}

	t_prob *find_slot(t_prob* prob, int *pcats, int parid) {
		if (prob == 0)
			prob = pProbs;
		if (parid >= numPars || pcats == 0)
			return prob;
		int parcat = pcats[parid];
		if (parcat < 0 && parcat >= numParCats[parid])
			return 0;
		if (parid == numPars - 1)
			return prob + parcat * pBlockSize[parid];
		return find_slot(prob + parcat * pBlockSize[parid], pcats, parid + 1);
	}

};

#endif /* PROBLIST_H_ */

