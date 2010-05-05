/*
 * rmixnet.h
 *
 *  Created on: Sep 21, 2009
 *      Author: Nikolay Balov
 */

#ifndef RMIXNET_H_
#define RMIXNET_H_

#include "mixnet.h"

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

class RMixNet : public CMixNet {
public:
	RMixNet() {
	}

	RMixNet(SEXP cnet);

	RMixNet(int nnodes, int maxpars, int maxcats = 2, const char **nodes = 0,
			const int * pnumpars = 0, const int **ppars = 0, const int *pcats = 0) : 
		CMixNet(nnodes, maxpars, maxcats, nodes, pnumpars, ppars, pcats) {
	}

	SEXP genRcatnet(const char * objectName);
	SEXP genRmixnet(const char * objectName);
	SEXP genProbList(int node, int paridx, int *pcats);
};

#endif /* RMIXNET_H_ */
