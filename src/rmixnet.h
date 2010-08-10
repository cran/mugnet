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
 * rmixnet.h
 *
 *  Created on: Sep 21, 2009
 *      Author: Nikolay Balov
 */

#ifndef RMIXNET_H_
#define RMIXNET_H_

#include "mixnet.h"
#include "poisson.h"

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

class RMixNet : public CMixNet {
public:
	RMixNet(SEXP cnet);

	SEXP genRcatnet(const char * objectName);
	SEXP genRmixnet(const char * objectName);
	SEXP genProbList(int node, int paridx, int *pcats);

	RMixNet(int nnodes, int maxpars, int maxcats = 2, const char **nodes = 0,
			const int * pnumpars = 0, const int **ppars = 0, const int *pcats = 0) : 
		CMixNet(nnodes, maxpars, maxcats, nodes, pnumpars, ppars, pcats) {
	}

};

#endif /* RMIXNET_H_ */
