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
 * utils.h
 *
 *  Created on: Nov 16, 2009
 *      Author: Nikolay Balov
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <float.h>
#include <time.h>

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#ifndef CATNET_PI
#define CATNET_PI	(double)3.14159265358979323846264338327950288
#endif
#ifndef CATNET_PI2
#define CATNET_PI2	(2*(double)CATNET_PI)
#endif

#define PI2	(2*(double)CATNET_PI)

#define MAX_MEM_ALLOC 4e9

void * CATNET_MALLOC(size_t nsize);
void CATNET_FREE(void *pMem);
void CATNET_MEM_ERR();

template<class t_elem>
void _quick_sort(t_elem *plist, int nlist){
	t_elem pivot;
	int j, nless, ngreater;
	if(nlist <= 1)
		return;
	t_elem *paux = (t_elem*)malloc(nlist*sizeof(t_elem));
	nless = 0;
	ngreater = nlist-1;
	pivot = plist[0];
	for(j = 1; j < nlist; j++) {
		if(plist[j] <= pivot) {
			paux[nless] = plist[j];
			nless++;
		}
		else {
			paux[ngreater] = plist[j];
			ngreater--;
		}
	}
	if(nless > 0)
		_quick_sort<t_elem>(paux, nless);
	if(nlist > ngreater + 1)
		_quick_sort<t_elem>(paux + ngreater + 1, nlist - ngreater - 1);
	paux[nless] = pivot;
	memcpy(plist, paux, nlist*sizeof(t_elem));
	free(paux);
	return;
}


template<class t_elem>
t_elem _gen_std_normal_var() {
	/* ISO C pseudo random generator */
	/* include stdlib.h and math.h */
	t_elem u, v;
	u = (t_elem)rand() / (t_elem)RAND_MAX;
	v = (t_elem)rand() / (t_elem)RAND_MAX;
	return(sqrt(-2*log(u)) * cos(CATNET_PI2*v));
}

#endif /* UTILS_H_ */

