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
 * utils.c
 *
 *  Created on: Nov 16, 2009
 *      Author: Nikolay Balov
 */

#include "utils.h"

extern size_t g_memcounter;
//extern FILE *g_hf;

unsigned int g_ask_mem_alloc = MAX_MEM_ALLOC;

void * CATNET_MALLOC(size_t nsize) { 
	void *pMem = 0;
	if(nsize <= 0)
		return 0;
	if(nsize > g_ask_mem_alloc) {
		Rprintf("A large chunk of memory is about to be allocated, %u. Continue? (y/n) ", (unsigned)nsize);
		char nch = getchar();
		//Rprintf("%c,%d, %c\n", nch, (int)nch, 'y');
		if(nch != 'y')
			error("Process stopped by user");
		g_ask_mem_alloc = nsize;
		// read the 0xa character
		nch = getchar();
	}
	pMem = malloc(nsize);
	if(!pMem)
		error("Insufficient memory");
	return pMem;
	g_memcounter += nsize;
	pMem = malloc(sizeof(int) + nsize);
	if(!pMem) {
		error("Insufficient memory");
		return 0;
	}
	*(int*)pMem = nsize;
	pMem = (void*)((int*)pMem + 1);
	//char str[128];
	//sprintf(str, "+%d    %d        %p\n", (int)nsize, (int)g_memcounter, pMem);
	//fprintf(g_hf,str);
	return pMem;
}

void CATNET_FREE(void *pMem) {
	if(!pMem)	
		return;
	free(pMem);
	return;
	pMem = (void*)((int*)pMem-1);
	size_t nsize = *((int*)pMem);
	g_memcounter -= nsize;	
	//char str[128];
	//sprintf(str, "-%d    %d        %p\n", (int)nsize, (int)g_memcounter, (char*)pMem+sizeof(int));
	//Rprintf(str);
	//fprintf(g_hf,str);
	free(pMem);
}

void CATNET_MEM_ERR() {
	// generate R-errorx
	error("Insufficient memory");
}

