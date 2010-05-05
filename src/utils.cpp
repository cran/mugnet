/*
 * utils.c
 *
 *  Created on: Nov 16, 2009
 *      Author: Nikolay Balov
 */

#include "utils.h"

extern size_t g_memcounter;
//extern FILE *g_hf;

int g_ask_mem_alloc = 1;

void * CATNET_MALLOC(size_t nsize) { 
	if(nsize <= 0)
		return 0;
	if(g_ask_mem_alloc && nsize > (1<<30)) {
		printf("A large chunk of memory is about to be allocated, %u. Continue? (y/n) ", (unsigned)nsize);
		char nch = getchar();
		//printf("%c,%d, %c\n", nch, (int)nch, 'y');
		if(nch != 'y')
			error("Process stopped by user");
		g_ask_mem_alloc = 0;
		// read the 0xa character
		nch = getchar();
	}
	return malloc(nsize);
	g_memcounter += nsize;
	void *pMem = malloc(sizeof(int) + nsize);
	if(!pMem) {
		error("Insufficient memory");
		return 0;
	}
	*(int*)pMem = nsize;
	pMem = (void*)((int*)pMem + 1);
	//char str[128];
	//sprintf(str, "+%d    %d        %p\n", (int)nsize, (int)g_memcounter, pMem);
	//fprintf(g_hf,str);
	//printf(str);
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
	//printf(str);
	//fprintf(g_hf,str);
	free(pMem);
}

void CATNET_MEM_ERR() {
	// generate R-errorx
	error("Insufficient memory");
}
