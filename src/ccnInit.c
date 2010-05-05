
#include "mixnet_rexport.h"

size_t g_memcounter;
int g_netcounter;
//FILE *g_hf;

static const R_CallMethodDef R_CallDef[] = {
        {"mgSearchOrderC", (DL_FUNC)&searchOrder, 12},
	{"mgPredictC", (DL_FUNC)&predict, 2},
	{"mgSampleC", (DL_FUNC)&sample, 2},
	{"mgNodeLoglikC", (DL_FUNC)&nodeLoglik, 4},
	{NULL, NULL, 0},
};

void R_init_mugnet(DllInfo *info)
{
	R_registerRoutines(info,NULL,R_CallDef,NULL,NULL);
	g_memcounter = 0;
	//g_hf = fopen("dump.txt", "wt");
}

void R_unload_mugnet(DllInfo *info)
{
	//fclose(g_hf);
}
