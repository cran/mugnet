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

#include "mixnet_rexport.h"

size_t g_memcounter;
int g_netcounter;
//FILE *g_hf;

static const R_CallMethodDef R_CallDef[] = {
        {"mgSearchOrderC", (DL_FUNC)&searchOrder, 14},
	{"mgPredictC", (DL_FUNC)&predict, 3},
	{"mgSampleC", (DL_FUNC)&sample, 3},
	{"mgNodeLoglikC", (DL_FUNC)&nodeLoglik, 5},
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
