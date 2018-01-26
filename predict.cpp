/*
 *  predict.cpp
 *  LGCPmarkGPU
 *
 *  Created by Timothy Johnson on 10/30/10.
 *  Copyright 2010 University of Michigan. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LGCP.h"
#include "randgen.h"

//#include <mach/mach_time.h>

extern int *midx;

void subQ(float *ans,SUB *sub,GP *gp,float cell_vol,int ntot,int nsubtypes,double *expyl)
{
	int i,j,isubtype;
	double *yprop,y;
	float likelihood_tst2(SUB *sub,GP *gp,double expyl,float cell_vol,int ntot);
	
	yprop = (double *)calloc(nsubtypes,sizeof(double));
//	expyl = (double *)calloc(nsubtypes,sizeof(double));
		

/*	for (i=0;i<ntot;i++) {
//		if (*(mask+i)) {
			j = *(midx+i);
			for (isubtype=0;isubtype<nsubtypes;isubtype++) {
				y = (double)(*(sub->mu) + gp[isubtype+1].Y.realp[j]);
				yprop[isubtype] += y*sub->ncell[j];
				expyl[isubtype] += exp(y); 	
			}
//		}
	}*/

	for (isubtype=0;isubtype<nsubtypes;isubtype++) 
		ans[isubtype] = likelihood_tst2(sub,&gp[isubtype+1],expyl[isubtype],cell_vol,ntot);

//		ans[isubtype] = (float)(yprop[isubtype] - expyl[isubtype]*cell_vol);

	free(yprop);
}

void predict_type(SUB *sub,GP *gp,float cell_vol,int ntot,int nsub,int nsubtypes)
{
	int isub,isubtype;
	float *ans;
	double *expyl;
	double like_pre(GP *gp,int ntot);

	expyl = (double *)calloc(nsubtypes,sizeof(double));
	for (isubtype=0;isubtype<nsubtypes;isubtype++)
		expyl[isubtype] = like_pre(&gp[isubtype+1],ntot);

	ans = (float *)calloc(nsubtypes,sizeof(float));
	
	for (isub=0;isub<nsub;isub++) {
		subQ(ans,&(sub[isub]),gp,cell_vol,ntot,nsubtypes,expyl);
		for (isubtype=0;isubtype<nsubtypes;isubtype++) {
			sub[isub].predict[isubtype] += expf(ans[isubtype] - ans[sub[isub].subtype]);
		}
	}
	free(ans);
	free(expyl);	 
}
