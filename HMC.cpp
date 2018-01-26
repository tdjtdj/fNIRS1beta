/*
 *  HMC.cpp
 *  LGCP_2D_SIM_STUDY
 *
 *  Created by Timothy Johnson on 3/11/14.
 *  Copyright 2014 University of Michigan. All rights reserved.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "randgen.h"
#include "cholesky.h"
#include "/usr/local/include/fftw3.h"
#include "fNIRS.h"

int MAXSTEPS = 50;


double likelihood(REP *rep,int P)
{
	int i;
    double loglike;
    void calculate_residuals(REP *rep,int P);
    
    loglike = 0;
    
    calculate_residuals(rep,P);

    for (i=0;i<P;i++)
        loglike += rep->precYstart*(rep->Y[i] - rep->Veta[i])*(rep->Y[i] - rep->Veta[i]);

    for (i=P;i<rep->dim_X[0];i++)
        loglike += rep->precY*rep->residuals[i]*rep->residuals[i];

    loglike *= 0.5;
    return loglike;
}

void compute_gradients(REP *rep,int P)
{
    int i,j;
    void calculate_residuals(REP *rep,int P);

    void Cx(fftw_plan corr_mat_fwd,fftw_plan corr_mat_inv,fftw_complex *ZYX,double *lambda,int mtot);
    
/*    for (i=0;i<rep->dist_length;i++) {
        rep->ZYX[i][0] = 0;
        rep->ZYX[i][1] = 0;
    }
    for (i=0;i<rep->dim_X[0];i++)
        rep->ZYX[i][0] = rep->gamma[i];
    
    Cx(rep->corr_mat_fwd,rep->corr_mat_inv,rep->ZYX,rep->eigen_val,rep->dist_length);*/
 
    calculate_residuals(rep,P);
 
    for (i=0;i<rep->dim_X[0]-P;i++)
        rep->ZYX[i][0] = rep->residuals[i+P];
 //   rep->ZYX[i][0] = rep->residuals2[i] - rep->Z[i];
 //   rep->ZYX[i][0] = rep->residuals2[i] - rep->phi*rep->ZYX[i][0];
    for (i=rep->dim_X[0]-P;i<rep->dist_length;i++)
        rep->ZYX[i][0] = 0;
    for (i=0;i<rep->dist_length;i++)
        rep->ZYX[i][1] = 0;

    Cx(rep->corr_mat_fwd,rep->corr_mat_inv,rep->ZYX,rep->eigen_val,rep->dist_length);
    
    for (i=0;i<rep->dim_X[0]-P;i++)
        rep->gradient[i] = rep->gamma[i] - rep->precY*rep->phi*rep->ZYX[i][0];
    for (i=rep->dim_X[0]-P;i<rep->dist_length;i++)
        rep->gradient[i] = rep->gamma[i];
    
}

void update_momentum_vector(REP *rep,double step_size)
{
	int i;
		
	for (i=0;i<rep->dist_length;i++)
		rep->momentum[i] -= step_size*rep->gradient[i];

}

void update_position_vectors(REP *rep,double *M,const int P,double step_size)
{
    int i;
    void calc_Z(REP *rep,const int P);
    void calW(double **W,double *Y,double *Z,double *Xb,double *Ve,const int nrow,const int ncol);
    void calAx(double *Ax,double **A,double *x,const int nrow,const int ncol);
    
    for (i=0;i<rep->dist_length;i++)
        rep->gamma[i] += step_size*rep->momentum[i]/M[i];
    
    /* calculate Z */
    calc_Z(rep,P);
    
    /* calculate W */
    calW(rep->W,rep->Y,rep->Z,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
 
    /* calculate Wdelta */
    calAx(rep->Wdelta,rep->W,rep->delta,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);

}

void compute_Hamiltonian(double *Hamiltonian,POP *pop,REP *rep,double *M)
{
	int i;
	double PE,KE;
	double prior,likelih;
	double quadform_momentum;

    prior = 0;
    quadform_momentum = 0;
    
    for (i=0;i<rep->dist_length;i++) {
        prior += rep->gamma[i]*rep->gamma[i];
        quadform_momentum += rep->momentum[i]*rep->momentum[i]/M[i];
    }
    quadform_momentum *= 0.5;
    prior *= 0.5;

    likelih = likelihood(rep,pop->P);
    
    PE = likelih + prior;
    
    KE = quadform_momentum;

	*Hamiltonian = KE + PE;
}

void HMC_update(POP *pop,SUB *sub,REP *rep,int iter,unsigned long *seed)
{
    int i,j,k,l,STEPS;
	double newlike,eps_div_2,*M_gamma;
	double Ham_new,Ham_old;
	double *save_gamma,*save_residuals,*save_Z;
    
    save_gamma = (double *)calloc(rep->dist_length,sizeof(double));
    save_Z = (double *)calloc(rep->dim_X[0]-pop->P,sizeof(double));
    save_residuals = (double *)calloc(rep->dim_X[0],sizeof(double));
 
    
    for (i=0;i<rep->dist_length;i++)
        save_gamma[i] = rep->gamma[i];
    for (i=0;i<rep->dim_X[0]-pop->P;i++)
        save_Z[i] = rep->Z[i];
    for (i=0;i<rep->dim_X[0];i++)
        save_residuals[i] = rep->residuals[i];
    
    STEPS = rpois((double)rep->maxsteps,seed);//runiform_n(40,seed) + MAXSTEPS - 20;
    eps_div_2 = rep->prop_sd[2]/2.0;
    
    M_gamma = (double *)calloc(rep->dist_length,sizeof(double));
    
    for (i=0;i<rep->dist_length;i++) {
        double tmp = exp(-rep->rho*pow(rep->dist[i],pop->expo));
        rep->ZYX[i][0] = tmp;
        rep->ZYX[i][1] = 0;
    }
    for (i=0;i<rep->dim_X[0]-pop->P;i++)
        M_gamma[i] = 2000;//rep->ZYX[0][0]*rep->precY*rep->phi*rep->phi + 1000;
    for (i=rep->dim_X[0]-pop->P;i<rep->dist_length;i++)
        M_gamma[i] = 1;

//    if (!(iter%100))
//        printf("%lf %lf %lf %lf\n",M_gamma[0],M_gamma[1],M_gamma[2],M_gamma[rep->dist_length-1]);
    
    for (i=0;i<rep->dist_length;i++)
        rep->momentum[i] = rnorm(0,sqrt(M_gamma[i]),seed);
    
//    calculate_residuals(rep,pop->P);
    compute_Hamiltonian(&Ham_old,pop,rep,M_gamma);
    
    compute_gradients(rep,pop->P);
    
    // Start Leap Frog with 1/2 steps for momentum variables
    update_momentum_vector(rep,eps_div_2);

    for (j=0;j<STEPS;j++) {
        
        // update position vectors by taking full steps
        
        update_position_vectors(rep,M_gamma,(const int)pop->P,rep->prop_sd[2]);
      
        // compute new gradients
        
        compute_gradients(rep,pop->P);

        // update momentum vectors by taking full steps
        
        update_momentum_vector(rep,rep->prop_sd[2]);
    }
    
    // subtract off half of the full step on momentum vectors to complete leap frog
    
    update_momentum_vector(rep,-eps_div_2);
    
//    calculate_residuals(rep,pop->P);
    compute_Hamiltonian(&Ham_new,pop,rep,M_gamma);
//printf("Ham_diff = %lf\n",Ham_old - Ham_new);
    if (log(kiss(seed)) < (Ham_old - Ham_new)) {
       (rep->accept[2])++;
    }
    else {
        for (i=0;i<rep->dist_length;i++)
            rep->gamma[i] = save_gamma[i];
        for (i=0;i<rep->dim_X[0]-pop->P;i++)
            rep->Z[i] = save_Z[i];
        for (i=0;i<rep->dim_X[0];i++)
            rep->residuals[i] = save_residuals[i];
    }
    (rep->attempt[2])++;
    
    free(save_gamma);
    free(save_Z);
    free(save_residuals);
    free(M_gamma);
}

