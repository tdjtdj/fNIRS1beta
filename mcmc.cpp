/*
 *  mcmc.c
 *  LGCP
 *
 *  Created by Timothy D. Johnson on 11/10/15.
 *  Copyright 2015 University of Michigan. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "randgen.h"
#include "cholesky.h"
#include <time.h>
#include "fNIRS.h"

int MAX_ITER =  12000;
int BURN_IN  =  2000;
extern sDLM *dlmStruc;
extern int maxP;

void adjust_acceptance(double x,double *X,double rate)
{
	double y;
	
	y = 1. + 2.*(x-rate)*(x-rate)*(x-rate);
	if (y < .8)
		y = .8;
	if (y > 1.2)
		y = 1.2;
	*X *= y;
}

void adjust_acceptance2(double x,double *X,double rate)
{
    double y;
    
    y = 1. + 10.*(x-rate)*(x-rate)*(x-rate);
    if (y < .8)
        y = .8;
    if (y > 1.2)
        y = 1.2;
    *X *= y;
}


void mcmc(POP *pop,unsigned long *seed) {
    int i,j,k,iter,isub,irep,Z_cnt;
    double ALPHA[3] = {4.9,4.9,.2};
    double **loglik;
    int nrow_loglik;
    char *C,*S;
    FILE *fout,*fout3,*fout4,*fout5,**fout_eta,**fout_delta,**fout_sub_beta,*fout_pop_beta;
    FILE **fout_nknots,**fout_knots,**fout_res,**fout_sub_precY,*fout_reprec,**fout_delta2;
    FILE *fout_ll,**fout_dlm,*fout_veta,*fout_wd;
    
    void draw_knot_locations(REP *rep,int sdegree,int *flag,unsigned long *seed);
    int knot_birth_death(REP *rep,POP *pop,const int sdegree,int iter,unsigned long *seed);
    void itoa(int n,char s[]);
 
    void DLMtst(REP *rep,int iter,int flag,unsigned long *seed);
    void DLM(REP *rep,int P,unsigned long *seed);
    
    void draw_precYstart(REP *rep,unsigned long *seed);
    void draw_preceta(REP *rep,double **Q,unsigned long *seed);
    void draw_reprec(POP *pop,SUB *sub,unsigned long *seed);
 
    void calculate_residuals(REP *rep,int P);
    void calculate_res1(REP *rep,int P);
    void calculate_res2(REP *rep,int P);
    void calculate_res3(REP *rep,int P);
    void calculate_res4(REP *rep,int P);
    void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    
    void draw_beta_eta(POP *pop,SUB *sub,REP *rep,unsigned long *seed);
    void draw_eta(REP *rep,double **Q,int P,unsigned long *seed);
    void draw_sub_beta(POP *pop,SUB *sub,unsigned long *seed);
    void draw_pop_beta(POP *pop,SUB *sub,unsigned long *seed);
     
    void DIC(POP *pop,REP *rep,unsigned long *seed);
/*    
//    int minP = 20;
    int minP = 2;
    fout_ll = fopen("ll.dat","w");
//    nrow_loglik = 10*20*(maxP-minP+1);
    nrow_loglik = 10*15*(maxP-minP+1);
    loglik = (double **)calloc(nrow_loglik,sizeof(double *));
    for (i=0;i<nrow_loglik;i++)
        loglik[i] = (double *)calloc(4,sizeof(double));
//    double DELTA[30] = {0.9990, 0.9991, 0.9992, 0.9993, 0.9994, 0.9995, 0.9996, 0.9997, 0.9998, 0.9999};
    double DELTA[10] = {0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99};
//    double BETA[32] = {0.9800, 0.981, 0.982, 0.983, 0.984, 0.985, 0.986, 0.987, 0.988, 0.989, 0.990, 0.991,
//0.992, 0.993, 0.994, 0.995, 0.996, 0.997, 0.998, 0.999};
    double BETA[15] = {0.7, 0.72, 0.74, 0.76, 0.78, 0.80, 0.82, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98};
 //   double beta,delta = 0.998 - 0.0001;
    int cnt = 0;
    for (int id=0;id<10;id++) {
 //        for (int ib=0;ib<20;ib++) {
         for (int ib=0;ib<15;ib++) {
//           for (int P=minP;P<=maxP;P++) {
            for (int P=minP;P<=maxP;P++) {
                loglik[cnt][0] = DELTA[id];
                loglik[cnt][1] = BETA[ib];
                loglik[cnt][2] = (double)P;
                cnt++;
             }
        }
    }
   
    printf("cnt = %d nrow_loglik = %d\n",cnt,nrow_loglik);*/
    int reps = 0;
    for (isub=0;isub<pop->N_SUBS;isub++)
            reps += pop->sub[isub].N_REPS;
    
    fout_eta = (FILE **)calloc(pop->N_SUBS*reps,sizeof(FILE *));
    fout_res = (FILE **)calloc(pop->N_SUBS*reps,sizeof(FILE *));
    fout_delta = (FILE **)calloc(pop->N_SUBS*reps,sizeof(FILE *));
    fout_delta2 = (FILE **)calloc(pop->N_SUBS*reps,sizeof(FILE *));
    fout_sub_precY = (FILE **)calloc(pop->N_SUBS*reps,sizeof(FILE *));
    fout_nknots = (FILE **)calloc(pop->N_SUBS*reps,sizeof(FILE *));
    fout_knots = (FILE **)calloc(pop->N_SUBS*reps,sizeof(FILE *));
    fout_sub_beta = (FILE **)calloc(pop->N_SUBS,sizeof(FILE *));
    fout_dlm = (FILE **)calloc(pop->N_SUBS*reps,sizeof(FILE *));;
    fout_veta = fopen("veta_draws.dat","w");
    fout_wd = fopen("wdelta_draws.dat","w");
    C = (char *)malloc(3*sizeof(char));
    S = (char *)malloc(100*sizeof(char));
    
    fout_reprec = fopen("pop_reprec.dat","w");
    fout_pop_beta = fopen("pop_betadraws.dat","w");
    for (isub=0;isub<pop->N_SUBS;isub++) {
        for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
            S = strcpy(S,"DLM_draws");
            itoa(isub,C);
            S = strcat(S,C);
            itoa(irep,C);
            S = strcat(S,C);
            S = strcat(S,".dat");
            //printf("S: %s\n",S);
            fout_dlm[isub*pop->N_SUBS+irep] = fopen(S,"w");
            S = strcpy(S,"etadraws");
            itoa(isub,C);
            S = strcat(S,C);
            itoa(irep,C);
            S = strcat(S,C);
            S = strcat(S,".dat");
            //printf("S: %s\n",S);
            fout_eta[isub*pop->N_SUBS+irep] = fopen(S,"w");
            S = strcpy(S,"resdraws");
            itoa(isub,C);
            S = strcat(S,C);
            itoa(irep,C);
            S = strcat(S,C);
            S = strcat(S,".dat");
            //printf("S: %s\n",S);
            fout_res[isub*pop->N_SUBS+irep] = fopen(S,"w");
            S = strcpy(S,"deltadraws");
            itoa(isub,C);
            S = strcat(S,C);
            itoa(irep,C);
            S = strcat(S,C);
            S = strcat(S,".dat");
            //printf("S: %s\n",S);
            fout_delta[isub*pop->N_SUBS+irep] = fopen(S,"w");
            S = strcpy(S,"delta2draws");
            itoa(isub,C);
            S = strcat(S,C);
            itoa(irep,C);
            S = strcat(S,C);
            S = strcat(S,".dat");
            //printf("S: %s\n",S);
            fout_delta2[isub*pop->N_SUBS+irep] = fopen(S,"w");
            S = strcpy(S,"nknotdraws");
            itoa(isub,C);
            S = strcat(S,C);
            itoa(irep,C);
            S = strcat(S,C);
            S = strcat(S,".dat");
            //printf("S: %s\n",S);
            fout_nknots[isub*pop->N_SUBS+irep] = fopen(S,"w");
            S = strcpy(S,"knotdraws");
            itoa(isub,C);
            S = strcat(S,C);
            itoa(irep,C);
            S = strcat(S,C);
            S = strcat(S,".dat");
            //printf("S: %s\n",S);
            fout_knots[isub*pop->N_SUBS+irep] = fopen(S,"w");
             S = strcpy(S,"sub_precYdraws");
            itoa(isub,C);
            S = strcat(S,C);
            itoa(irep,C);
            S = strcat(S,C);
            S = strcat(S,".dat");
            //printf("S: %s\n",S);
            fout_sub_precY[isub*pop->N_SUBS+irep] = fopen(S,"w");
        }
        S = strcpy(S,"sub_betadraws");
        itoa(isub,C);
        S = strcat(S,C);
        S = strcat(S,".dat");
        //printf("S: %s\n",S);
        fout_sub_beta[isub] = fopen(S,"w");
    }
    free(C);
    free(S);
    int sdegree = 4;
    int aaa = 0;
    FILE *foutfit = fopen("fit.dat","w");
    for (iter=0;iter<=MAX_ITER;iter++) {

        for (isub=0;isub<pop->N_SUBS;isub++) {
 
            for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
//printf("A\n");fflush(stdout);
                if (iter == 0)
//                    draw_eta(&(pop->sub[isub].rep[irep]),pop->Q2,pop->P,seed);
                    draw_beta_eta(pop,&(pop->sub[isub]),&(pop->sub[isub].rep[irep]),seed);
//printf("B %d\n",pop->sub[isub].rep[irep].nKnots);fflush(stdout);
                aaa = knot_birth_death(&(pop->sub[isub].rep[irep]),pop,sdegree,iter,seed);
//printf("C\n");fflush(stdout);
                for (i=0;i<pop->sub[isub].rep[irep].nKnots-sdegree*2;i++) {
                    int flag = 0;
                    draw_knot_locations(&(pop->sub[isub].rep[irep]),sdegree,&flag,seed);
                    if (flag) {
//                        draw_beta_eta(pop,&(pop->sub[isub]),&(pop->sub[isub].rep[irep]),pop->P,seed); 
                        if ((pop->sub[isub].rep[irep].dim_W[1] > 0) && (iter > 1000000000)) {
                            DLMtst(&(pop->sub[isub].rep[irep]),iter,0,seed);
                            DLMtst(&(pop->sub[isub].rep[irep]),iter,1,seed);
                            DLMtst(&(pop->sub[isub].rep[irep]),iter,2,seed);
                            DLM(&(pop->sub[isub].rep[irep]),pop->sub[isub].rep[irep].P,seed);
                        }
//                        else {
//                            DLMtst(&(pop->sub[isub].rep[irep]),&delta,&beta,&(pop->P),iter,1,seed);                        
//                        }
//                        draw_preceta(&(pop->sub[isub].rep[irep]),pop->Q2,seed);
//                        draw_precYstart(&(pop->sub[isub].rep[irep]),pop->P,seed);
                    }
                }
                draw_beta_eta(pop,&(pop->sub[isub]),&(pop->sub[isub].rep[irep]),seed);
                if ((pop->sub[isub].rep[irep].dim_W[1] > 0) && (iter > -1)) {
                    DLMtst(&(pop->sub[isub].rep[irep]),iter,0,seed);
                    DLMtst(&(pop->sub[isub].rep[irep]),iter,1,seed);
                    DLMtst(&(pop->sub[isub].rep[irep]),iter,2,seed);
                    DLM(&(pop->sub[isub].rep[irep]),pop->sub[isub].rep[irep].P,seed);
                }
//                else {
//                    DLMtst(&(pop->sub[isub].rep[irep]),&delta,&beta,&(pop->P),iter,1,seed);
//                    DLM(&(pop->sub[isub].rep[irep]),pop->P,delta,beta,seed);
//                }
                draw_preceta(&(pop->sub[isub].rep[irep]),pop->Q2,seed);
                draw_precYstart(&(pop->sub[isub].rep[irep]),seed);
        
                
                
                fprintf(fout_nknots[isub*pop->N_SUBS+irep],"%d ",pop->sub[isub].rep[irep].nKnots-8);
                if (!(iter%100)) {
                for (i=0;i<pop->sub[isub].rep[irep].dim_V[0];i++)
                    fprintf(fout_veta,"%lf ",pop->sub[isub].rep[irep].Veta[i]);
                fprintf(fout_veta,"\n");
                 for (i=0;i<pop->sub[isub].rep[irep].dim_W[0];i++)
                    fprintf(fout_wd,"%lf ",pop->sub[isub].rep[irep].Wdelta[i]);
                fprintf(fout_wd,"\n");
                }
                if ((iter>BURN_IN)) {
                    for (i=pop->sub[isub].rep[irep].P;i<pop->sub[isub].rep[irep].dim_W[0];i++) {
                        for (j=0;j<pop->sub[isub].rep[irep].dim_W[1];j++) {
                            pop->sub[isub].rep[irep].mdelta[i*pop->sub[isub].rep[irep].dim_W[1]+j] += pop->sub[isub].rep[irep].delta[i*pop->sub[isub].rep[irep].dim_W[1]+j];
                            pop->sub[isub].rep[irep].mdelta2[i*pop->sub[isub].rep[irep].dim_W[1]+j] += pop->sub[isub].rep[irep].delta[i*pop->sub[isub].rep[irep].dim_W[1]+j]*pop->sub[isub].rep[irep].delta[i*pop->sub[isub].rep[irep].dim_W[1]+j];
                        }
                    }
 
                     for (i=0;i<pop->sub[isub].rep[irep].dim_X[0];i++) {
                            pop->sub[isub].rep[irep].md_Y[i] += 1./sqrt(pop->sub[isub].rep[irep].d_Y[i]);
                            pop->sub[isub].rep[irep].sd_Y[i] += 1./pop->sub[isub].rep[irep].d_Y[i];
                     }
                    
                   
 //                   for (i=pop->P;i<pop->sub[isub].rep[irep].dim_X[0];i++)
 //                       fprintf(fout_sub_precY[isub*pop->N_SUBS+irep],"%lf ",pop->sub[isub].rep[irep].d_Y[i]);
//                    fprintf(fout_sub_precY[isub*pop->N_SUBS+irep],"\n");fflush(NULL);
                   for (i=0;i<pop->sub[isub].rep[irep].dim_V[1];i++)
                        fprintf(fout_eta[isub*pop->N_SUBS+irep],"%lf ",pop->sub[isub].rep[irep].eta[i]);
                    fprintf(fout_eta[isub*pop->N_SUBS+irep],"\n");
                   // if (iter == MAX_ITER) {
                        for (i=0;i<pop->sub[isub].rep[irep].dim_V[0];i++)
                            pop->sub[isub].rep[irep].mean_Y[i] += pop->sub[0].rep[irep].Y[i];
                   // }
                    
                    double TTT=0;
                    for (i=pop->sub[isub].rep[irep].P;i<pop->sub[isub].rep[irep].dim_X[0];i++)
                        TTT += pop->sub[isub].rep[irep].d_Y[i];

                    DIC(pop,&(pop->sub[isub].rep[irep]),seed);
                    for (i=0;i<pop->sub[isub].rep[irep].dim_X[0];i++) {
                        pop->sub[isub].rep[irep].mean_res[i] += pop->sub[isub].rep[irep].residuals[i];
                        pop->sub[isub].rep[irep].mean_d_Y[i] += pop->sub[isub].rep[irep].d_Y[i];
                        pop->sub[isub].rep[irep].mean_fit[i] += pop->sub[isub].rep[irep].Y[i] - pop->sub[isub].rep[irep].residuals[i];
                    }
                                       
                     
                    for (i=0;i<pop->sub[isub].rep[irep].dim_V[0];i++)
                        pop->sub[isub].rep[irep].mVeta[i] += pop->sub[isub].rep[irep].Veta[i];
                    for (i=0;i<pop->sub[isub].rep[irep].dim_W[0];i++)
                        pop->sub[isub].rep[irep].mWdelta[i] += pop->sub[isub].rep[irep].Wdelta[i];
                    
                    fprintf(fout_dlm[isub*pop->N_SUBS+irep],"%lf %lf %d %d\n",pop->sub[isub].rep[irep].df_delta1,pop->sub[isub].rep[irep].df_delta2,pop->sub[isub].rep[irep].P,aaa);
                    
                    for (i=4;i<pop->sub[isub].rep[irep].nKnots-4;i++)
                        fprintf(fout_knots[isub*pop->N_SUBS+irep],"%lf ",pop->sub[isub].rep[irep].knots[i]);
                        
                 }
               
                if (!(iter%100)) {
                    printf("iter = %6d\t Sub %d, Rep %d \t int knots = %d",iter,isub,irep,pop->sub[isub].rep[irep].nKnots-2*4);
                    printf("\t precY0 = %10.6lf \t preceta = %10.6lf\n",pop->sub[isub].rep[irep].d_Y[pop->sub[isub].rep[irep].dim_X[0]-1],pop->sub[isub].rep[irep].preceta);
                    printf("\t df_delta1 = %10.6lf \t df_delta2 = %10.6lf \t P = %3d\n",pop->sub[isub].rep[irep].df_delta1,pop->sub[isub].rep[irep].df_delta2,pop->sub[isub].rep[irep].P);
                  //  printf("\t\t precdelta = %10.6lf, preceta = %10.6lf\n\n",pop->sub[isub].rep[irep].precdelta,pop->sub[isub].rep[irep].preceta);
                    
                 }
                 if (!(iter%100) && (iter <= BURN_IN) && (iter > 1)) {
                    double rt;
//                    if (iter > 2000) {
                    rt = (double)pop->sub[isub].rep[irep].accept[0]/(double)pop->sub[isub].rep[irep].attempt[0];
                    adjust_acceptance2(rt,&(pop->sub[isub].rep[irep].prop_sd[0]),0.35);
                    pop->sub[isub].rep[irep].accept[0] = pop->sub[isub].rep[irep].attempt[0] = 0;
                    printf("\t\t rt = %lf proposal = %lf\n",rt,pop->sub[isub].rep[irep].prop_sd[0]);
//                    }
                    rt = (double)pop->sub[isub].rep[irep].accept[1]/(double)pop->sub[isub].rep[irep].attempt[1];
                    adjust_acceptance2(rt,&(pop->sub[isub].rep[irep].prop_sd[1]),0.35);
                    pop->sub[isub].rep[irep].accept[1] = pop->sub[isub].rep[irep].attempt[1] = 0;
                    printf("\t\t rt = %lf proposal = %lf\n",rt,pop->sub[isub].rep[irep].prop_sd[1]);
   
                    rt = (double)pop->sub[isub].rep[irep].accept[3]/(double)pop->sub[isub].rep[irep].attempt[3];
                    adjust_acceptance2(rt,&(pop->sub[isub].rep[irep].prop_sd[3]),0.35);
                    pop->sub[isub].rep[irep].accept[3] = pop->sub[isub].rep[irep].attempt[3] = 0;
                    printf("\t\t rt = %lf proposal = %lf\n",rt,pop->sub[isub].rep[irep].prop_sd[3]);
                 }
                  fflush(NULL);
            }
//            draw_sub_beta(pop,&(pop->sub[isub]),seed);
         
            if ((iter>BURN_IN)) {
                for (i=0;i<pop->Nb*(pop->Ns);i++)
                    fprintf(fout_sub_beta[isub],"%lf ",pop->sub[isub].beta[i]);
                fprintf(fout_sub_beta[isub],"\n");
/*                for (i=0;i<pop->sub[isub].N_REPS;i++)
                    for (j=0;j<3;j++)
                        fprintf(fout_sub_precY[isub],"%lf ",pop->sub[isub].rep[i].precY[j]);
                fprintf(fout_sub_precY[isub],"\n");*/
            }
        }
 //       printf("here1\n");fflush(NULL);
        if (pop->non_parm || pop->GRP)
            draw_reprec(pop,pop->sub,seed);
 //       printf("here2\n");fflush(NULL);
        if (pop->GRP) {
            draw_pop_beta(pop,pop->sub,seed);
        }
        if (!(iter%100)) {
            double nmean[25];
            for (k=0;k<pop->Ns;k++) {
                printf("Stimulus %d\n",k);     
                if (pop->GRP) {
                    for (i=0;i<pop->Nb;i++)
                        printf("\t %g ",pop->re_prec[i+k*pop->Nb]);
                    printf("\n\n");
                }
                for (i=0;i<pop->N_SUBS;i++) {
                    for (j=0;j<pop->Nb;j++) {
                        printf("\t %10.6lf \t ",pop->sub[i].beta[j+pop->Nb*k]);
                    }
                    printf("\n");
                }
                printf("\n");
                if (pop->GRP) {
                    for (i=0;i<pop->Nb;i++) 
                        printf("\t %10.6lf \t ",pop->beta[i + pop->Nb*k]);
                    printf("\n\n");fflush(stdout);
                }
            }
        }
        if ((iter>BURN_IN)) {
            for (i=0;i<pop->Nb*pop->Ns;i++)
                fprintf(fout_pop_beta,"%lf ",pop->beta[i]);
            fprintf(fout_pop_beta,"\n");
            for (i=0;i<pop->Ns;i++)
                for (j=0;j<pop->Nb;j++)
                fprintf(fout_reprec,"%lf ",pop->re_prec[j+i*pop->Nb]);
            fprintf(fout_reprec,"\n");fflush(NULL);
            
        }
    }
    fclose(fout_veta);
    fclose(fout_wd);
    fclose(fout_pop_beta);
    fclose(fout_reprec);
    fclose(foutfit);
    for (isub=0;isub<pop->N_SUBS;isub++) {
        for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
            fclose(fout_dlm[isub*pop->N_SUBS+irep]);
            fclose(fout_eta[isub*pop->N_SUBS+irep]);
            fclose(fout_nknots[isub*pop->N_SUBS+irep]);
            fclose(fout_knots[isub*pop->N_SUBS+irep]);
        }
        fclose(fout_sub_beta[isub]);
    }


   for (isub=0;isub<pop->N_SUBS;isub++) {
        for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
            for (i=pop->sub[isub].rep[irep].P;i<pop->sub[isub].rep[irep].dim_W[0];i++) {
                for (j=0;j<pop->sub[isub].rep[irep].dim_W[1];j++) {
                    pop->sub[isub].rep[irep].mdelta[i*pop->sub[isub].rep[irep].dim_W[1]+j] /= (double)(MAX_ITER-BURN_IN);
                    fprintf(fout_delta[isub*pop->N_SUBS+irep],"%lf ",pop->sub[isub].rep[irep].mdelta[i*pop->sub[isub].rep[irep].dim_W[1]+j]);
                    fprintf(fout_delta2[isub*pop->N_SUBS+irep],"%lf ",sqrt(pop->sub[isub].rep[irep].mdelta2[i*pop->sub[isub].rep[irep].dim_W[1]+j]/(MAX_ITER-BURN_IN)- pop->sub[isub].rep[irep].mdelta[i*pop->sub[isub].rep[irep].dim_W[1]+j]*pop->sub[isub].rep[irep].mdelta[i*pop->sub[isub].rep[irep].dim_W[1]+j]));
                }
                fprintf(fout_delta[isub*pop->N_SUBS+irep],"\n");                  
                fprintf(fout_delta2[isub*pop->N_SUBS+irep],"\n");  
            }                
            fclose(fout_delta[isub*pop->N_SUBS+irep]);
            fclose(fout_delta2[isub*pop->N_SUBS+irep]);
        }
    }
    
   for (isub=0;isub<pop->N_SUBS;isub++) {
        for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
            for (i=0;i<pop->sub[isub].rep[irep].dim_X[0];i++) {
                    pop->sub[isub].rep[irep].md_Y[i] /= (double)(MAX_ITER-BURN_IN);
                    fprintf(fout_sub_precY[isub*pop->N_SUBS+irep],"%lf ",pop->sub[isub].rep[irep].md_Y[i]);
 //                   fprintf(fout_delta2[isub*pop->N_SUBS+irep],"%lf ",sqrt(pop->sub[isub].rep[irep].mdelta2[i*pop->sub[isub].rep[irep].dim_W[1]+j]/(MAX_ITER-BURN_IN)- pop->sub[isub].rep[irep].mdelta[i*pop->sub[isub].rep[irep].dim_W[1]+j]*pop->sub[isub].rep[irep].mdelta[i*pop->sub[isub].rep[irep].dim_W[1]+j]));
             }                
               fprintf(fout_sub_precY[isub*pop->N_SUBS+irep],"\n");                  
  //              fprintf(fout_d_Y2[isub*pop->N_SUBS+irep],"\n");  
            fclose(fout_sub_precY[isub*pop->N_SUBS+irep]);
    //        fclose(fout_delta2[isub*pop->N_SUBS+irep]);
        }
    }
    
printf("here\n");fflush(NULL);
    C = (char *)malloc(3*sizeof(char));
    S = (char *)malloc(100*sizeof(char));
      
      
    for (isub=0;isub<pop->N_SUBS;isub++) {
        for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
            S = strcpy(S,"Wdelta");
            itoa(isub,C);
            S = strcat(S,C);
            itoa(irep,C);
            S = strcat(S,C);
            S = strcat(S,".dat");
            fout = fopen(S,"w");
            for (i=0;i<pop->sub[isub].rep[irep].dim_W[0];i++)
                fprintf(fout,"%lf ",pop->sub[isub].rep[irep].mWdelta[i]/(double)((MAX_ITER-BURN_IN)));
            fclose(fout);
        }
    }
    fout = fopen("lastres.dat","w");
    for (i=0;i<pop->sub[0].rep[0].dim_X[0];i++)
        fprintf(fout,"%lf ",pop->sub[0].rep[0].residuals[i]);
    fprintf(fout,"\n");
    fclose(fout);
    
    for (isub=0;isub<pop->N_SUBS;isub++) {
        for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
            S = strcpy(S,"Veta");
            itoa(isub,C);
            S = strcat(S,C);
            itoa(irep,C);
            S = strcat(S,C);
            S = strcat(S,".dat");
            fout = fopen(S,"w");
            for (i=0;i<pop->sub[isub].rep[irep].dim_V[0];i++)
                fprintf(fout,"%lf ",pop->sub[isub].rep[irep].mVeta[i]/(double)((MAX_ITER-BURN_IN)));
            fclose(fout);
        }
    }
    pop->ED /= ((MAX_ITER-BURN_IN));
    for (isub=0;isub<pop->N_SUBS;isub++) {
        for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
            for (i=0;i<pop->sub[isub].rep[irep].dim_X[0];i++) {
                pop->sub[isub].rep[irep].mean_res[i] /= ((MAX_ITER-BURN_IN));
                pop->sub[isub].rep[irep].mean_fit[i] /= ((MAX_ITER-BURN_IN));
                pop->sub[isub].rep[irep].mean_d_Y[i] /= ((MAX_ITER-BURN_IN));
            }
            pop->sub[isub].rep[irep].mean_precY /= ((MAX_ITER-BURN_IN));
       }
    }
    double DE=0;
    double mlogsqrt2pi = -0.9189385332046727;
    for (isub=0;isub<pop->N_SUBS;isub++) {
        for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
            for (i=pop->sub[isub].rep[irep].P;i<pop->sub[isub].rep[irep].dim_X[0];i++)
                DE += mlogsqrt2pi + sqrt(pop->sub[isub].rep[irep].mean_d_Y[i]) - 0.5*(pop->sub[isub].rep[irep].mean_res[i])*(pop->sub[isub].rep[irep].mean_res[i]*pop->sub[isub].rep[irep].mean_d_Y[i]);
        }
    }
    DE *= -2;
    printf("DIC = %lf pD = %lf\n",2*pop->ED - DE,pop->ED - DE);
    fout = fopen("DIC.dat","w");
    fprintf(fout, "DIC = %15.3lf, pD = %15.3lf\n",2*pop->ED - DE,pop->ED - DE);
    fclose(fout);
    for (isub=0;isub<pop->N_SUBS;isub++) {
        for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
            S = strcpy(S,"mean_res");
            itoa(isub,C);
            S = strcat(S,C);
            itoa(irep,C);
            S = strcat(S,C);
            S = strcat(S,".dat");
            fout = fopen(S,"w");
            for (i=0;i<pop->sub[isub].rep[irep].dim_X[0];i++)
                fprintf(fout,"%lf ",pop->sub[isub].rep[irep].mean_res[i]);
            fclose(fout);
        }
    }

    for (isub=0;isub<pop->N_SUBS;isub++) {
        for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
            S = strcpy(S,"mean_fit");
            itoa(isub,C);
            S = strcat(S,C);
            itoa(irep,C);
            S = strcat(S,C);
            S = strcat(S,".dat");
            fout = fopen(S,"w");
            for (i=0;i<pop->sub[isub].rep[irep].dim_X[0];i++)
                fprintf(fout,"%lf ",pop->sub[isub].rep[irep].mean_fit[i]);
            fclose(fout);
        }
    }

    for (isub=0;isub<pop->N_SUBS;isub++) {
        for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
            S = strcpy(S,"mean_Y");
            itoa(isub,C);
            S = strcat(S,C);
            itoa(irep,C);
            S = strcat(S,C);
            S = strcat(S,".dat");
            fout = fopen(S,"w");
            for (i=0;i<pop->sub[isub].rep[irep].dim_V[0];i++)
                fprintf(fout,"%lf ",pop->sub[isub].rep[irep].mean_Y[i]/=(double)(MAX_ITER-BURN_IN));
            fclose(fout);
        }
    }

    
    free(C);
    free(S);
    
  }

void draw_precYstart(REP *rep,unsigned long *seed) {
    double ALPHA=0.001,BETA=0.001;
    double a=0, b=0;
    
    a = (double)rep->P/2. + ALPHA;
    
    b = 0;
    for (int i=0;i<rep->P;i++)
        b += (rep->Y[i]-rep->Veta[i])*(rep->Y[i]-rep->Veta[i]);
    
    rep->precYstart = rgamma(a,b,seed);
    for (int i=0;i<rep->P;i++)
        rep->d_Y[i] = rep->precYstart;

}

void calH(double *H,double *A,double *delta,const int nrow,const int ncol,const int P) {
    int t,j,k;
    double *tmp;
    
    for (j=0;j<nrow;j++) 
        for (k=0;k<ncol;k++)
            H[j*ncol+k] = 0;
    
    for (t=P;t<nrow;t++) 
        for (j=0;j<P;j++) 
            for (k=0;k<ncol;k++) 
                H[t*ncol+k] += A[(t-j-1)*ncol+k]*delta[t*P+j];
}

void calApVinvX(double *ApX,double *A,double *X,double *Vinv,const int nrow,const int ncol) {
    int i,j;
    
    for (i=0;i<ncol;i++)
        ApX[i] = 0;
        
    for (j=0;j<nrow;j++) 
        for (i=0;i<ncol;i++)
            ApX[i] += A[j*ncol+i]*Vinv[j]*X[j];  
}

void calApVinvA(double *ApA,double *A,double *Vinv,const int nrow,const int ncol) {
// only for diagonal matrices, V
    int i,j,k;
        
    for (i=0;i<ncol;i++)
        for (j=0;j<ncol;j++)
            ApA[i*ncol+j] = 0;
   
     for (k=0;k<nrow;k++) {
        double V = Vinv[k];
        double *tmp = &A[k*ncol];
        for (i=0;i<ncol;i++) {
            double *tmpout = &ApA[i*ncol];
            for (j=0;j<ncol;j++)
                tmpout[j] += V*tmp[i]*tmp[j];
        }
    }
}

void calAx(double *Ax,double *A,double *x,const int nrow,const int ncol) {
    int i,j;
    
    for (i=0;i<nrow;i++) {
        Ax[i] = 0;
        for (j=0;j<ncol;j++)
            Ax[i] += A[i*ncol+j]*x[j];
    }
    
}

void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol) {
    int i,j;
    
     for (i=0;i<nrow;i++) {
        Ax[i] = 0;
        if (i>= ncol) {
        for (j=0;j<ncol;j++)
            Ax[i] += A[i*ncol+j]*x[i*ncol+j];}
    }
    
}

void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P) {
    int t,i,j;
        
    for (i=0;i<nrow;i++) {
        int min = (i < ncol) ? i:ncol;
        for (j=0;j<min;j++) {
            W[i*ncol+j] = Y[i-j-1] - Xb[i-j-1] - Ve[i-j-1];
       }
    }
}

void calculate_residuals(REP *rep,int P) {
    int i,j;
    
    double tmp = 0;
    for (i=0;i<P;i++)
        rep->residuals[i] = rep->Y[i] - rep->Veta[i];
    for (i=P;i<rep->dim_X[0];i++)
        rep->residuals[i] = rep->Y[i] - (rep->Veta[i] + rep->Wdelta[i] + rep->Xbeta[i]);
    
}

void calculate_res1(REP *rep,int P) {
    int i,j;
    
    for (i=0;i<P;i++)
        rep->residuals1[i] = 0;
    for (i=P;i<rep->dim_X[0];i++) {
        rep->residuals1[i] = rep->Y[i] - rep->Wdelta[i] - rep->Veta[i];
    }
    
}

void calculate_res2(REP *rep,int P) {

    for (int i=0;i<P;i++)
        rep->residuals2[i] = 0;
    for (int i=P;i<rep->dim_X[0];i++) {
        rep->residuals2[i] = rep->Y[i] - rep->Wdelta[i] - rep->Xbeta[i] - rep->Veta[i];
    }
}

void calculate_res3(REP *rep,int P) {
//    for (int i=0;i<P;i++)
//        rep->residuals3[i] = 0;

    for (int i=0;i<rep->dim_X[0];i++) {
        rep->residuals3[i] = rep->Y[i] -  rep->Xbeta[i] - rep->Veta[i];
    }
}

void calculate_res4(REP *rep,int P) {
    for (int i=0;i<P;i++) {
        rep->residuals4[i] = rep->Y[i];
    }
   for (int i=P;i<rep->dim_V[0];i++)
        rep->residuals4[i] = rep->Y[i] - rep->Xbeta[i] - rep->Wdelta[i];
}

void calculate_res5(REP *rep,int P) {
    for (int i=0;i<P;i++) {
        rep->residuals5[i] = rep->Y[i];
    }
   for (int i=P;i<rep->dim_V[0];i++)
        rep->residuals5[i] = rep->Y[i] - rep->Wdelta[i];
}

void draw_preceta(REP *rep,double **Q,unsigned long *seed) {
    double ALPHA=0, BETA=0;
    double a, b;
    ALPHA = 1;
    BETA =  1;
    b = 0;
    
    for (int i=0;i<rep->dim_V[1];i++)
//        for (int j=0;j<rep->dim_V[1];j++)
            b += rep->eta[i]*rep->eta[i];
 //           b += rep->eta[i]*Q[i][j]*rep->eta[j];
    
    a = 0.5*rep->dim_V[1] + ALPHA;
    
    b = 0.5*b + BETA;
    
    rep->preceta = rgamma(a,b,seed);
}

void draw_reprec(POP *pop,SUB *sub,unsigned long *seed) {
    int N;
    double ALPHA, BETA,tmp;
    
    ALPHA = 1;
    BETA  = 1;
    
    N = pop->N_SUBS;
    for (int is=0;is<pop->Ns;is++) {
        for (int j=0;j<pop->Nb;j++) {
            tmp = 0;
            for (int i=0;i<pop->N_SUBS;i++) {
                tmp += (sub[i].beta[j+is*pop->Nb]-pop->beta[j+is*pop->Nb])*(sub[i].beta[j+is*pop->Nb]-pop->beta[j+is*pop->Nb]);
            }
            double tmp2 = rgamma(ALPHA + 0.5*(double)N,BETA + 0.5*tmp,seed);
            pop->re_prec[j+is*pop->Nb] = tmp2;
        }
    }
}

double tden(double x,double mean,double var,double df) {
    double tmp,value = 0;
    double pi = 3.14159;
    
    tmp = (x-mean)*(x-mean)/var;
    value = lgamma(0.5*(df+1)) - lgamma(0.5*df) - 0.5*log(df*pi*var) - 0.5*(df+1)*log1p(tmp/df);
    return value;
}

void calW2(double *W,double *mres,const int nrow,const int ncol,const int P) {
    int t,i,j;
        
    for (i=0;i<nrow;i++) {
        int min = (i < ncol) ? i:ncol;
        for (j=0;j<min;j++) {
            W[i*ncol+j] = mres[i-j-1];
       }
    }
}

void DLMtst(REP *rep,int iter,int flag,unsigned long *seed) {
    double n, n0=1,S0=1,d;
    double Q,*A;
    double f,e,S_old,*W;
   
    int P = rep->P;
    W = rep->W;
    A = (double *)calloc(P,sizeof(double));

    calculate_res3(rep,P);

    int prop_P;
    double new_log_dens,old_log_dens;
    int sign;
    double u;
    if (flag == 2) {
        double *CDF1,*PDF1;
        CDF1 = (double *)calloc((maxP-1+1),sizeof(double));
        PDF1 = (double *)calloc((maxP-1+1),sizeof(double));
        double lfact = 0;
        for (int i=1;i<=maxP;i++) {
            lfact += log((double)i);
            PDF1[i-1] =  exp(-P + log(pow((double)P,(double)i)) - lfact);
            CDF1[i-1] = PDF1[i-1];
        }
        for (int i=2;i<=maxP;i++)
            CDF1[i-1] += CDF1[i-2];
        
        for (int i=1;i<=maxP;i++)
            PDF1[i-1] /= CDF1[maxP-1];
        for (int i=1;i<=maxP;i++)
            CDF1[i-1] /= CDF1[maxP-1];
        
        prop_P = trunc_rpois(CDF1,1,maxP,seed);
        new_log_dens = log(PDF1[prop_P-1]);
 
        lfact = 0;
        for (int i=1;i<=maxP;i++) {
            lfact += log((double)i);
            PDF1[i-1] =  exp(-prop_P + log(pow((double)prop_P,(double)i)) - lfact);
            CDF1[i-1] = PDF1[i-1];
        }
        for (int i=2;i<=maxP;i++)
            CDF1[i-1] += CDF1[i-2];
        
        for (int i=1;i<=maxP;i++)
            PDF1[i-1] /= CDF1[maxP-1];
        for (int i=1;i<=maxP;i++)
            CDF1[i-1] /= CDF1[maxP-1];
 
        old_log_dens = log(PDF1[P-1]);
  
        free(PDF1);
        free(CDF1);
 //       printf("prop_P = %3d \t P = %3d \t",prop_P,P);fflush(stdout);
    }
    else {
        new_log_dens = old_log_dens = 0;
        prop_P = P;
    }


    int Pmax = (P > prop_P) ? P:prop_P;
    
    // Calculate old log likelihood
    for (int i=0;i<P*P;i++)
        dlmStruc[P-1].C[i] = 0;
    for (int i=0;i<P;i++)
        dlmStruc[P-1].C[i+P*i] = 1;
                   
    S0 = S0/n0;
    
    n = n0;
    d = n*S0;
    dlmStruc[P-1].S = d/n;

    for (int i=0;i<P;i++)
        dlmStruc[P-1].m[i] = 0;
    
    double ll=0;
    for (int t=P;t<rep->dim_X[0];t++) {
        for (int i=0;i<P;i++)
            dlmStruc[t].a[i] = dlmStruc[t-1].m[i];
                
        f = 0;
        for (int i=0;i<P;i++)
            f += W[t*P+i]*dlmStruc[t].a[i];
    
        for (int i=0;i<P;i++)
            for (int j=0;j<P;j++) 
                dlmStruc[t].R[i*P+j] = dlmStruc[t-1].C[i*P+j]/(rep->df_delta1);
 
        Q = dlmStruc[t-1].S;
        for (int i=0;i<P;i++)
            for (int j=0;j<P;j++)
                Q += W[t*P+i]*dlmStruc[t].R[i*P+j]*W[t*P+j];
            
        for (int i=0;i<P;i++) {
            A[i] = 0;
            for (int j=0;j<P;j++)
                A[i] += dlmStruc[t].R[i*P+j]*W[t*P+j];
            A[i] /= Q; 
        }                
                           
        e = rep->residuals3[t] - f;        

        if (t >= Pmax)
            ll += tden(rep->residuals3[t],f,Q,n);
  
        S_old = d/n;
        n =  n*(rep->df_delta2) + 1;
        d = (rep->df_delta2)*d + dlmStruc[t-1].S*e*e/Q;
        dlmStruc[t].S = d/n;
                    
        double tmp = dlmStruc[t].S/dlmStruc[t-1].S;
        for (int i=0;i<P;i++) {
            for (int j=0;j<P;j++) {
                dlmStruc[t].C[i*P+j] = tmp*(dlmStruc[t].R[i*P+j] - A[i]*A[j]*Q);
            }
        }
                
        for (int i=0;i<P;i++)
            dlmStruc[t].m[i] = dlmStruc[t].a[i] + A[i]*e;
    }
    free(A);
    double old_loglik = ll;
    

    // Calculate new log likelihood
    A = (double *)calloc(prop_P,sizeof(double));

    double prop_delta,prop_beta;
    double low,high;
    double range = rep->prop_sd[0];
    low = rep->df_delta1 - range;
    low = (low < 0.8) ? 0.8:low;
    high = rep->df_delta1 + range;
    high = (high > 0.99999) ? 0.99999:high;
    
    prop_delta = runif_atob(seed,low,high);
    
    range = rep->prop_sd[1];
    low = rep->df_delta2 - range;
    low = (low < 0.7) ? 0.7:low;
    high = rep->df_delta2 + range;
    high = (high > 0.99999) ? 0.99999:high;
    
    prop_beta = runif_atob(seed,low,high);
     
    if (flag == 0) {
        prop_beta = rep->df_delta2;
    }   
    else if (flag == 1)
        prop_delta = rep->df_delta1;
    else if (flag == 2) {
        prop_delta = rep->df_delta1;
        prop_beta = rep->df_delta2;
    }
        
        
    for (int i=0;i<prop_P*prop_P;i++)
        dlmStruc[prop_P-1].C[i] = 0;
    for (int i=0;i<prop_P;i++)
        dlmStruc[prop_P-1].C[i+prop_P*i] = 1;
 
    S0 = 1;
    n0 = 1; 
    S0 = S0/n0;
    
    n = n0;
    d = n*S0;
    dlmStruc[prop_P-1].S = d/n;

    for (int i=0;i<prop_P;i++)
        dlmStruc[prop_P-1].m[i] = 0;
    
    if (prop_P != P) {
        calW(rep->W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)prop_P,P);
        W = rep->W;
    }
    ll=0;
    for (int t=prop_P;t<rep->dim_X[0];t++) {
        for (int i=0;i<prop_P;i++)
            dlmStruc[t].a[i] = dlmStruc[t-1].m[i];
                
        f = 0;
        for (int i=0;i<prop_P;i++)
            f += W[t*prop_P+i]*dlmStruc[t].a[i];
    
        for (int i=0;i<prop_P;i++)
            for (int j=0;j<prop_P;j++) 
                dlmStruc[t].R[i*prop_P+j] = dlmStruc[t-1].C[i*prop_P+j]/prop_delta;
 
        Q = dlmStruc[t-1].S;
        for (int i=0;i<prop_P;i++)
            for (int j=0;j<prop_P;j++)
                Q += W[t*prop_P+i]*dlmStruc[t].R[i*prop_P+j]*W[t*prop_P+j];
            
        for (int i=0;i<prop_P;i++) {
            A[i] = 0;
            for (int j=0;j<prop_P;j++)
                A[i] += dlmStruc[t].R[i*prop_P+j]*W[t*prop_P+j];
            A[i] /= Q; 
        }                
                           
        e = rep->residuals3[t] - f;        

        if (t >= Pmax)
            ll += tden(rep->residuals3[t],f,Q,n);

        S_old = d/n;
        n =  n*prop_beta + 1;
        d = prop_beta*d + dlmStruc[t-1].S*e*e/Q;
        dlmStruc[t].S = d/n;
                    
        double tmp = dlmStruc[t].S/dlmStruc[t-1].S;
        for (int i=0;i<prop_P;i++) {
            for (int j=0;j<prop_P;j++) {
                dlmStruc[t].C[i*prop_P+j] = tmp*(dlmStruc[t].R[i*prop_P+j] - A[i]*A[j]*Q);
            }
        }
                
        for (int i=0;i<prop_P;i++)
            dlmStruc[t].m[i] = dlmStruc[t].a[i] + A[i]*e;
    }
    double new_loglik = ll;
// if (flag==2)   printf("log_prop_ratio = %lf \n",new_loglik - old_loglik- new_log_dens + old_log_dens); 
    if (log(kiss(seed)) < (new_loglik - old_loglik) - (new_log_dens - old_log_dens)) {
        rep->df_delta1 = prop_delta;
        rep->df_delta2 = prop_beta;
        rep->P = prop_P;
        rep->dim_W[1] = rep->P;
        if (flag < 2)
            (rep->accept[flag])++;
    }
    if (flag < 2)
        (rep->attempt[flag])++;
        
    calW(rep->W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1],rep->P);

// if(flag==2)   printf("delta = %lf \t beta = %lf \t P = %d\n",*delta,*beta,*PP);
    free(A);
}

void DLM(REP *rep,int P,unsigned long *seed) {
    double n, n0=1,S0=1,d;
    double *C0,Q,*A,*Var;
    double f,e,S_old;
    int err;
    clock_t  start,end;
    
    start = clock();
    
    A = (double *)calloc(P,sizeof(double));
    C0 = (double *)calloc(P*P,sizeof(double));
//     for (int i=0;i<P;i++)
//        C0[i] = (double *)calloc(P,sizeof(double));
    for (int i=0;i<P;i++)
        C0[i+P*i] = 1.0;//0.01;
        
    Var = (double *)calloc(P*P,sizeof(double));
//    for (int i=0;i<P;i++)
//        Var[i] = (double *)calloc(P,sizeof(double));
 
    for (int i=0;i<P;i++)
        for (int j=0;j<P;j++)
            dlmStruc[P-1].C[i*P+j] = 0;
    for (int i=0;i<P;i++)
            dlmStruc[P-1].C[i*P+i] = C0[i+P*i];
                   
    S0 = S0/n0;
    
    n = n0;
    d = n*S0;
    dlmStruc[P-1].S = d/n;

 /*   m0[0] = 0.5;
    m0[1] = 0.2;
    m0[2] = 0.1;*/
    for (int i=0;i<P;i++)
        dlmStruc[P-1].m[i] = 0;
    
    calculate_res3(rep,P);
    
    double ll=0;
    for (int t=P;t<rep->dim_X[0];t++) {
        for (int i=0;i<P;i++)
            dlmStruc[t].a[i] = dlmStruc[t-1].m[i];
             
        f = 0;
        for (int i=0;i<P;i++)
            f += rep->W[t*P+i]*dlmStruc[t].a[i];
    
        
        for (int i=0;i<P;i++)
            for (int j=0;j<P;j++) 
                dlmStruc[t].R[i*P+j] = dlmStruc[t-1].C[i*P+j]/rep->df_delta1;
 
        Q = dlmStruc[t-1].S;
        for (int i=0;i<P;i++)
            for (int j=0;j<P;j++)
                Q += rep->W[t*P+i]*dlmStruc[t].R[i*P+j]*rep->W[t*P+j];
        
        for (int i=0;i<P;i++) {
            A[i] = 0;
            for (int j=0;j<P;j++)
                A[i] += dlmStruc[t].R[i*P+j]*rep->W[t*P+j];
            A[i] /= Q; 
        }                
                           
        e = rep->residuals3[t] - f;        
        
        S_old = d/n;
        n =  n*rep->df_delta2 + 1;
        d = rep->df_delta2*d + dlmStruc[t-1].S*e*e/Q;
        dlmStruc[t].S = d/n;
                    
        double tmp = dlmStruc[t].S/dlmStruc[t-1].S;
        for (int i=0;i<P;i++) {
            for (int j=0;j<P;j++) {
                dlmStruc[t].C[i*P+j] = tmp*(dlmStruc[t].R[i*P+j] - A[i]*A[j]*Q);
                Var[i*P+j] = dlmStruc[t].C[i*P+j];
            }
        }
                
        for (int i=0;i<P;i++)
            dlmStruc[t].m[i] = dlmStruc[t].a[i] + A[i]*e;
 
        if (t == rep->dim_X[0]-1) {
            err = cholesky_decomp2vec(Var,P);
            if (err) {  // err = 1 means P is SPD
                err = rmvtvec(&(rep->delta[t*P]),Var,P,n,dlmStruc[t].m,seed);
            }
            else {
                printf("error in DLM, var is not SPD\n");
                exit(0);
            }
            rep->d_Y[t] = rgamma(n/2.,d/2.,seed);
        }     
    }
//    printf("ll = %lf\n",ll);
//    fprintf(fout,"%lf ",ll);
    double *C,*m;
    m = (double *)calloc(P,sizeof(double));
//    C = (double *)calloc(P*P,sizeof(double));
//    for (int i=0;i<P;i++)
//        C[i] = (double *)calloc(P,sizeof(double));
        
    n = rep->dim_X[0];
    for (int t=rep->dim_X[0]-2;t>=P;t--) {
        
        for (int i=0;i<P;i++)
            m[i] = dlmStruc[t].m[i] + rep->df_delta1*(rep->delta[(t+1)*P + i] - dlmStruc[t+1].a[i]);


        // compute mean
/*        err = cholesky_decomp2vec(dlmStruc[t+1].R,P);
        
        for (int i=0;i<P;i++)
            mean[i] = rep->delta[(t+1)*P + i] - dlmStruc[t+1].a[i];
        if (err) {  // err = 1 means P is SPD
            err = forward_substitution2vec(dlmStruc[t+1].R,mean,P);
            err = cholesky_backsub2vec(dlmStruc[t+1].R,mean,P);
        }
        else {printf("backward smoothing error\n");exit(0);}
    
        for (int i=0;i<P;i++) {
            m[i] = dlmStruc[t].m[i];
            for (int j=0;j<P;j++)
                m[i] += dlmStruc[t].C[i*P+j]*mean[j];
        }    */
        
        // compute var
 /*       for (int i=0;i<P;i++) {
            for (int j=0;j<P;j++) 
                mean[j] = dlmStruc[t].C[i*P+j];
            
            err = forward_substitution2vec(dlmStruc[t+1].R,mean,P);
            err = cholesky_backsub2vec(dlmStruc[t+1].R,mean,P);
            
            for (int j=0;j<P;j++) {
                C[i*P+j] = dlmStruc[t].C[i*P+j];
                Var[i*P+j] = mean[j];
            }
        }
        
        for (int i=0;i<P;i++)
            for (int j=0;j<P;j++)
                for (int k=0;k<P;k++)
                   C[i*P+j] -= dlmStruc[t].C[i*P+k]*Var[j*P+k];
                                            
    */
        
/*        double delta2 = delta*delta;
        for (int i=0;i<P;i++) {
            for (int j=0;j<P;j++) {
                C[i*P+j] = dlmStruc[t].C[i*P+j] - delta2*dlmStruc[t+1].R[i*P+j];
            }
        }       
        
        double tmp = dlmStruc[t+1].S/dlmStruc[t].S;
        for (int i=0;i<P;i++)
            for (int j=0;j<P;j++)
                C[i*P+j] *= tmp;*/
                
        double tmp = 1-rep->df_delta1;//dlmStruc[t].S/dlmStruc[t-1].S*(1-delta); // THIS MAY BE WRONG, now corrected
// printf("%lf \n",dlmStruc[t].S/dlmStruc[t-1].S);
        for (int i=0;i<P;i++)
            for (int j=0;j<P;j++)
                Var[i*P+j] = dlmStruc[t].C[i*P+j]*tmp;
                
               
        // draw new delta[t]
        err = cholesky_decomp2vec(Var,P);
        if (err) {  // err = 1 means C is SPD
            err = rmvtvec(&(rep->delta[t*P]),Var,P,n,m,seed);
        }
        else {
           printf("error in DLM last loop, var is not SPD\n");
           exit(0);
        }
        
        // draw new prec[t]
        rep->d_Y[t] = rep->df_delta2*rep->d_Y[t+1];
        n = n - 1;    
        rep->d_Y[t] += rgamma(0.5*(1-rep->df_delta2)*n,0.5*n*dlmStruc[t].S,seed);
    }
    
    free(m);
//    free(C);
    
    calWdelta(rep->Wdelta,rep->W,rep->delta,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
    calculate_residuals(rep,P);

//     for (int i=0;i<P;i++)
//        free(C0[i]);
    free(C0);
    free(Var);
    free(A);
    
       end = clock();
    
 //   printf("%.6lf \n",double(end-start)/CLOCKS_PER_SEC);

}

void draw_beta_eta(POP *pop,SUB *sub,REP *rep,unsigned long *seed) {
    int i,j,ncol;
    double *M,*YY,*V,*VpV,*J,*X,*betaeta;
        
    /* calculate H(X,delta) */
    /* calculate J(V,delta) */
  
    ncol = (rep->dim_X[1] + rep->dim_V[1]);
    
    betaeta = (double *)calloc(ncol,sizeof(double));

    V = (double *)calloc(rep->dim_V[0]*ncol,sizeof(double));
//    for (i=0;i<rep->dim_V[0];i++)
//        V[i] = (double *)calloc(ncol,sizeof(double));
    
    J = (double *)calloc(rep->dim_V[0]*ncol,sizeof(double));
//    for (i=0;i<rep->dim_V[0];i++)
//        J[i] = (double *)calloc(ncol,sizeof(double));
   
    X = (double *)calloc(rep->dim_V[0]*ncol,sizeof(double));
//    for (i=0;i<rep->dim_V[0];i++)
//        X[i] = (double *)calloc(ncol,sizeof(double));

    for (i=0;i<rep->dim_V[0];i++)
        for (j=0;j<rep->dim_X[1];j++)
            X[i*ncol+j] = rep->X[i*rep->dim_X[1]+j];
     for (i=0;i<rep->dim_V[0];i++)
        for (j=rep->dim_X[1];j<ncol;j++)
            X[i*ncol+j] = rep->V[i*rep->dim_V[1]+j-rep->dim_X[1]];
    int cnt = 0;
    for (i=0;i<rep->dim_X[1];i++)
        betaeta[i] = sub->beta[i];
    for (i=rep->dim_X[1];i<ncol;i++)
        betaeta[i] = rep->eta[i-rep->dim_X[1]];
 
     calH(J,X,rep->delta,(const int)rep->dim_V[0],(const int)ncol,(const int)rep->P);
  
    /* calculate H(X,delta)beta */
    /* calculate J(V,delta)eta */
    double *Jbe;
    Jbe = (double *)calloc(rep->dim_V[0],sizeof(double));
    calAx(Jbe,J,betaeta,(const int)rep->dim_V[0],(const int)ncol);
    
    free(betaeta);
    /* calculate Y - Wdelta - J(V,delta)eta - H(X,delta)beta*/
    
    calculate_res5(rep,rep->P);

    for (i=0;i<rep->dim_V[0];i++)
        rep->residuals5[i] -= Jbe[i];
       
    free(Jbe);
   
    M = (double *)calloc(ncol,sizeof(double));
    VpV = (double *)calloc(ncol*ncol,sizeof(double));
//    for (i=0;i<ncol;i++)
//        VpV[i] = (double *)calloc(ncol,sizeof(double));
   
  
    /* Add X and -J(V,delta)*/
    
    for (i=0;i<rep->dim_V[0];i++)
        for (j=0;j<ncol;j++)
            V[i*ncol+j] = X[i*ncol+j] - J[i*ncol+j];

            
    calApVinvX(M,V,rep->residuals5,rep->d_Y,rep->dim_V[0],ncol);
   
    /* calculate  (X - J)'(X - J) */
    
    calApVinvA(VpV,V,rep->d_Y,(const int)rep->dim_V[0],(const int)ncol);
   
    for (i=rep->dim_X[1];i<ncol;i++)
        VpV[i*ncol+i] += rep->preceta;

    if (pop->non_parm || pop->GRP) {
        for (int is=0;is<pop->Ns;is++) {
            for (int i=0;i<pop->Nb;i++) {
//                for (int j=0;j<pop->Nb;j++) {
                    M[i+is*pop->Nb] += pop->re_prec[j+is*pop->Nb]*pop->beta[i+is*pop->Nb];
                    VpV[(i+is*pop->Nb)*ncol+i+is*pop->Nb] += pop->re_prec[j+is*pop->Nb];
//                }
            }
        }
    }
   
    int err = cholesky_decomp2vec(VpV,ncol);
    
    double *mean;
    mean = (double *)calloc(ncol,sizeof(double));
   
   if (err) {  // err = 1 means P is SPD
        err = forward_substitution2vec(VpV,M,ncol);
        err = cholesky_backsub2vec(VpV,M,ncol);
        err = rmvnorm3vec(mean,VpV,ncol,M,seed,1);
    }
    else {
        printf("error in draw_beta_eta, precision is not SPD\n");
        exit(0);
    }

    /*  Will need to change this so that each rep has its own beta, for cases with multiple reps per subject. And r.e. about sub->beta and r.e. about pop->beta. CODE changes required */

    for (i=0;i<rep->dim_X[1];i++)
        sub->beta[i] = mean[i];
    for (i=rep->dim_X[1];i<ncol;i++)
        rep->eta[i-rep->dim_X[1]] = mean[i];
        
    /* calculate Veta and Xbeta */
    
    calAx(rep->Veta,rep->V,rep->eta,(const int)rep->dim_V[0],(const int)rep->dim_V[1]);
    calAx(rep->Xbeta,rep->X,sub->beta,(const int)rep->dim_X[0],(const int)rep->dim_X[1]);

    /* update W and Wdelta */
    calW(rep->W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1],rep->P);
    calWdelta(rep->Wdelta,rep->W,rep->delta,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);

//     for (i=0;i<rep->dim_V[0];i++)
//        free(J[i]);
    free(J);
 //   for (i=0;i<rep->dim_V[0];i++)
 //       free(X[i]);
    free(X);
//    for (i=0;i<ncol;i++)
//        free(VpV[i]);
    free(VpV);
    free(M);
    free(mean);
 //   for (i=0;i<rep->dim_V[0];i++)
 //       free(V[i]);
    free(V);
  }

void draw_sub_beta(POP *pop,SUB *sub,unsigned long *seed) {
    double *X,*P,*M,*mean;
    REP *rep;
    
    int N;
    N = pop->Ns*pop->Nb;
    mean = (double *)calloc(N,sizeof(double));
    M = (double *)calloc(N,sizeof(double));
    P = (double *)calloc(N*N,sizeof(double));
            
    for (int i=0;i<sub->N_REPS;i++) {
        rep = &(sub->rep[i]);
        /* calculate H(X,delta) */
 
        calH(rep->H,rep->X,rep->delta,(const int)rep->dim_X[0],(const int)rep->dim_X[1],(const int)sub->rep[i].P);

        /* calculate H(X,delta)beta */
    
        calAx(rep->Hbeta,rep->H,sub->beta,(const int)rep->dim_X[0],(const int)rep->dim_X[1]);

        /* calculate Y - Veta - Wdelta - H(X,delta)beta */
    
        calculate_res1(rep,sub->rep[i].P);
    
        for (int j=sub->rep[i].P;j<rep->dim_X[0];j++)
            rep->residuals1[j] -= rep->Hbeta[j];

        /* Subtract H(X,delta) from X */
        X = (double *)calloc(rep->dim_X[0]*rep->dim_X[1],sizeof(double));
 //       for (int j=0;j<rep->dim_X[0];j++)
 //           X[j] = (double *)calloc(rep->dim_X[1],sizeof(double));
         

        for (int j=sub->rep[i].P;j<rep->dim_X[0];j++)
            for (int k=0;k<rep->dim_X[1];k++)
                X[j*rep->dim_X[1]+k] = rep->X[j*rep->dim_X[1]+k] - rep->H[j*rep->dim_X[1]+k];
    
        calApVinvX(M,X,rep->residuals1,rep->d_Y,rep->dim_X[0],rep->dim_X[1]);
    
        for (int j=0;j<rep->dim_X[1];j++)
            mean[j] = M[j];
            
        /* calculate  (X - H)'V^(-1)(X - H) */
    
        calApVinvA(rep->XpX,X,rep->d_Y,(const int)rep->dim_X[0],(const int)rep->dim_X[1]);
   
 //       for (int j=0;j<rep->dim_X[0];j++)
 //           free(X[j]);
        free(X);
    
    
        for (int j=0;j<rep->dim_X[1];j++)
            for (int k=0;k<rep->dim_X[1];k++)
                P[j*rep->dim_X[1]+k] += rep->XpX[j*rep->dim_X[1]+k];
    }
    
    if (pop->non_parm || pop->GRP) {
        for (int is=0;is<pop->Ns;is++) {
            for (int i=0;i<pop->Nb;i++) {
                for (int j=0;j<pop->Nb;j++) {
                    mean[i+is*pop->Nb] += pop->re_prec[j+is*pop->Nb]*pop->Q[i+is*pop->Nb][j+is*pop->Nb]*pop->beta[j+is*pop->Nb];
                    P[(i+is*pop->Nb)*rep->dim_X[1] + (j+is*pop->Nb)] += pop->re_prec[j+is*pop->Nb]*pop->Q[i+is*pop->Nb][j+is*pop->Nb];
                }
            }
        }
    }
    int err = cholesky_decomp2vec(P,rep->dim_X[1]);
    
    if (err) {  // err = 1 means P is SPD
        err = forward_substitution2vec(P,mean,N);
        err = cholesky_backsub2vec(P,mean,N);
        err = rmvnorm3vec(sub->beta,P,N,mean,seed,1);
    }
    else {
        printf("error in draw_sub_beta, precision is not SPD\n");
        exit(0);
    }
    
    free(P);
    free(M);
    free(mean);
    
    for (int i=0;i<sub->N_REPS;i++) {
        rep = &(sub->rep[i]);
        
      /* calculate Xbeta */
        calAx(rep->Xbeta,rep->X,sub->beta,(const int)rep->dim_X[0],(const int)rep->dim_X[1]);

        /* update W and Wdelta */
    
        calW(rep->W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1],(const int)rep->P);
        calWdelta(rep->Wdelta,rep->W,rep->delta,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);

        //    printf("exiting sub_draw_beta\n");fflush(NULL);
    }
}


void draw_pop_beta(POP *pop,SUB *sub,unsigned long *seed) {
    double *P,*M;
    
    int N = pop->Ns*pop->Nb;
    M = (double *)calloc(N,sizeof(double));
    P = (double *)calloc(N*N,sizeof(double));
            
    for (int is=0;is<pop->Ns;is++) {
        for (int isub=0;isub<pop->N_SUBS;isub++) {
            for (int j=0;j<pop->Nb;j++) {
 //               for (int k=0;k<pop->Nb;k++) {
                    M[j+is*pop->Nb] += pop->re_prec[j+is*pop->Nb]*sub[isub].beta[j+is*pop->Nb];
                    P[(j+is*pop->Nb)*N + (j+is*pop->Nb)] += pop->re_prec[j+is*pop->Nb];
 //               }
            }
        }
    }
    for (int i=0;i<N;i++)
        P[i*N+i] += 0.01;
    int err = cholesky_decomp2vec(P,N);
    
    if (err) {  // err = 1 means P is SPD
        err = forward_substitution2vec(P,M,N);
        err = cholesky_backsub2vec(P,M,N);
        err = rmvnorm3vec(pop->beta,P,N,M,seed,1);
    }
    else {
        printf("error in draw_pop_beta, precision is not SPD\n");
        exit(0);
    }
    
    free(P);
    free(M);
}


void DIC(POP *pop,REP *rep,unsigned long *seed) {
    int i;
    double tmp;
    double mlogsqrt2pi = -0.9189385332046727;
    
    
    tmp = 0;
    for (i=rep->P;i<rep->dim_X[0];i++)
            tmp += mlogsqrt2pi + sqrt(rep->d_Y[i]) - 0.5*rep->d_Y[i]*(rep->residuals[i])*(rep->residuals[i]);

    pop->ED += -2*tmp;
}
