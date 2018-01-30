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

extern int MAX_ITER;
extern int BURN_IN;
extern sDLM *dlmStruc;
extern int maxP;
extern FILE *flog;

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
    double **loglik;
    int nrow_loglik;
    FILE *fout;
    
    void draw_knot_locations(REP *rep,int sdegree,int *flag,unsigned long *seed);
    int knot_birth_death(REP *rep,POP *pop,const int sdegree,int iter,unsigned long *seed);
    void itoa(int n,char s[]);
 
    void DLMtst(REP *rep,int iter,int flag,unsigned long *seed);
    void DLM(REP *rep,int P,unsigned long *seed);
    
    void draw_precYstart(REP *rep,unsigned long *seed);
    void draw_preceta(REP *rep,unsigned long *seed);
    void draw_reprec(POP *pop,SUB *sub,unsigned long *seed);
 
    void calculate_residuals(REP *rep,int P);
    void calculate_res1(REP *rep,int P);
    void calculate_res2(REP *rep,int P);
    void calculate_res3(REP *rep,int P);
    void calculate_res4(REP *rep,int P);
    void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    
    void draw_beta_eta(POP *pop,SUB *sub,REP *rep,unsigned long *seed,int iter);
    void draw_pop_beta(POP *pop,SUB *sub,unsigned long *seed);
    
    void DIC(POP *pop,REP *rep,unsigned long *seed);

    pop->ED = 0;
    int reps = 0;
    for (isub=0;isub<pop->N_SUBS;isub++)
            reps += pop->sub[isub].N_REPS;
       
    int sdegree = 4;
    int aaa = 0;

    for (iter=0;iter<=MAX_ITER;iter++) {
        if (!(iter%100)) printf("%d",iter);fflush(stdout);
        if (!(iter%20) && (iter%100)) printf(".");fflush(stdout);
        for (isub=0;isub<pop->N_SUBS;isub++) {
 
            for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
//printf("A\n");fflush(stdout);
                if (iter == 0)
                    draw_beta_eta(pop,&(pop->sub[isub]),&(pop->sub[isub].rep[irep]),seed,iter);
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
//                        draw_preceta(&(pop->sub[isub].rep[irep]),seed);
//                        draw_precYstart(&(pop->sub[isub].rep[irep]),pop->P,seed);
                    }
                }
                draw_beta_eta(pop,&(pop->sub[isub]),&(pop->sub[isub].rep[irep]),seed,iter);
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
                draw_preceta(&(pop->sub[isub].rep[irep]),seed);
                draw_precYstart(&(pop->sub[isub].rep[irep]),seed);
        
                
                
                fprintf(pop->sub[isub].rep[irep].fout_nknots,"%d ",pop->sub[isub].rep[irep].nKnots-8);

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
                    
                   
                   for (i=0;i<pop->sub[isub].rep[irep].dim_V[1];i++)
                        fprintf(pop->sub[isub].rep[irep].fout_eta,"%lf ",pop->sub[isub].rep[irep].eta[i]);
                    fprintf(pop->sub[isub].rep[irep].fout_eta,"\n");
                   
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
                    
                    fprintf(pop->sub[isub].rep[irep].fout_dlm,"%lf %lf %d\n",pop->sub[isub].rep[irep].df_delta1,pop->sub[isub].rep[irep].df_delta2,pop->sub[isub].rep[irep].P);
                    
                    for (i=4;i<pop->sub[isub].rep[irep].nKnots-4;i++)
                        fprintf(pop->sub[isub].rep[irep].fout_knots,"%lf ",pop->sub[isub].rep[irep].knots[i]);
                        
                 }
               
                if (!(iter%100)) {
                    fprintf(flog,"iter = %6d\t Sub %d, Rep %d \t int knots = %d",iter,isub,irep,pop->sub[isub].rep[irep].nKnots-2*4);
                    fprintf(flog,"\t precY0 = %10.6lf \t preceta = %10.6lf\n",pop->sub[isub].rep[irep].d_Y[pop->sub[isub].rep[irep].dim_X[0]-1],pop->sub[isub].rep[irep].preceta);
                    fprintf(flog,"\t df_delta1 = %10.6lf \t df_delta2 = %10.6lf \t P = %3d\n",pop->sub[isub].rep[irep].df_delta1,pop->sub[isub].rep[irep].df_delta2,pop->sub[isub].rep[irep].P);
                  //  printf("\t\t precdelta = %10.6lf, preceta = %10.6lf\n\n",pop->sub[isub].rep[irep].precdelta,pop->sub[isub].rep[irep].preceta);
                    
                 }
                 if (!(iter%100) && (iter <= BURN_IN) && (iter > 1)) {
                    double rt;
//                    if (iter > 2000) {
                    rt = (double)pop->sub[isub].rep[irep].accept[0]/(double)pop->sub[isub].rep[irep].attempt[0];
                    adjust_acceptance2(rt,&(pop->sub[isub].rep[irep].prop_sd[0]),0.35);
                    pop->sub[isub].rep[irep].accept[0] = pop->sub[isub].rep[irep].attempt[0] = 0;
                    fprintf(flog,"\t\t rt = %lf proposal = %lf\n",rt,pop->sub[isub].rep[irep].prop_sd[0]);
//                    }
                    rt = (double)pop->sub[isub].rep[irep].accept[1]/(double)pop->sub[isub].rep[irep].attempt[1];
                    adjust_acceptance2(rt,&(pop->sub[isub].rep[irep].prop_sd[1]),0.35);
                    pop->sub[isub].rep[irep].accept[1] = pop->sub[isub].rep[irep].attempt[1] = 0;
                    fprintf(flog,"\t\t rt = %lf proposal = %lf\n",rt,pop->sub[isub].rep[irep].prop_sd[1]);
   
                    rt = (double)pop->sub[isub].rep[irep].accept[3]/(double)pop->sub[isub].rep[irep].attempt[3];
                    adjust_acceptance2(rt,&(pop->sub[isub].rep[irep].prop_sd[3]),0.35);
                    pop->sub[isub].rep[irep].accept[3] = pop->sub[isub].rep[irep].attempt[3] = 0;
                    fprintf(flog,"\t\t rt = %lf proposal = %lf\n",rt,pop->sub[isub].rep[irep].prop_sd[3]);
                 }
                  fflush(NULL);
            }
         
            if ((iter>BURN_IN)) {
                for (i=0;i<pop->Nb*(pop->Ns);i++)
                    fprintf(pop->sub[isub].fout_beta,"%lf ",pop->sub[isub].beta[i]);
                fprintf(pop->sub[isub].fout_beta,"\n");
            }
        }
        if (pop->non_parm || pop->GRP)
            draw_reprec(pop,pop->sub,seed);
        if (pop->GRP) {
            draw_pop_beta(pop,pop->sub,seed);
        }
        if (!(iter%100)) {
            double nmean[25];
            for (k=0;k<pop->Ns;k++)
                fprintf(flog,"Stimulus %d\n",k);
            if (pop->GRP) {
                for (i=0;i<pop->Ns;i++)
                    fprintf(flog,"\t %g ",pop->re_prec[i]);
                fprintf(flog,"\n\n");
            }
            for (i=0;i<pop->N_SUBS;i++) {
                for (j=0;j<pop->Nb*pop->Ns;j++) {
                    fprintf(flog,"\t %10.6lf \t ",pop->sub[i].beta[j]);
                }
                fprintf(flog,"\n");
            }
            fprintf(flog,"\n");
            if (pop->GRP) {
                for (i=0;i<pop->Ncov*pop->Nb*pop->Ns;i++)
                    fprintf(flog,"\t %10.6lf \t ",pop->beta[i]);
                fprintf(flog,"\n\n");fflush(stdout);
            }
        }

        if ((iter>BURN_IN)) {
            for (i=0;i<pop->Ncov*pop->Nb*pop->Ns;i++)
                fprintf(pop->fout_beta,"%lf ",pop->beta[i]);
            fprintf(pop->fout_beta,"\n");
            for (i=0;i<pop->Ns;i++)
                fprintf(pop->fout_reprec,"%lf ",pop->re_prec[i]);
            fprintf(pop->fout_reprec,"\n");fflush(NULL);
            
        }
    }
    fclose(pop->fout_beta);
    fclose(pop->fout_reprec);

    for (isub=0;isub<pop->N_SUBS;isub++) {
        for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
            fclose(pop->sub[isub].rep[irep].fout_dlm);
            fclose(pop->sub[isub].rep[irep].fout_eta);
            fclose(pop->sub[isub].rep[irep].fout_nknots);
            fclose(pop->sub[isub].rep[irep].fout_knots);
        }
        fclose(pop->sub[isub].fout_beta);
    }

   for (isub=0;isub<pop->N_SUBS;isub++) {
        for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
            for (i=pop->sub[isub].rep[irep].P;i<pop->sub[isub].rep[irep].dim_W[0];i++) {
                for (j=0;j<pop->sub[isub].rep[irep].dim_W[1];j++) {
                    pop->sub[isub].rep[irep].mdelta[i*pop->sub[isub].rep[irep].dim_W[1]+j] /= (double)(MAX_ITER-BURN_IN);
                    fprintf(pop->sub[isub].rep[irep].fout_delta,"%lf ",pop->sub[isub].rep[irep].mdelta[i*pop->sub[isub].rep[irep].dim_W[1]+j]);
                }
                fprintf(pop->sub[isub].rep[irep].fout_delta,"\n");                  
            }                
            fclose(pop->sub[isub].rep[irep].fout_delta);
        }
    }
    
   for (isub=0;isub<pop->N_SUBS;isub++) {
        for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
            for (i=0;i<pop->sub[isub].rep[irep].dim_X[0];i++) {
                    pop->sub[isub].rep[irep].md_Y[i] /= (double)(MAX_ITER-BURN_IN);
                    fprintf(pop->sub[isub].rep[irep].fout_prec,"%lf ",pop->sub[isub].rep[irep].md_Y[i]);
            }                
            fprintf(pop->sub[isub].rep[irep].fout_prec,"\n");                  
            fclose(pop->sub[isub].rep[irep].fout_prec);
        }
    }     
      
    for (isub=0;isub<pop->N_SUBS;isub++) {
        for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
            for (i=0;i<pop->sub[isub].rep[irep].dim_W[0];i++)
                fprintf(pop->sub[isub].rep[irep].fout_wdelta,"%lf ",pop->sub[isub].rep[irep].mWdelta[i]/(double)((MAX_ITER-BURN_IN)));
            fclose(pop->sub[isub].rep[irep].fout_wdelta);
        }
    }
    
    for (isub=0;isub<pop->N_SUBS;isub++) {
        for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
           for (i=0;i<pop->sub[isub].rep[irep].dim_V[0];i++)
                fprintf(pop->sub[isub].rep[irep].fout_veta,"%lf ",pop->sub[isub].rep[irep].mVeta[i]/(double)((MAX_ITER-BURN_IN)));
            fclose(pop->sub[isub].rep[irep].fout_veta);
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
    fprintf(flog,"DIC = %lf pD = %lf\n",2*pop->ED - DE,pop->ED - DE);
    fout = fopen("./log/DIC.log","w");
    fprintf(fout, "DIC = %15.3lf, pD = %15.3lf\n",2*pop->ED - DE,pop->ED - DE);
    fclose(fout);
    
    for (isub=0;isub<pop->N_SUBS;isub++) {
        for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
            for (i=0;i<pop->sub[isub].rep[irep].dim_X[0];i++)
                fprintf(pop->sub[isub].rep[irep].fout_res,"%lf ",pop->sub[isub].rep[irep].mean_res[i]);
            fclose(pop->sub[isub].rep[irep].fout_res);
        }
    }

    for (isub=0;isub<pop->N_SUBS;isub++) {
        for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
           for (i=0;i<pop->sub[isub].rep[irep].dim_X[0];i++)
                fprintf(pop->sub[isub].rep[irep].fout_fit,"%lf ",pop->sub[isub].rep[irep].mean_fit[i]);
            fclose(pop->sub[isub].rep[irep].fout_fit);
        }
    }
        
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

void draw_preceta(REP *rep,unsigned long *seed) {
    double ALPHA=0, BETA=0;
    double a, b;
    ALPHA = 1.1;
    BETA =  1.1;
    b = 0;
    
    for (int i=0;i<rep->dim_V[1];i++)
            b += rep->eta[i]*rep->eta[i];
    
    a = 0.5*rep->dim_V[1] + ALPHA;
    
    b = 0.5*b + BETA;
    
    rep->preceta = rgamma(a,b,seed);
}

void draw_reprec(POP *pop,SUB *sub,unsigned long *seed) {
    int N,M;
    double ALPHA, BETA,tmp;
    double *Xb;

    ALPHA = 1.1;//10;
    BETA  = 1.1;//10;
        
    Xb = (double *)calloc(pop->Ns*pop->Nb,sizeof(double));
    for (int is=0;is<pop->Ns;is++) {
        tmp = 0;
        for (int i=0;i<pop->N_SUBS;i++) {
            M = sub[i].dim_X[0];
            N = sub[i].dim_X[1];
            calAx(Xb,sub[i].X,pop->beta,(const int)M,(const int)N);
            for (int j=0;j<M;j++)
                tmp += (sub[i].beta[j]-Xb[j])*(sub[i].beta[j]-Xb[j]);
        }
        double tmp2 = rgamma(ALPHA + 0.5*(double)pop->N_SUBS*pop->Nb,BETA + 0.5*tmp,seed);
        pop->re_prec[is] = tmp2;
//        printf("%lf\n",pop->re_prec[is]);
    }
    free(Xb);
}

/*void draw_reprec(POP *pop,SUB *sub,unsigned long *seed) {
    int N,M;
    double ALPHA, BETA,*tmp;
    double *Xb;
    ALPHA = 1.1;
    BETA  = 1.1;
    
    M = pop->Ns*pop->Nb;
    N = pop->Ncov*M;
    Xb = (double *)calloc(M,sizeof(double));
    tmp = (double *)calloc(M,sizeof(double));
    

    for (int isub=0;isub<pop->N_SUBS;isub++) {
        calAx(Xb,sub[isub].X,pop->beta,(const int)M,(const int)N);
        for (int i=0;i<M;i++)
            tmp[i] += (sub[isub].beta[i]-Xb[i])*(sub[isub].beta[i]-Xb[i]);
    }
    for (int i=0;i<M;i++)
        pop->re_prec[i] = rgamma(ALPHA + 0.5*(double)pop->N_SUBS,BETA + 0.5*tmp[i],seed);
    
    free(Xb);
    free(tmp);

}*/

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

void draw_beta_eta(POP *pop,SUB *sub,REP *rep,unsigned long *seed,int iter) {
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
   
  
    /* Add X and -J(V,delta)*/
    
    for (i=0;i<rep->dim_V[0];i++)
        for (j=0;j<ncol;j++)
            V[i*ncol+j] = X[i*ncol+j] - J[i*ncol+j];

            
    calApVinvX(M,V,rep->residuals5,rep->d_Y,rep->dim_V[0],ncol);
   
    /* calculate  (X - J)'(X - J) */
    
    calApVinvA(VpV,V,rep->d_Y,(const int)rep->dim_V[0],(const int)ncol);
   
    for (i=rep->dim_X[1];i<ncol;i++)
        VpV[i*ncol+i] += rep->preceta;
    free(V);
    V = (double *)calloc(ncol*ncol,sizeof(double));
    for (i=0;i<ncol*ncol;i++)
        V[i] = VpV[i];
        
    if (pop->non_parm || pop->GRP) {
        int Nrow = rep->dim_X[1];
        int Ncol = pop->Ncov*Nrow;
        double *Xb = (double *)calloc(Nrow,sizeof(double));
        calAx(Xb,sub->X,pop->beta,(const int)Nrow,(const int)Ncol);
 /*       for (int i=0;i<Nrow;i++)
            M[i] += pop->re_prec[i]*Xb[i];
        for (int i=0;i<Nrow;i++)
            VpV[i*ncol + i] += pop->re_prec[i];*/
         for (int is=0;is<pop->Ns;is++) {
            for (int i=0;i<pop->Nb;i++) {
                    M[i+is*pop->Nb] += pop->re_prec[is]*Xb[i+is*pop->Nb];
                    VpV[(i+is*pop->Nb)*ncol+i+is*pop->Nb] += pop->re_prec[is];
            }
        }
        free(Xb);
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
        fprintf(flog,"iter = %d\n",iter);
        for (i=0;i<pop->Ns;i++)
            fprintf(flog,"re_prec = %lf\n",pop->re_prec[i]);
        fprintf(flog,"rep->preceta = %lf\n",rep->preceta);
        fprintf(flog,"\n");
        for (i=0;i<ncol;i++) {
            for (j=0;j<ncol;j++)
                fprintf(flog,"%lf ",V[i*ncol+j]);
            fprintf(flog,"\n");
        }
        fprintf(flog,"\n");
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

    free(J);
    free(X);
    free(VpV);
    free(M);
    free(mean);
    free(V);
  }


void draw_pop_beta(POP *pop,SUB *sub,unsigned long *seed) {
    double *P,*M,*tmpM,*tmpP,*reprec;
    int N = pop->Ncov*pop->Ns*pop->Nb;
    int SB = pop->Ns*pop->Nb;
    M = (double *)calloc(N,sizeof(double));
    tmpM = (double *)calloc(N,sizeof(double));
    P = (double *)calloc(N*N,sizeof(double));
    tmpP = (double *)calloc(N*N,sizeof(double));
    reprec = (double *)calloc(SB,sizeof(double));
    for (int i=0;i<pop->Ns;i++)
        for (int j=0;j<pop->Nb;j++)
            reprec[i*pop->Nb+j] = pop->re_prec[i];

    for (int isub=0;isub<pop->N_SUBS;isub++) {
        calApVinvX(tmpM,sub[isub].X,sub[isub].beta,reprec,(const int)SB,(const int)N);
        for (int i=0;i<N;i++)
            M[i] += tmpM[i];
    }

    for (int isub=0;isub<pop->N_SUBS;isub++) {
        calApVinvA(tmpP,sub[isub].X,reprec,(const int)SB,(const int)N);
        for (int i=0;i<N;i++)
            for (int j=0;j<N;j++)
                P[i*N+j] += tmpP[i*N+j];
    }
    
/*    for (int is=0;is<pop->Ns;is++) {
        for (int isub=0;isub<pop->N_SUBS;isub++) {
            for (int j=0;j<pop->Nb;j++) {
                P[(j+is*pop->Nb)*N + (j+is*pop->Nb)] += pop->re_prec[j+is*pop->Nb];
            }
        }
    }*/

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
    free(reprec);
    free(P);
    free(M);
    free(tmpM);
    free(tmpP);
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
