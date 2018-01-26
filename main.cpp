#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "cholesky.h"
#include "randgen.h"
#include "fNIRS.h"

double M;
int PREC = 64;
sDLM *dlmStruc;
int maxP = 80;
int main (int argc, const char * argv[]) {
// command line arguments
//  1 - setup file name

    int low_drift,PPP;
    int *dim_design,*dim_HRF;
    double **HRF;
    int i,j,k,N,ell;
	unsigned long *seed;
    double yin;
	char *S,*c;
	FILE *fseed,*fout,*fsetup;
    POP *pop;
    
    double *DCT_basis(int N,double T,int period,int *K);
    double *create_bspline_basis(int N,double TIME,double freq,int int_knots,int *dim_HRF,double *knots,int type);
    double **canonical_HRF(int T,double freq,int *dim_HRF,int type);
    double *convolve(double **design,double **hrf,int *dim_design,int *dim_hrf);
    void setup_fft(REP *rep,const int P);
    void compute_dist_mat(double *dist,int dim,int *m,int *n,double freq);
    void itoa(int n,char s[]);
    
    void mcmc(POP *pop,unsigned long *seed);
    void proj(double *out,double *v,double *u,int length);
   
	M = exp(-PREC*log(2.0));  /* PREC should be set to the compiler precisions 32 or 64 bit */

    seed = (unsigned long *)calloc(3,sizeof(unsigned long));
	fseed = fopen("seed.dat","r");
	for (i=0;i<3;i++)
		int ifs = fscanf(fseed,"%lu ",&(seed[i]));
	fclose(fseed);

    dim_HRF = (int *)calloc(2,sizeof(int));
    dim_design = (int *)calloc(2,sizeof(int));

    pop = (POP *)calloc(1,sizeof(POP));
    pop->expo = 1.99;
    pop->ED = 0;
    pop->priorprec_beta = 0;
    
    S = (char *)malloc(600*sizeof(char));
    char *C;
    int isub=0;
    int irep=0;
    int ifreq=0;
    int istim=0;
    int is=0;
    int idata=0;
    int ides=0;
   
    fsetup = fopen(argv[1],"r");

    while (fgets(S,600,fsetup)) {
        C = strtok(S," ");
        if (!strcmp(C,"NSUBS")) {printf("NSUBS\n");fflush(stdout);
            while ((C = strtok(NULL," ")) != NULL) {
                pop->N_SUBS = atoi(C);
            }
            pop->sub = (SUB *)calloc(pop->N_SUBS,sizeof(SUB));
        }
        else if (!strcmp(C,"GROUP_Analysis")) {printf("GROUP_Analysis\n");fflush(stdout);
             while ((C = strtok(NULL," ")) != NULL) {
                pop->GRP = atoi(C);
            }
        }
        else if (!strcmp(C,"SUB_Replicates")) {printf("SUB_Replicates\n");fflush(stdout);
            while ((C = strtok(NULL," ")) != NULL) {
                printf("C = %s\n",C);
                pop->sub[isub].N_REPS = atoi(C);
            }
            pop->sub[isub].rep = (REP *)calloc(pop->sub[isub].N_REPS,sizeof(REP));
            printf("%d %d\n",isub,pop->sub[isub].N_REPS);
            isub++;
        }
        else if (!strcmp(C,"AR_Degree")) {printf("AR_Degree\n");fflush(stdout);
            while ((C = strtok(NULL," ")) != NULL) {
                PPP = atoi(C);
            }
        }
        else if (!strcmp(C,"SUB_Freq")) {printf("SUB_Freq\n");fflush(stdout);
            while ((C = strtok(NULL," ")) != NULL) {
                pop->sub[ifreq].freq = atof(C);
            }
            ifreq++;
        }
        else if (!strcmp(C,"POP_Stim")) {printf("POP_Stim\n");fflush(stdout);
            while ((C = strtok(NULL," ")) != NULL) {
                pop->Ns = atoi(C);
            }
         }
        else if (!strcmp(C,"SUB_EV")) {printf("SUB_EV\n");fflush(stdout);
            int start = 0;
            int ie = 0;
            int ir = 0;
            while ((C = strtok(NULL," ")) != NULL) {
                if (start) {printf("is = %d, ir = %d, ie = %d\n",is,ir,ie);fflush(stdout);
                    pop->sub[is].rep[ir].Ne[ie] = atoi(C);
                    ie++;
                    if (ie == pop->Ns) {
                        ie = 0;
                        ir++;
                    }
                }
                if (!strcmp(C,"=")) {
                    for (i=0;i<pop->sub[is].N_REPS;i++)
                        pop->sub[is].rep[i].Ne = (int *)calloc(pop->Ns,sizeof(int));
                    start = 1;
                }
            }printf("ir = %d\n",ir);fflush(stdout);
            is++;
        }
        else if (!strcmp(C,"SUB_Data")) {printf("Sub_Data\n");fflush(stdout);
            int start = 0;
            int ir = 0;
            while ((C = strtok(NULL," ")) != NULL) {
                if (start) {
                    pop->sub[idata].rep[ir].dataname = (char *)calloc(100,sizeof(char));
                    strcpy(pop->sub[idata].rep[ir].dataname,C);
                    char *CC;
                    if ((CC = strstr(C,"\n")) != NULL) {
                        CC = strtok(C,"\n");
                        printf("CC = %s %d %d\n",CC,idata,ir);fflush(stdout);
                        strcpy(pop->sub[idata].rep[ir].dataname,CC);
                    }
                    ir++;
                }
                if (!strcmp(C,"="))
                    start = 1;
            }
            idata++;
        }
        else if (!strcmp(C,"SUB_Design")) {printf("SUB_Design\n");fflush(stdout);
            int start = 0;
            int ir = 0;
            while ((C = strtok(NULL," ")) != NULL) {
                if (start) {
                    pop->sub[ides].rep[ir].designname = (char *)calloc(100,sizeof(char));
                    strcpy(pop->sub[ides].rep[ir].designname,C);
                    char *CC;
                    if ((CC = strstr(C,"\n")) != NULL) {
                        CC = strtok(C,"\n");
                        strcpy(pop->sub[ides].rep[ir].designname,CC);
                    }
                    ir++;
                }
                if (!strcmp(C,"="))
                    start = 1;
            }
            ides++;
        }
        else if (!strcmp(C,"Canonical_HRF")) {printf("Canonical_HRF\n");fflush(stdout);
            while ((C = strtok(NULL," ")) != NULL)
                pop->non_parm = atoi(C);
        }
        else if (!strcmp(C,"Num_Basis")) {printf("Num_Basis\n");fflush(stdout);
            while ((C = strtok(NULL," ")) != NULL)
                pop->Nb = atoi(C);
        }
    }
    if (pop->N_SUBS == 1)
        pop->GRP = 0;

    pop->beta = (double *)calloc(pop->Nb*pop->Ns,sizeof(double));
    pop->re_prec = (double *)calloc(pop->Nb*pop->Ns,sizeof(double));
    for (i=0;i<pop->Ns;i++)
        for (j=0;j<pop->Nb;j++)
            pop->re_prec[j+i*pop->Nb] = 1;
    
    printf("N_SUBS = %d\n",pop->N_SUBS);
    printf("Ns = %d\n",pop->Ns);
    for (i=0;i<pop->N_SUBS;i++) {
        printf("Reps SUB %d = %d\n",i,pop->sub[i].N_REPS);
        printf("freq SUB %d = %lf\n",i,pop->sub[i].freq);
    }
    for (i=0;i<pop->N_SUBS;i++)
        for (j=0;j<pop->sub[i].N_REPS;j++)
            for (k=0;k<pop->Ns;k++) {
                printf("pop.sub[%d].rep[%d].Ne[%d] = %d\n",i,j,k,pop->sub[i].rep[j].Ne[k]);
            }
    for (i=0;i<pop->N_SUBS;i++)
        for (j=0;j<pop->sub[i].N_REPS;j++)
            printf("pop.sub[%d].rep[%d].dataname = %s\n",i,j,pop->sub[i].rep[j].dataname);
    for (i=0;i<pop->N_SUBS;i++)
        for (j=0;j<pop->sub[i].N_REPS;j++)
            printf("pop.sub[%d].rep[%d].designname = %s\n",i,j,pop->sub[i].rep[j].designname);
     printf("P = %d\n",PPP);
    printf("non_parm = %d\n",pop->non_parm);
    printf("Nb = %d\n",pop->Nb);
    fclose(fsetup);

/*** OPEN TIME SERIES FILE AND READ DATA ***/
    int int_in;
    for (isub=0;isub<pop->N_SUBS;isub++) {
        for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
            pop->sub[isub].rep[irep].P = PPP;
            pop->sub[isub].rep[irep].df_delta1 = 0.975;
            pop->sub[isub].rep[irep].df_delta2 = 0.975;

            printf("%s\n",pop->sub[isub].rep[irep].dataname);fflush(stdout);
            fout = fopen(pop->sub[isub].rep[irep].dataname,"r");
             N = 0;
            while (fscanf(fout,"%lf ",&(yin)) != EOF)
                N++;
           
           pop->sub[isub].rep[irep].Y = (double *)calloc(N,sizeof(double));
            rewind(fout);
            
            for (i=0;i<N;i++)
                int ifs = fscanf(fout,"%lf ",&(pop->sub[isub].rep[irep].Y[i]));
            
            fclose(fout);
            
    //        for (i=0;i<N;i++)
    //            pop->sub[isub].rep[irep].Y[i] -= 3*sin(2.*3.14159*(double)i/(N-1));
            // Rescale TS
 
            double sdY = 0;
            double meanY = 0;
            for (i=0;i<N;i++)
                meanY += pop->sub[isub].rep[irep].Y[i];
            meanY /= N;
            for (i=0;i<N;i++)
                sdY += (pop->sub[isub].rep[irep].Y[i]-meanY)*(pop->sub[isub].rep[irep].Y[i]-meanY);
            sdY = sqrt(sdY/(N-1));
            for (i=0;i<N;i++)
                pop->sub[isub].rep[irep].Y[i] = (pop->sub[isub].rep[irep].Y[i] - meanY)/sdY;

            /********************************************/
            /*** CONVOLVE DESIGN MATRIX WITH HRF ***/
            switch (pop->non_parm) {
                case 0:
                    HRF = canonical_HRF(35,pop->sub[isub].freq,dim_HRF,1);
                    pop->Nb = dim_HRF[1];
                    break;
                case 1: default:
                    pop->sub[isub].rep[irep].nKnots = pop->Nb-2;
//                    pop->sub[isub].rep[irep].nKnots = pop->Nb;
//                    pop->sub[isub].rep[irep].nKnots = pop->Nb+4;
                    pop->sub[isub].rep[irep].knots = (double *)calloc(pop->sub[isub].rep[irep].nKnots+6,sizeof(double));
                   // HRF = create_bspline_basis(30,pop->sub[isub].freq,pop->Nb-2,dim_HRF,0);
                    // HRF = create_bspline_basis(30,pop->sub[isub].freq,pop->Nb,dim_HRF,1);
                    // HRF = create_bspline_basis(30,pop->sub[isub].freq,pop->Nb+2,dim_HRF,2);
 //                   HRF = create_bspline_basis(35*pop->sub[isub].freq+1,35,pop->sub[isub].freq,pop->sub[isub].rep[irep].nKnots,dim_HRF,pop->sub[isub].rep[irep].knots,0);
 //                   HRF = create_bspline_basis(30*pop->sub[isub].freq,30,pop->sub[isub].freq,pop->sub[isub].rep[irep].nKnots,dim_HRF,pop->sub[isub].rep[irep].knots,1);
 //                     HRF = create_bspline_basis(35*pop->sub[isub].freq+1,35,pop->sub[isub].freq,pop->sub[isub].rep[irep].nKnots,dim_HRF,pop->sub[isub].rep[irep].knots,3);
                    free(pop->sub[isub].rep[irep].knots);
                    break;
            }
 
            fout = fopen("HRF.dat","w");
            for (i=0;i<dim_HRF[0];i++) {
                for (j=0;j<dim_HRF[1];j++)
                    fprintf(fout,"%lf ",HRF[i][j]);
                fprintf(fout,"\n");
            }
            fclose(fout);
            dim_design[0] = N;
            
            fout = fopen(pop->sub[isub].rep[irep].designname,"r");
            printf("%s\n",pop->sub[isub].rep[irep].designname);fflush(stdout);
            int cnt=0;
            while (fscanf(fout,"%lf ",&yin) != EOF)
                cnt++;
            
            rewind(fout);
            
            dim_design[1] = cnt/(N);
        
            double **design2;
            design2 = (double **)calloc(dim_design[0],sizeof(double *));
            for (i=0;i<dim_design[0];i++)
                design2[i] = (double *)calloc(dim_design[1],sizeof(double));
            for (i=0;i<dim_design[0];i++)
                for (j=0;j<dim_design[1];j++)
                    int ifs = fscanf(fout,"%lf ",&(design2[i][j]));
        
            fclose(fout);

            pop->sub[isub].rep[irep].X = convolve(design2,HRF,dim_design,dim_HRF);
            
/*            if (pop->non_parm == 10) { //  orthonormalize bases
                double tmp;
                double *out = (double *)calloc(dim_design[0],sizeof(double));
                double **v = (double **)calloc(pop->Nb,sizeof(double *));
                for (i=0;i<pop->Nb;i++)
                    v[i] = (double *)calloc(dim_design[0],sizeof(double));
                
                for (ell=0;ell<pop->Ns;ell++) {
                    for (i=0;i<dim_design[0];i++) {
                        int idx=0;
                        for (j=ell*pop->Nb;j<(ell+1)*pop->Nb;j++) {
                            v[idx][i] = pop->sub[isub].rep[irep].X[i*j];
                            idx++;
                        }
                    }
                    
                    for (i=0;i<pop->Nb;i++) {  // modified (or stable) Gram-Schmidt orthonormalization
                        tmp = 0;
                        for (k=0;k<dim_design[0];k++)
                            tmp += v[i][k]*v[i][k];
                        tmp = sqrt(tmp);
                        for (k=0;k<dim_design[0];k++)
                            v[i][k] /= tmp;
 
                        for (j=(i+1);j<pop->Nb;j++) {
                            proj(out,v[j],v[i],dim_design[0]);
                            for (k=0;k<dim_design[0];k++)
                                v[j][k] -= out[k];
                           
                        }
                    }
                    
                    for (i=0;i<dim_design[0];i++) {
                        int idx=0;
                        for (j=ell*pop->Nb;j<(ell+1)*pop->Nb;j++) {
                            pop->sub[isub].rep[irep].X[i][j] = v[idx][i];
                            idx++;
                        }
                    }
                }
                
                free(out);
                for (i=0;i<pop->Nb;i++)
                    free(v[i]);
                free(v);
            }*/
            
            for (i=0;i<dim_design[0];i++)
                free(design2[i]);
            free(design2);
            
            pop->sub[isub].rep[irep].dim_X = (int *)calloc(2,sizeof(int));
            pop->sub[isub].rep[irep].dim_X[0] = dim_design[0];
            pop->sub[isub].rep[irep].dim_X[1] = dim_design[1]*dim_HRF[1];

            
            double max = 0;
            double tmp = 0;
            if (pop->non_parm > -1) {
                for (j=0;j<dim_design[1]*dim_HRF[1];j++) { // normalize
                    tmp = 0;
                    for (i=0;i<dim_design[0];i++)
                       tmp = (tmp > pop->sub[isub].rep[irep].X[i*dim_design[1]*dim_HRF[1]+j]) ? tmp:pop->sub[isub].rep[irep].X[i*dim_design[1]*dim_HRF[1]+j];
 //                       tmp += pop->sub[isub].rep[irep].X[i][j]*pop->sub[isub].rep[irep].X[i][j];
 //                   tmp = sqrt(tmp);
                    max = (tmp > max) ? tmp:max;
                }
                for (j=0;j<dim_design[1]*dim_HRF[1];j++) { // normalize
                    for (i=0;i<dim_design[0];i++)
                        pop->sub[isub].rep[irep].X[i*dim_design[1]*dim_HRF[1]+j] /= max;
                }
            }
            printf("max = %lf\n",max);fflush(stdout);
            C = (char *)malloc(3*sizeof(char));
            S = strcpy(S,"X");
            itoa(isub,C);
            S = strcat(S,C);
            itoa(irep,C);
            S = strcat(S,C);
            free(C);
            S = strcat(S,".dat");
            //printf("S: %s\n",S);
            fout = fopen(S,"w");
            for (i=0;i<pop->sub[isub].rep[irep].dim_X[0];i++) {
                for (j=0;j<pop->sub[isub].rep[irep].dim_X[1];j++) {
                    fprintf(fout,"%.20lf ",pop->sub[isub].rep[irep].X[i*dim_design[1]*dim_HRF[1]+j]);
                }
                fprintf(fout,"\n");
            }
            fclose(fout);
      
            for (i=0;i<dim_HRF[0];i++)
                free(HRF[i]);
            free(HRF);
           
  
            /*** CREATE DISCRETE COSINE BASIS TO REMOVE LOW FREQUENCY DRIFT ***/
            
            double T = (double)(N-1)/(double)pop->sub[0].freq;
            int K=0;
            
            pop->sub[isub].rep[irep].dim_V = (int *)calloc(2,sizeof(int));
            
            pop->sub[isub].rep[irep].nKnots = 15;
            pop->sub[isub].rep[irep].knots = (double *)calloc(pop->sub[isub].rep[irep].nKnots+6,sizeof(double));
            
            low_drift=1;
              switch (low_drift) {
                case 0:
                    pop->sub[isub].rep[irep].V = DCT_basis(N,T,80,&K);
                    pop->sub[isub].rep[irep].dim_V[0] = N;
                    pop->sub[isub].rep[irep].dim_V[1] = K+1;
                    break;
                case 1: default:
                    pop->sub[isub].rep[irep].V = create_bspline_basis(N,T,pop->sub[isub].freq,pop->sub[isub].rep[irep].nKnots,pop->sub[isub].rep[irep].dim_V,pop->sub[isub].rep[irep].knots,0);
                    pop->sub[isub].rep[irep].nKnots += 6;                  
                    break;
            }


            printf("NKNOTS = %d, dim_V[1] = %d\n",pop->sub[isub].rep[irep].nKnots,pop->sub[isub].rep[irep].dim_V[1]);
            for (i=0;i<pop->sub[isub].rep[irep].nKnots;i++)
                printf("%lf ",pop->sub[isub].rep[irep].knots[i]);
                printf("\n");
             if (low_drift) {
 /*               double tmp;
                double *out = (double *)calloc(pop->sub[isub].rep[irep].dim_V[0],sizeof(double));
                double *v = (double *)calloc(pop->sub[isub].rep[irep].dim_V[0],sizeof(double));
                double *u = (double *)calloc(pop->sub[isub].rep[irep].dim_V[0],sizeof(double));
                for (int j=1;j<pop->sub[isub].rep[irep].dim_V[1];j++) {  // orthogonalize using Gram-Schmidt
                    for (int i=0;i<pop->sub[isub].rep[irep].dim_V[0];i++)
                        v[i] = pop->sub[isub].rep[irep].V[i][j];
                    for (int k=0;k<j;k++) {
                        for (int i=0;i<pop->sub[isub].rep[irep].dim_V[0];i++)
                            u[i] = pop->sub[isub].rep[irep].V[i][k];
                        proj(out,v,u,pop->sub[isub].rep[irep].dim_V[0]);
                        for (int i=0;i<pop->sub[isub].rep[irep].dim_V[0];i++)
                            pop->sub[isub].rep[irep].V[i][j] -= out[i];
                    }
                }
                
                free(out);
                free(v);
                free(u);*/
                
        /*        for (int j=0;j<pop->sub[isub].rep[irep].dim_V[1];j++) { // normalize
                    tmp = 0;
                    max = 0;
                    for (int i=0;i<pop->sub[isub].rep[irep].dim_V[0];i++)
                        tmp += pop->sub[isub].rep[irep].V[i][j]*pop->sub[isub].rep[irep].V[i][j];
                    tmp = sqrt(tmp);
                    max = (tmp > max) ? tmp:max;
                }
                for (int j=0;j<pop->sub[isub].rep[irep].dim_V[1];j++) {
                    for (int i=0;i<pop->sub[isub].rep[irep].dim_V[0];i++)
                        pop->sub[isub].rep[irep].V[i][j] /= max;
                }*/
            }
            C = (char *)malloc(3*sizeof(char));
            S = strcpy(S,"V");
            itoa(isub,C);
            S = strcat(S,C);
            itoa(irep,C);
            S = strcat(S,C);
            free(C);
            S = strcat(S,".dat");
            //printf("S: %s\n",S);
            fout = fopen(S,"w");
   
            for (i=0;i<pop->sub[isub].rep[irep].dim_V[0];i++) {
                for (j=0;j<pop->sub[isub].rep[irep].dim_V[1];j++)
                    fprintf(fout,"%lf ",pop->sub[isub].rep[irep].V[i*pop->sub[isub].rep[irep].dim_V[1]+j]);
                fprintf(fout,"\n");
            }
            fclose(fout);
            
            
            /*** AR TERMS ***/
            
  //          N = N - pop->P;
            
            pop->sub[isub].rep[irep].dim_W = (int *)calloc(2,sizeof(int));
            pop->sub[isub].rep[irep].dim_W[0] = N;
            pop->sub[isub].rep[irep].dim_W[1] = PPP;
            
//            pop->sub[isub].rep[irep].W = (double *)calloc(N*PPP,sizeof(double));
            pop->sub[isub].rep[irep].W = (double *)calloc(N*maxP,sizeof(double)); // needed to select DLM hyperpriors
            
            for (i=PPP;i<N;i++) {
                for (j=0;j<PPP;j++) {
                    pop->sub[isub].rep[irep].W[i*PPP + PPP-j-1] = pop->sub[isub].rep[irep].Y[i+j-PPP];
//                    pop->sub[isub].rep[irep].W[i][PPP-j-1] = pop->sub[isub].rep[irep].Y[i+j-PPP];
                }
            }

            fout = fopen("AR_matrix.dat","w");
            for (i=0;i<N;i++) {
                for (j=0;j<PPP;j++)
                    fprintf(fout,"%lf ",pop->sub[isub].rep[irep].W[i*PPP+j]);
                fprintf(fout,"\n");
            }
            fclose(fout);
            
            C = (char *)malloc(3*sizeof(char));
            S = strcpy(S,"W");
            itoa(isub,C);
            S = strcat(S,C);
            itoa(irep,C);
            S = strcat(S,C);
            free(C);
            S = strcat(S,".dat");
            //printf("S: %s\n",S);
            fout = fopen(S,"w");
    
            for (i=0;i<pop->sub[isub].rep[irep].dim_W[0];i++) {
                for (j=0;j<pop->sub[isub].rep[irep].dim_W[1];j++)
                    fprintf(fout,"%lf ",pop->sub[isub].rep[irep].W[i*PPP+j]);
                fprintf(fout,"\n");
            }
            fclose(fout);
            
            
            pop->sub[isub].rep[irep].XpX = (double *)calloc(pop->sub[isub].rep[irep].dim_X[1]*pop->sub[isub].rep[irep].dim_X[1],sizeof(double));
//            for (i=0;i<pop->sub[isub].rep[irep].dim_X[1];i++)
//                pop->sub[isub].rep[irep].XpX[i] = (double *)calloc(pop->sub[isub].rep[irep].dim_X[1],sizeof(double));
            
            for (k=0;k<pop->sub[isub].rep[irep].dim_X[0];k++)
                for (i=0;i<pop->sub[isub].rep[irep].dim_X[1];i++)
                    for (j=0;j<pop->sub[isub].rep[irep].dim_X[1];j++)
                        pop->sub[isub].rep[irep].XpX[i*pop->sub[isub].rep[irep].dim_X[1]+j] += pop->sub[isub].rep[irep].X[k*pop->sub[isub].rep[irep].dim_X[1]+i]*pop->sub[isub].rep[irep].X[k*pop->sub[isub].rep[irep].dim_X[1]+j];
             
            pop->sub[isub].rep[irep].prbs = (double *)calloc(pop->sub[isub].rep[irep].dim_V[0],sizeof(double));
            
            meanY = sdY = 0;
            for (i=1;i<pop->sub[isub].rep[irep].dim_V[0];i++) {
//                pop->sub[isub].rep[irep].prbs[i] = pop->sub[isub].rep[irep].Y[i] - pop->sub[isub].rep[irep].Y[i-1];
//                pop->sub[isub].rep[irep].prbs[i] = (fabs(pop->sub[isub].rep[irep].Y[i] - pop->sub[isub].rep[irep].Y[i-1]));
//                pop->sub[isub].rep[irep].prbs[i] *= pop->sub[isub].rep[irep].prbs[i];
                pop->sub[isub].rep[irep].prbs[i]  = 1;
 //               meanY += pop->sub[isub].rep[irep].prbs[i];
            }
 /*           meanY /= pop->sub[isub].rep[irep].dim_V[0]-1;
            for (i=1;i<pop->sub[isub].rep[irep].dim_V[0];i++)
                sdY += (pop->sub[isub].rep[irep].prbs[i]-meanY)*(pop->sub[isub].rep[irep].prbs[i]-meanY);
            sdY = sqrt(sdY/(pop->sub[isub].rep[irep].dim_V[0]-1));
            for (i=1;i<pop->sub[isub].rep[irep].dim_V[0];i++)
                pop->sub[isub].rep[irep].prbs[i] = (pop->sub[isub].rep[irep].prbs[i] - meanY)/sdY;
            int cntg=0;
            for (i=0;i<pop->sub[isub].rep[irep].dim_V[0];i++) {
                if (fabs(pop->sub[isub].rep[irep].prbs[i]) > 2.5)
                    cntg++;
            }
            int *large;
            large = (int *)calloc(cntg,sizeof(int));
            int cntr=0;
            for (i=0;i<pop->sub[isub].rep[irep].dim_V[0];i++) {
                if (fabs(pop->sub[isub].rep[irep].prbs[i]) > 2.5) {
                    large[cntr] = i;
                    cntr++;
                }
            }
             for (i=1;i<pop->sub[isub].rep[irep].dim_V[0];i++)
                pop->sub[isub].rep[irep].prbs[i] = 0.2/(pop->sub[isub].rep[irep].dim_V[0]-1 - cntg);
            for (i=0;i<cntg;i++)
                pop->sub[isub].rep[irep].prbs[large[i]] = 0.8/(cntg);
     */
            for (i=1;i<pop->sub[isub].rep[irep].dim_V[0];i++)
                pop->sub[isub].rep[irep].prbs[i] += pop->sub[isub].rep[irep].prbs[i-1];
            for (i=1;i<pop->sub[isub].rep[irep].dim_V[0];i++)
                pop->sub[isub].rep[irep].prbs[i] /= pop->sub[isub].rep[irep].prbs[pop->sub[isub].rep[irep].dim_V[0]-1];
 /*           for (i=0;i<pop->sub[isub].rep[irep].dim_V[0];i++)
                printf("%lf\n",pop->sub[isub].rep[irep].prbs[i]);*/
            fout = fopen("prbs.dat","w");
            for (i=0;i<pop->sub[isub].rep[irep].dim_V[0];i++)
                fprintf(fout,"%lf\n",pop->sub[isub].rep[irep].prbs[i]);
            fclose(fout); 

            pop->sub[isub].rep[irep].residuals = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0],sizeof(double));
            pop->sub[isub].rep[irep].residuals1 = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0],sizeof(double));
            pop->sub[isub].rep[irep].residuals2 = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0],sizeof(double));
            pop->sub[isub].rep[irep].residuals3 = (double *)calloc(pop->sub[isub].rep[irep].dim_W[0],sizeof(double));
            pop->sub[isub].rep[irep].residuals4 = (double *)calloc(pop->sub[isub].rep[irep].dim_V[0],sizeof(double));
            pop->sub[isub].rep[irep].residuals5 = (double *)calloc(pop->sub[isub].rep[irep].dim_V[0],sizeof(double));
            pop->sub[isub].rep[irep].eta = (double *)calloc(pop->sub[isub].rep[irep].dim_V[1],sizeof(double));
            pop->sub[isub].rep[irep].Veta = (double *)calloc(pop->sub[isub].rep[irep].dim_V[0],sizeof(double));
            pop->sub[isub].rep[irep].Jeta = (double *)calloc(pop->sub[isub].rep[irep].dim_V[0],sizeof(double));
            pop->sub[isub].rep[irep].J = (double *)calloc(pop->sub[isub].rep[irep].dim_V[0]*pop->sub[isub].rep[irep].dim_V[1],sizeof(double));
 //           for (i=0;i<pop->sub[isub].rep[irep].dim_V[0];i++)
 //               pop->sub[isub].rep[irep].J[i] = (double *)calloc(pop->sub[isub].rep[irep].dim_V[1],sizeof(double));
            pop->sub[isub].rep[irep].mVeta = (double *)calloc(pop->sub[isub].rep[irep].dim_V[0],sizeof(double));
 
            pop->sub[isub].rep[irep].delta = (double *)calloc(pop->sub[isub].rep[irep].dim_W[0]*maxP,sizeof(double));
  //          for (i=0;i<pop->sub[isub].rep[irep].dim_W[0];i++)
  //              pop->sub[isub].rep[irep].delta[i] = (double *)calloc(pop->sub[isub].rep[irep].dim_W[1],sizeof(double));

           pop->sub[isub].rep[irep].mdelta = (double *)calloc(pop->sub[isub].rep[irep].dim_W[0]*maxP,sizeof(double));
  //          for (i=0;i<pop->sub[isub].rep[irep].dim_W[0];i++)
  //              pop->sub[isub].rep[irep].mdelta[i] = (double *)calloc(pop->sub[isub].rep[irep].dim_W[1],sizeof(double));
           pop->sub[isub].rep[irep].mdelta2 = (double *)calloc(pop->sub[isub].rep[irep].dim_W[0]*maxP,sizeof(double));
  //          for (i=0;i<pop->sub[isub].rep[irep].dim_W[0];i++)
  //              pop->sub[isub].rep[irep].mdelta2[i] = (double *)calloc(pop->sub[isub].rep[irep].dim_W[1],sizeof(double));

/*            pop->sub[isub].rep[irep].delta[pop->P-1][0] = 0.6;
            pop->sub[isub].rep[irep].delta[pop->P-1][1] = 0.2;
            pop->sub[isub].rep[irep].delta[pop->P-1][2] = 0.15;
            pop->sub[isub].rep[irep].delta[pop->P-1][3] = 0.10;*/
 //           pop->sub[isub].rep[irep].delta[pop->P-1][4] = -0.02;
 //           pop->sub[isub].rep[irep].delta[pop->P-1][5] = -0.04;
 //           pop->sub[isub].rep[irep].delta[pop->P-1][6] = 0.05;
 //           pop->sub[isub].rep[irep].delta[pop->P-1][7] = -0.07;
            pop->sub[isub].rep[irep].Wdelta = (double *)calloc(pop->sub[isub].rep[irep].dim_W[0],sizeof(double));
            pop->sub[isub].rep[irep].mWdelta = (double *)calloc(pop->sub[isub].rep[irep].dim_W[0],sizeof(double));

            pop->sub[isub].rep[irep].Xbeta = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0],sizeof(double));
            pop->sub[isub].rep[irep].Hbeta = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0],sizeof(double));
            pop->sub[isub].rep[irep].H = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0]*pop->sub[isub].rep[irep].dim_X[1],sizeof(double));
//            for (i=0;i<pop->sub[isub].rep[irep].dim_X[0];i++)
//                pop->sub[isub].rep[irep].H[i] = (double *)calloc(pop->sub[isub].rep[irep].dim_X[1],sizeof(double));
            
            pop->sub[isub].rep[irep].precY = (double *)calloc(3,sizeof(double));
            pop->sub[isub].rep[irep].precY[0] = 1;
            pop->sub[isub].rep[irep].precY[1] = 1./16.;
            pop->sub[isub].rep[irep].precY[2] = 1./64.;
            pop->sub[isub].rep[irep].d_Y = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0],sizeof(double));
            pop->sub[isub].rep[irep].md_Y = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0],sizeof(double));
            pop->sub[isub].rep[irep].sd_Y = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0],sizeof(double));
            for (i=PPP;i<pop->sub[isub].rep[irep].dim_X[0];i++)
                pop->sub[isub].rep[irep].d_Y[i] = 1;
             
            pop->sub[isub].rep[irep].K = 1;
            pop->sub[isub].rep[irep].preceta = 1;
            pop->sub[isub].rep[irep].precdelta = 100.;
            pop->sub[isub].rep[irep].LY = 1./16.;
            pop->sub[isub].rep[irep].precYstart = 1;
            pop->sub[isub].rep[irep].mean_precY = 0;
            pop->sub[isub].rep[irep].mean_res = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0],sizeof(double));
            pop->sub[isub].rep[irep].mean_d_Y = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0],sizeof(double));
            pop->sub[isub].rep[irep].mean_fit = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0],sizeof(double));
            pop->sub[isub].rep[irep].mean_Y = (double *)calloc(pop->sub[isub].rep[irep].dim_V[0],sizeof(double));
           
            pop->sub[isub].rep[irep].prop_sd = (double *)calloc(4,sizeof(double));
            pop->sub[isub].rep[irep].prop_sd[0] = 0.005;
            pop->sub[isub].rep[irep].prop_sd[1] = 0.05;
            pop->sub[isub].rep[irep].prop_sd[2] = 0.001;
            pop->sub[isub].rep[irep].prop_sd[3] = 10;
            pop->sub[isub].rep[irep].maxsteps = 50;
            
            pop->sub[isub].rep[irep].accept = (int *)calloc(4,sizeof(int));
            pop->sub[isub].rep[irep].attempt = (int *)calloc(4,sizeof(int));
            
        }
        
        pop->sub[isub].beta = (double *)calloc(pop->Nb*(pop->Ns),sizeof(double));
  
        // int lpbeta = pop->Ns*pop->Nb;
        int lpbeta = pop->Ns;

 //       pop->sub[isub].precbeta = 1.;

    }
    printf("here\n");fflush(NULL);
    double **D;
    D = (double **)calloc(dim_HRF[1],sizeof(double *));
    for (i=0;i<dim_HRF[1];i++)
        D[i] = (double *)calloc(dim_HRF[1],sizeof(double));
    for (i=0;i<dim_HRF[1];i++) {
        if (i<=1) {
            D[i][i] = 1;
        }
         else {
            D[i][i-2] = -1;
            D[i][i-1] = 2;
        }
    }
/*    for (i=0;i<dim_HRF[1];i++) {
        if (i==0) {
            D[i][0] = -2;
            D[i][i+1] = 1;
        }
        else if (i==dim_HRF[1]-1) {
            D[i][i-1] = 1;
            D[i][i] = -2;
        }
        else {
            D[i][i] = -2;
            D[i][i-1] = D[i][i+1] = 1;
        }
    }*/
    
    pop->Q = (double **)calloc(dim_HRF[1]*pop->Ns,sizeof(double *));
    for (i=0;i<dim_HRF[1]*pop->Ns;i++)
        pop->Q[i] = (double *)calloc(dim_HRF[1]*pop->Ns,sizeof(double));
 
   if (pop->non_parm) {
        for (k=0;k<dim_HRF[1];k++) {
            for (i=0;i<dim_HRF[1];i++) {
                for (j=0;j<dim_HRF[1];j++) {
                    pop->Q[i][j] += D[k][i]*D[k][j];
                }
            }
        }
        for (k=1;k<pop->Ns;k++) {
            for (i=0;i<dim_HRF[1];i++)
                for (j=0;j<dim_HRF[1];j++)
                    pop->Q[i+k*dim_HRF[1]][j+k*dim_HRF[1]] = pop->Q[i][j]; 
        }
    }
    else {
        for (i=0;i<dim_HRF[1]*pop->Ns;i++)
            pop->Q[i][i] = 1;
    }
    for (i=0;i<dim_HRF[1];i++)
        free(D[i]);
    free(D);

/*    for (i=0;i<dim_HRF[1]*pop->Ns;i++) {
        for (j=0;j<dim_HRF[1]*pop->Ns;j++) {
            printf("%2d ",(int)pop->Q[i][j]);
        }
        printf("\n");
    }*/

    free(dim_HRF);
    free(dim_design);

    int dimV = pop->sub[0].rep[0].dim_V[1]+5;
    D = (double **)calloc(dimV,sizeof(double *));
    for (i=0;i<dimV;i++)
        D[i] = (double *)calloc(dimV,sizeof(double));
    for (i=0;i<dimV;i++) {
        if (i<=1) {
            D[i][i] = 1;
        }
         else {
            D[i][i-2] = -1;
            D[i][i-1] = 2;
        }
    }
    
    pop->Q2 = (double **)calloc(dimV,sizeof(double *));
    for (i=0;i<dimV;i++)
        pop->Q2[i] = (double *)calloc(dimV,sizeof(double));
    for (k=0;k<dimV;k++) {
        for (i=0;i<dimV;i++) {
            for (j=0;j<dimV;j++) {
                pop->Q2[i][j] += D[k][i]*D[k][j];
            }
        }
    }
    for (i=0;i<dimV;i++)
        free(D[i]);
    free(D);

    if (low_drift) {
        for (i=0;i<dimV;i++) {
            for (j=0;j<dimV;j++) {
                pop->Q2[i][j] = 0;
            }
        }
        for (i=0;i<dimV;i++)
            pop->Q2[i][i] = 1;
        
    }
 /*   for (i=0;i<dimV;i++) {
        for (j=0;j<dimV;j++) {
            printf("%2d ",(int)pop->Q2[i][j]);
        }
        printf("\n");
    }*/
    
    int max = 0;
    for (i=0;i<pop->N_SUBS;i++)
        for (j=0;j<pop->sub[i].N_REPS;j++)
            max = (max > pop->sub[i].rep[j].dim_X[0]) ? max:pop->sub[i].rep[j].dim_X[0];
 
     dlmStruc = (sDLM *)calloc(max,sizeof(sDLM));
    for (i=0;i<max;i++)
        dlmStruc[i].a = (double *)calloc(maxP,sizeof(double));
    for (i=0;i<max;i++)
        dlmStruc[i].m = (double *)calloc(maxP,sizeof(double));
    for (i=0;i<max;i++) {
        dlmStruc[i].R = (double *)calloc(maxP*maxP,sizeof(double));
//        for (j=0;j<PPP;j++)
//            dlmStruc[i].R[j] = (double *)calloc(PPP,sizeof(double));
    }
    for (i=0;i<max;i++) {
        dlmStruc[i].C = (double *)calloc(maxP*maxP,sizeof(double));
//        for (j=0;j<PPP;j++)
//            dlmStruc[i].C[j] = (double *)calloc(PPP,sizeof(double));
    }
/*** CALL MCMC ***/
    printf("ENTERING MCMC\n");fflush(NULL);
    mcmc(pop,seed);
    printf("EXITING MCMC\n");fflush(NULL);
/*****************/

/*** SAVE NEW SEED FOR NEXT RUN ***/
    
    fseed = fopen("seed.dat","w");
	for (i=0;i<3;i++)
		fprintf(fseed,"%lu ",seed[i]);
	fclose(fseed);
	free(seed);
    return 0;
  
    free(S);
}

void itoa(int n,char s[])
{
    int i,sign;
    void reverse(char s[]);
    
    if ((sign = n) < 0)
        n = -n;
    i = 0;
    do {
        s[i++] = n % 10 + '0';
    } while ((n /= 10) > 0);
    if (sign < 0)
        s[i++] = '-';
    s[i] = '\0';
    reverse(s);
}

void reverse(char s[])
{
    int c,i,j;
    
    for (i=0,j = strlen(s)-1;i<j;i++,j--) {
        c = s[i];
        s[i] = s[j];
        s[j] = c;
    }
}
