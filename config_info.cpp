#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fNIRS.h"

int max_dimX;
extern int maxP;
extern int MAX_ITER;
extern int BURN_IN;
extern FILE *flog,*fseed;
extern sDLM *dlmStruc;

void load_config_info(POP *pop,const char *config_file,unsigned long *seed)
{
    int N;
    int isub=0;
    int irep=0;
    int ifreq=0;
    int istim=0;
    int is=0;
    int idata=0;
    int ides=0;
    double yin;
    char *S,*C,*covarMat;
    SUB *sub;
    REP *rep;
    FILE *fconfig,*fcovar=NULL;

    void center_covars(double *X,const int nrow,const int ncol); 
   
    S = (char *)malloc(1000*sizeof(char));
    fconfig = fopen(config_file,"r");

    while (fgets(S,1000,fconfig)) {
        C = strtok(S," ");
        if (!strcmp(C,"GROUP_Analysis")) {
            while ((C = strtok(NULL," ")) != NULL) {
                pop->GRP = atoi(C);
            }
        }
        else if (!strcmp(C,"Expected_Knots")) {
            while ((C = strtok(NULL," ")) != NULL) {
                pop->knots = atoi(C);
            }
        }
        else if (!strcmp(C,"SEED_Matrix")) {
            int start = 0;
            char *seedMat;
            while ((C = strtok(NULL," ")) != NULL) {
                if (start) {
                    seedMat = (char *)calloc(300,sizeof(char));
                    strcpy(seedMat,C);
                    char *CC;
                    if ((CC = strstr(C,"\n")) != NULL) {
                        CC = strtok(C,"\n");
                        strcpy(seedMat,CC);
                    }
                }
                if (!strcmp(C,"="))
                    start = 1;
            }
    	    fseed = fopen(seedMat,"r+");
	        for (int i=0;i<3;i++)
		        int ifs = fscanf(fseed,"%lu ",&(seed[i]));
	        free(seedMat);
        }
        else if (!strcmp(C,"COVAR_Matrix")) {
            int start = 0;
            while ((C = strtok(NULL," ")) != NULL) {
                if (start) {
                    covarMat = (char *)calloc(300,sizeof(char));
                    strcpy(covarMat,C);
                    char *CC;
                    if ((CC = strstr(C,"\n")) != NULL) {
                        CC = strtok(C,"\n");
                        strcpy(covarMat,CC);
                    }
                }
                if (!strcmp(C,"="))
                    start = 1;
            }
        }
        else if (!strcmp(C,"MAX_ITER")) {
            while ((C = strtok(NULL," ")) != NULL) {
                MAX_ITER = atoi(C);
            }
        }
        else if (!strcmp(C,"BURN_IN")) {
            while ((C = strtok(NULL," ")) != NULL) {
                BURN_IN = atoi(C);
            }
        }
        else if (!strcmp(C,"NSUBS")) {
            while ((C = strtok(NULL," ")) != NULL) {
                pop->N_SUBS = atoi(C);
            }
            pop->sub = (SUB *)calloc(pop->N_SUBS,sizeof(SUB));
        }
        else if (!strcmp(C,"SUB_Replicates")) {
            while ((C = strtok(NULL," ")) != NULL) {
                pop->sub[isub].N_REPS = atoi(C);
            }
            pop->sub[isub].rep = (REP *)calloc(pop->sub[isub].N_REPS,sizeof(REP));
            isub++;
        }
        else if (!strcmp(C,"SUB_Freq")) {
            while ((C = strtok(NULL," ")) != NULL) {
                pop->sub[ifreq].freq = atof(C);
            }
            ifreq++;
        }
        else if (!strcmp(C,"POP_Stim")) {
            while ((C = strtok(NULL," ")) != NULL) {
                pop->Ns = atoi(C);
            }
         }
 /*       else if (!strcmp(C,"SUB_EV")) {
            int start = 0;
            int ie = 0;
            int ir = 0;
            while ((C = strtok(NULL," ")) != NULL) {
                if (start) {
                    pop->sub[is].rep[ir].Ne[ie] = atoi(C);
                    ie++;
                    if (ie == pop->Ns) {
                        ie = 0;
                        ir++;
                    }
                }
                if (!strcmp(C,"=")) {
                    for (int i=0;i<pop->sub[is].N_REPS;i++)
                        pop->sub[is].rep[i].Ne = (int *)calloc(pop->Ns,sizeof(int));
                    start = 1;
                }
            }
            is++;
        }*/
        else if (!strcmp(C,"SUB_Data")) {
            int start = 0;
            int ir = 0;
            while ((C = strtok(NULL," ")) != NULL) {
                if (start) {
                    pop->sub[idata].rep[ir].dataname = (char *)calloc(100,sizeof(char));
                    strcpy(pop->sub[idata].rep[ir].dataname,C);
                    char *CC;
                    if ((CC = strstr(C,"\n")) != NULL) {
                        CC = strtok(C,"\n");
                        strcpy(pop->sub[idata].rep[ir].dataname,CC);
                    }
                    ir++;
                }
                if (!strcmp(C,"="))
                    start = 1;
            }
            idata++;
        }
        else if (!strcmp(C,"SUB_Design")) {
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
        else if (!strcmp(C,"Include_temporal_derivative")) {
            while ((C = strtok(NULL," ")) != NULL) {
                if (atoi(C) == 1)
                    pop->Nb = 2;
                else 
                    pop->Nb = 1;
            }
        }
    }
    pop->non_parm = 0;
    fclose(fconfig);
    free(S);

    if (BURN_IN >= MAX_ITER) {
        printf("Error in configuration file %s\n",config_file);
        printf("\t BURN_IN must be less than MAX_ITER. Exiting program...\n");
        exit(0);   
    }

    if (pop->N_SUBS == 1)
        pop->GRP = 0;
     
 
   if (pop->GRP) {
        fcovar = fopen(covarMat,"r");
        if (fcovar == NULL) {
            printf("A covariance matrix is missing.\n");
            printf("\t When doing a group level analysis, a covariance matrix is required.\n");
            printf("\t At a minimum, the intercept term is required. Exiting program...\n");
            exit(0);
        }
        N = 0;
        while (fscanf(fcovar,"%lf ",&(yin)) != EOF)
            N++;
        if ((N%pop->N_SUBS)) {
            printf("N = %d N_SUBS = %d mod = %d\n",N,pop->N_SUBS,N%pop->N_SUBS);
            printf("Error in %s \n",covarMat);
            printf("\t The covariance matrix must be balanced; that is, no missing data.\n");
            printf("\t Number of elements, N = %d, is not a multiple of the number of subjects, N_SUB = %d. Exiting program...\n",N,pop->N_SUBS);
            exit(0);
        }           
        pop->Ncov = N/pop->N_SUBS;
        rewind(fcovar);
        
        double *covars;
        covars = (double *)calloc(N,sizeof(double));
        for (isub=0;isub<pop->N_SUBS;isub++)
            for (int j=0;j<pop->Ncov;j++)
                int ifs = fscanf(fcovar,"%lf ",&covars[isub*pop->Ncov + j]);
        center_covars(covars,(const int)pop->N_SUBS,(const int)pop->Ncov);
        for (isub=0;isub<pop->N_SUBS;isub++) {
            sub = &(pop->sub[isub]);
            sub->dim_X = (int *)calloc(2,sizeof(int));
            sub->dim_X[0] = pop->Nb*pop->Ns;
            sub->dim_X[1] = pop->Ncov*pop->Nb*pop->Ns;
            sub->X = (double *)calloc(sub->dim_X[0]*sub->dim_X[1],sizeof(double));
            for (int i=0;i<sub->dim_X[0];i++)
                for (int j=0;j<pop->Ncov;j++) {
                    sub->X[(j+i*pop->Ncov)+i*sub->dim_X[1]] = covars[isub*pop->Ncov+j];
                }
            for (int i=0;i<sub->dim_X[0];i++) {
                for (int j=0;j<sub->dim_X[1];j++)
                    fprintf(flog,"%lf ",sub->X[i*sub->dim_X[1]+j]);
                fprintf(flog,"\n");
            }
            fprintf(flog,"\n");
        }
        free(covars);
        fclose(fcovar);
    }
    else {
        pop->Ncov = 1;
    }
     
    pop->beta = (double *)calloc(pop->Ncov*pop->Nb*pop->Ns,sizeof(double));
    pop->re_prec = (double *)calloc(pop->Ns,sizeof(double));
    for (int i=0;i<pop->Ns;i++)
        pop->re_prec[i] = 1;

    for (int i=0;i<pop->N_SUBS;i++) {
        fprintf(flog,"Reps SUB %d = %d\n",i,pop->sub[i].N_REPS);
        fprintf(flog,"freq SUB %d = %lf\n",i,pop->sub[i].freq);
    }
/*    for (int i=0;i<pop->N_SUBS;i++)
        for (int j=0;j<pop->sub[i].N_REPS;j++)
            for (int k=0;k<pop->Ns;k++) {
                fprintf(flog,"pop.sub[%d].rep[%d].Ne[%d] = %d\n",i,j,k,pop->sub[i].rep[j].Ne[k]);
            }*/
    for (int i=0;i<pop->N_SUBS;i++) {
        sub = &(pop->sub[i]);
        for (int j=0;j<sub->N_REPS;j++) {
            rep = &(sub->rep[j]);
            fprintf(flog,"pop.sub[%d].rep[%d].dataname = %s\n",i,j,rep->dataname);
        }
    }
    for (int i=0;i<pop->N_SUBS;i++) {
        sub = &(pop->sub[i]);
        for (int j=0;j<sub->N_REPS;j++) {
            rep = &(sub->rep[j]);
            fprintf(flog,"pop.sub[%d].rep[%d].designname = %s\n",i,j,rep->designname);
        }
    }
    fprintf(flog,"non_parm = %d\n",pop->non_parm);
    fprintf(flog,"Nb = %d\n",pop->Nb);
}

void load_data_structs(POP *pop,int PPP) 
{
    int N;
    int *dim_design,*dim_HRF;
    double yin,**HRF;
    char *S,*C;
    SUB *sub;
    REP *rep;
    FILE *fout;
    
    void itoa(int n,char s[]);
    double *create_bspline_basis(int N,double TIME,double freq,int int_knots,int *dim_HRF,double *knots,int type);
    double **canonical_HRF(int T,double freq,int *dim_HRF,int type);
    double *convolve(double **design,double **hrf,int *dim_design,int *dim_hrf);
    double *subsample_data(double *old_TS,int len,int *sublen,double true_freq,double samp_freq);
    double *subsample_design(double *mat,int len,int ncol,int *sublen,double true_freq,double samp_freq);
    void standardize_data(double *x,const int len);       

    S = (char *)calloc(300,sizeof(char));
    C = (char *)calloc(300,sizeof(char));
    dim_design = (int *)calloc(2,sizeof(int));
    dim_HRF = (int *)calloc(2,sizeof(int));
 
    pop->fout_reprec = fopen("./log/pop_reprec.log","w");
    pop->fout_beta = fopen("./log/pop_beta.log","w");

 /*** OPEN TIME SERIES FILE AND READ DATA ***/
 
    for (int isub=0;isub<pop->N_SUBS;isub++) {
        sub = &(pop->sub[isub]);
        for (int irep=0;irep<sub->N_REPS;irep++) {
            rep = &(sub->rep[irep]);
            
            rep->P = PPP;
            rep->df_delta1 = 0.999;
            rep->df_delta2 = 0.98;
      
            fprintf(flog,"%s\n",rep->dataname);fflush(flog);
            fout = fopen(rep->dataname,"r");
 
            N = 0;
            while (fscanf(fout,"%lf ",&(yin)) != EOF)
                N++;
           
            rep->Y = (double *)calloc(N,sizeof(double));
            rewind(fout);
            
            for (int i=0;i<N;i++)
                int ifs = fscanf(fout,"%lf ",&(rep->Y[i]));
            
            fclose(fout);
            
            fout = fopen(rep->designname,"r");
            fprintf(flog,"%s\n",rep->designname);fflush(flog);
            int cnt=0;
            while (fscanf(fout,"%lf ",&yin) != EOF)
                cnt++;
            
            rewind(fout);
            
            dim_design[0] = N;
            dim_design[1] = cnt/(N);
        
            double **design;
            design = (double **)calloc(N,sizeof(double *));
            for (int i=0;i<N;i++)
                design[i] = (double *)calloc(dim_design[1],sizeof(double));
            for (int i=0;i<N;i++)
                for (int j=0;j<dim_design[1];j++)
                    int ifs = fscanf(fout,"%lf ",&(design[i][j]));
        
            fclose(fout);

            // CONVOLVE DESIGN MATRIX WITH HRF
            HRF = canonical_HRF(35,sub->freq,dim_HRF,pop->Nb-1);
            pop->Nb = dim_HRF[1];
 
            rep->X = convolve(design,HRF,dim_design,dim_HRF);
  
            for (int i=0;i<N;i++)
                free(design[i]);
            free(design);
 
            // subsample data and design
            
            int sfreq = 10;
            if (sub->freq > sfreq) {
                int subN;

                double *XXX,*Y;        
                XXX = subsample_design(rep->X,N,dim_design[1]*dim_HRF[1],&subN,sub->freq,sfreq); 
                free(rep->X);
                rep->X = XXX;
 

                Y = subsample_data(rep->Y,N,&subN,sub->freq,sfreq);

                free(rep->Y);
                rep->Y = Y;
                sub->freq = sub->freq/round(sub->freq/sfreq);
                N = subN;
            }
            
            dim_design[0] = N;
            rep->dim_X = (int *)calloc(2,sizeof(int));
            rep->dim_X[0] = N;
            rep->dim_X[1] = dim_design[1]*dim_HRF[1];

            // standardize TS

            standardize_data(rep->Y,(const int)N);

            // rescale design matrix
            
            double max = 0;
            double tmp = 0;

            for (int j=0;j<rep->dim_X[1];j++) { // normalize
                tmp = 0;
                for (int i=0;i<rep->dim_X[0];i++)
                    tmp = (tmp > rep->X[i*rep->dim_X[1]+j]) ? tmp:rep->X[i*rep->dim_X[1]+j];
 //                  tmp += rep->X[i*rep->dim_X[1]+j]*rep->X[i*rep->dim_X[1]+j];
 //                  tmp = sqrt(tmp);
                max = (tmp > max) ? tmp:max;
            }
            for (int j=0;j<rep->dim_X[1];j++) { // normalize
                for (int i=0;i<rep->dim_X[0];i++)
                    rep->X[i*rep->dim_X[1]+j] /= max;
            }
            
            fprintf(flog,"max = %lf\n",max);fflush(stdout);      
           
  
            /*** CREATE Initial B-Spline BASIS TO REMOVE LOW FREQUENCY DRIFT ***/
            
            double T = (double)(N-1)/(double)pop->sub[0].freq;
            int K=0;
            
            pop->sub[isub].rep[irep].dim_V = (int *)calloc(2,sizeof(int));
            
            pop->sub[isub].rep[irep].nKnots = pop->knots + 2;
            pop->sub[isub].rep[irep].knots = (double *)calloc(pop->sub[isub].rep[irep].nKnots+6,sizeof(double));
            
            pop->sub[isub].rep[irep].V = create_bspline_basis(N,T,pop->sub[isub].freq,pop->sub[isub].rep[irep].nKnots,pop->sub[isub].rep[irep].dim_V,pop->sub[isub].rep[irep].knots,0);
            pop->sub[isub].rep[irep].nKnots += 6;                  
 

            fprintf(flog,"NKNOTS = %d, dim_V[1] = %d\n",pop->sub[isub].rep[irep].nKnots,pop->sub[isub].rep[irep].dim_V[1]);
            for (int i=0;i<pop->sub[isub].rep[irep].nKnots;i++)
                fprintf(flog,"%lf ",pop->sub[isub].rep[irep].knots[i]);
            fprintf(flog,"\n");
           
            
            /*** initialize TVAR TERMS ***/
            
            
            pop->sub[isub].rep[irep].dim_W = (int *)calloc(2,sizeof(int));
            pop->sub[isub].rep[irep].dim_W[0] = N;
            pop->sub[isub].rep[irep].dim_W[1] = PPP;
            
            pop->sub[isub].rep[irep].W = (double *)calloc(N*maxP,sizeof(double)); // needed to select DLM hyperpriors
            pop->sub[isub].rep[irep].tableP = (int *)calloc(maxP+1,sizeof(int));
            pop->sub[isub].rep[irep].md1 = 0;
            pop->sub[isub].rep[irep].md2 = 0;
            
            for (int i=PPP;i<N;i++) {
                for (int j=0;j<PPP;j++) {
                    pop->sub[isub].rep[irep].W[i*PPP + PPP-j-1] = pop->sub[isub].rep[irep].Y[i+j-PPP];
                }
            }
            
            /*** intialize remaining data ***/
            
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

            pop->sub[isub].rep[irep].mVeta = (double *)calloc(pop->sub[isub].rep[irep].dim_V[0],sizeof(double));
 
            pop->sub[isub].rep[irep].delta = (double *)calloc(pop->sub[isub].rep[irep].dim_W[0]*maxP,sizeof(double));
            pop->sub[isub].rep[irep].m = (double *)calloc(pop->sub[isub].rep[irep].dim_W[0]*maxP,sizeof(double));

            pop->sub[isub].rep[irep].mdelta = (double *)calloc(pop->sub[isub].rep[irep].dim_W[0]*maxP,sizeof(double));

            pop->sub[isub].rep[irep].mdelta2 = (double *)calloc(pop->sub[isub].rep[irep].dim_W[0]*maxP,sizeof(double));
 
            pop->sub[isub].rep[irep].Wdelta = (double *)calloc(pop->sub[isub].rep[irep].dim_W[0],sizeof(double));
            pop->sub[isub].rep[irep].Wm = (double *)calloc(pop->sub[isub].rep[irep].dim_W[0],sizeof(double));
            pop->sub[isub].rep[irep].mWdelta = (double *)calloc(pop->sub[isub].rep[irep].dim_W[0],sizeof(double));

            pop->sub[isub].rep[irep].Xbeta = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0],sizeof(double));
            pop->sub[isub].rep[irep].mXbeta = (double *)calloc(pop->sub[isub].rep[irep].dim_V[0],sizeof(double));
            pop->sub[isub].rep[irep].Hbeta = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0],sizeof(double));
            pop->sub[isub].rep[irep].H = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0]*pop->sub[isub].rep[irep].dim_X[1],sizeof(double));
            
            pop->sub[isub].rep[irep].precY = (double *)calloc(3,sizeof(double));
            pop->sub[isub].rep[irep].precY[0] = 1;
            pop->sub[isub].rep[irep].precY[1] = 1./16.;
            pop->sub[isub].rep[irep].precY[2] = 1./64.;
            pop->sub[isub].rep[irep].d_Y = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0],sizeof(double));
            pop->sub[isub].rep[irep].md_Y = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0],sizeof(double));
            pop->sub[isub].rep[irep].sd_Y = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0],sizeof(double));
            for (int i=0;i<pop->sub[isub].rep[irep].dim_X[0];i++)
                pop->sub[isub].rep[irep].d_Y[i] = 1;
             
            pop->sub[isub].rep[irep].K = 1;
            pop->sub[isub].rep[irep].preceta = 1;
            pop->sub[isub].rep[irep].precdelta = 100.;
            pop->sub[isub].rep[irep].LY = 1./16.;
            pop->sub[isub].rep[irep].precYstart = 1;

            pop->sub[isub].rep[irep].mean_res = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0],sizeof(double));
            pop->sub[isub].rep[irep].mean_d_Y = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0],sizeof(double));
            pop->sub[isub].rep[irep].mean_fit = (double *)calloc(pop->sub[isub].rep[irep].dim_X[0],sizeof(double));
           
            pop->sub[isub].rep[irep].prop_sd = (double *)calloc(4,sizeof(double));
            pop->sub[isub].rep[irep].prop_sd[0] = 0.005;
            pop->sub[isub].rep[irep].prop_sd[1] = 0.005;
            pop->sub[isub].rep[irep].prop_sd[2] = 0.001;
            pop->sub[isub].rep[irep].prop_sd[3] = 100;
            
            pop->sub[isub].rep[irep].accept = (int *)calloc(4,sizeof(int));
            pop->sub[isub].rep[irep].attempt = (int *)calloc(4,sizeof(int));
            
            /*** open log files for output ***/
            
            C = strcpy(C,pop->sub[isub].rep[irep].dataname);
            C = strtok(C,".");
            
            S = strcpy(S,"./log/");
            S = strcat(S,C);
            S = strcat(S,"_eta.log");
            pop->sub[isub].rep[irep].fout_eta = fopen(S,"w");
          
            S = strcpy(S,"./log/");
            S = strcat(S,C);
            S = strcat(S,"_DLM.log");
            pop->sub[isub].rep[irep].fout_dlm = fopen(S,"w+");

            S = strcpy(S,"./log/");
            S = strcat(S,C);
            S = strcat(S,"_res.log");
            pop->sub[isub].rep[irep].fout_res = fopen(S,"w");

            S = strcpy(S,"./log/");
            S = strcat(S,C);
            S = strcat(S,"_delta.log");
            pop->sub[isub].rep[irep].fout_delta = fopen(S,"w");

            S = strcpy(S,"./log/");
            S = strcat(S,C);
            S = strcat(S,"_nknots.log");
            pop->sub[isub].rep[irep].fout_nknots = fopen(S,"w");

            S = strcpy(S,"./log/");
            S = strcat(S,C);
            S = strcat(S,"_knots.log");
            pop->sub[isub].rep[irep].fout_knots = fopen(S,"w");

            S = strcpy(S,"./log/");
            S = strcat(S,C);
            S = strcat(S,"_prec.log");
            pop->sub[isub].rep[irep].fout_prec = fopen(S,"w");
            
            S = strcpy(S,"./log/");
            S = strcat(S,C);
            S = strcat(S,"_wdelta.log");
            pop->sub[isub].rep[irep].fout_wdelta = fopen(S,"w");
            
            S = strcpy(S,"./log/");
            S = strcat(S,C);
            S = strcat(S,"_veta.log");
            pop->sub[isub].rep[irep].fout_veta = fopen(S,"w");
 
            S = strcpy(S,"./log/");
            S = strcat(S,C);
            S = strcat(S,"_fit.log");
            pop->sub[isub].rep[irep].fout_fit = fopen(S,"w");
 
            S = strcpy(S,"./log/");
            S = strcat(S,C);
            S = strcat(S,"_X.log");
            pop->sub[isub].rep[irep].fout_X = fopen(S,"w");
          
            int nrow = pop->sub[isub].rep[irep].dim_X[0];
            int ncol = pop->sub[isub].rep[irep].dim_X[1];
            for (int i=0;i<nrow;i++) {
                for (int j=0;j<ncol;j++) {
                    fprintf(pop->sub[isub].rep[irep].fout_X,"%.20lf ",pop->sub[isub].rep[irep].X[i*ncol+j]);
                }
                fprintf(pop->sub[isub].rep[irep].fout_X,"\n");
            }
            fclose(pop->sub[isub].rep[irep].fout_X);
  
            S = strcpy(S,"./log/");
            S = strcat(S,C);
            S = strcat(S,"_Y.log");
            pop->sub[isub].rep[irep].fout_Y = fopen(S,"w");
          
            nrow = pop->sub[isub].rep[irep].dim_X[0];
            for (int i=0;i<nrow;i++) 
                    fprintf(pop->sub[isub].rep[irep].fout_Y,"%.20lf ",pop->sub[isub].rep[irep].Y[i]);
            fclose(pop->sub[isub].rep[irep].fout_Y);
  
            S = strcpy(S,"./log/");
            S = strcat(S,C);
            S = strcat(S,"_HRF.log");
            pop->sub[isub].rep[irep].fout_HRF = fopen(S,"w");
 
            nrow = dim_HRF[0];
            ncol = dim_HRF[1];
            for (int i=0;i<nrow;i++) {
                for (int j=0;j<ncol;j++)
                    fprintf(pop->sub[isub].rep[irep].fout_HRF,"%lf ",HRF[i][j]);
                fprintf(pop->sub[isub].rep[irep].fout_HRF,"\n");
            }
            fclose(pop->sub[isub].rep[irep].fout_HRF);
            
            for (int i=0;i<nrow;i++)
                free(HRF[i]);
            free(HRF);

            S = strcpy(S,"./log/");
            S = strcat(S,C);
            S = strcat(S,"_Xbeta.log");
            pop->sub[isub].rep[irep].fout_Xbeta = fopen(S,"w");
        }
        
        S = strcpy(S,"./log/");
        S = strcat(S,C);
        S = strcat(S,"_beta.log");
        pop->sub[isub].fout_beta = fopen(S,"w");
 
        pop->sub[isub].beta = (double *)calloc(pop->Nb*(pop->Ns),sizeof(double));
    }
 
    free(C);
    free(S);
    free(dim_design);
    free(dim_HRF);
    
    /*** create DLM structures ***/
    
    int max_dimX = 0;
    for (int i=0;i<pop->N_SUBS;i++)
        for (int j=0;j<pop->sub[i].N_REPS;j++)
            max_dimX = (max_dimX > pop->sub[i].rep[j].dim_X[0]) ? max_dimX:pop->sub[i].rep[j].dim_X[0];
 
    dlmStruc = (sDLM *)calloc(max_dimX,sizeof(sDLM));
    for (int i=0;i<max_dimX;i++)
        dlmStruc[i].a = (double *)calloc(maxP,sizeof(double));
    for (int i=0;i<max_dimX;i++)
        dlmStruc[i].m = (double *)calloc(maxP,sizeof(double));
    for (int i=0;i<max_dimX;i++)
        dlmStruc[i].R = (double *)calloc(maxP*maxP,sizeof(double));
    for (int i=0;i<max_dimX;i++)
        dlmStruc[i].C = (double *)calloc(maxP*maxP,sizeof(double));
     for (int i=0;i<max_dimX;i++)
        dlmStruc[i].S = 1.;
    for (int i=0;i<max_dimX;i++)
        dlmStruc[i].n = 1.;
    for (int i=0;i<max_dimX;i++)
        dlmStruc[i].d = 1.;

}


double *subsample_data(double *old_TS,int len,int *sublen,double true_freq,double samp_freq)
{
    int sub_samp;
    double *new_TS;
    
    sub_samp = (int)round(true_freq/samp_freq);
    *sublen = 0;
    for (int i=0;i<len;i+=sub_samp)
        (*sublen)++;
    new_TS = (double *)calloc(*sublen,sizeof(double));
    int n = 0;
    for (int i=0;i<len;i+=sub_samp) {
        new_TS[n] = old_TS[i];
        n++;
    }
    return new_TS;        
}

double *subsample_design(double *mat,int len,int ncol,int *sublen,double true_freq,double samp_freq)
{
    double *ts,*subsampMat;
    double *new_TS;
    
    ts = (double *)calloc(len,sizeof(double));
    for (int j=0;j<ncol;j++) {
        for (int i=0;i<len;i++)
            ts[i] = mat[i*ncol+j];
        new_TS = subsample_data(ts,len,sublen,true_freq,samp_freq);
        if (j==0) {
            subsampMat = (double *)calloc(*sublen*ncol,sizeof(double));
        }
        for (int i=0;i<*sublen;i++)
            subsampMat[i*ncol+j] = new_TS[i];
        free(new_TS);
    }
    free(ts);
    return subsampMat;
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
