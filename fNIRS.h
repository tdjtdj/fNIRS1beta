
typedef struct dlm{
    double S;
    double *a;
    double *m;
    double *C;
    double *R;
} sDLM;

typedef struct replication{
    char *dataname;
    char *designname;
    char *motionname;
    double K;
    double *pi_Y;
    double pi_eta;
    double pi_delta;
    int *Ne;
    int *dim_X;
    int *dim_W;
    int *dim_V;
    int P;
    double df_delta1;
    double df_delta2;
    double *d_Y;
    double *md_Y;
    double *sd_Y;
    double *Y;
    double *precY;
    double LY;
    double precYstart;
    double *eta;
    double preceta;
    double *delta;
    double precdelta;
    double *Xbeta;
    double *Wdelta;
    double *mWdelta;
    double *X;
    double *W;
    double *V;
    double *Veta;
    double *mVeta;
    double *lambda;
    double *H;
    double *J;
    double *Hbeta;
    double *Jeta;
      
    double *mdelta;
    double *mdelta2; 
    double *mean_res;
    double *mean_d_Y;
    double *mean_fit;
    double mean_precY;
    double *residuals;
    double *residuals1;
    double *residuals2;
    double *residuals3;
    double *residuals4;
    double *residuals5;
        
    double *prop_sd;
    int *attempt;
    int *accept;
    int maxsteps;
    int nKnots;
    double *knots;
    
    FILE *fout_eta;
    FILE *fout_dlm;
    FILE *fout_res;
    FILE *fout_delta;
    FILE *fout_nknots;
    FILE *fout_knots;
    FILE *fout_prec;
    FILE *fout_wdelta;
    FILE *fout_veta;
    FILE *fout_fit;
} REP;


typedef struct subdata{
    int N_REPS;
    double freq;
    REP *rep;
    double *beta;  // subject level effects  (Nb*Ns)
    double *conv_beta;  // converted subject level effects
//    double precbeta;  // event level r.e. precision, one for each stimulus
    double *conv_precbeta;  // r.e. precision
    int *dim_X;
    double *X; // (Nb*Ns)x(Ncov*Nb*Ns)  covariate information for each subject
    
    FILE *fout_beta;
} SUB;


typedef struct popdata{
    int GRP;
    int N_SUBS;
    int Ncov; // number of covariates, defaults to 1 for intercept
    int Nb;     // number of basis in the HRF
    int Ns;      // number of stimuli (or events)
    int non_parm; // canonical HRF = 0, non-parm B-SPlINE = 1
    int knots; // number of bspline bases
    double ED;
    SUB *sub;
    double *beta;  // population level effects, (Ncov*Nb*Ns)
    double *re_prec; // sub level random effects, (Ns)x(Ns) diagonal
    
    FILE *fout_reprec;
    FILE *fout_beta;
} POP;
