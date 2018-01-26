
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
    double *XpX;
    double *lambda;
    double *H;
    double *J;
    double *Hbeta;
    double *Jeta;
      
    double *mdelta;
    double *mdelta2; 
    double *mean_Y;
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
    double *prbs;
} REP;


typedef struct subdata{
    int N_REPS;
    double freq;
    REP *rep;
    double *beta;  // subject level effects
    double *conv_beta;  // converted subject level effects
//    double precbeta;  // event level r.e. precision, one for each stimulus
    double *conv_precbeta;  // r.e. precision
} SUB;


typedef struct popdata{
    int GRP;
    int N_SUBS;
    int Nb;     // number of basis in the HRF
    int Ns;      // number of stimuli (or events)
    int non_parm; // canonical HRF = 0, non-parm B-SPlINE = 1
    double expo;
    double ED;
    double **Q; // smoothing prec matrix for stimulus level r.e.
    double **Q2; // smoothing prec matrix for low level drift
    SUB *sub;
    double *beta;  // population level effects, one for each stimulus
    double *re_prec; // sub level random effects, one for each stimulus
    double priorprec_beta;
} POP;
