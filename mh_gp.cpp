/*
 *  mh_gp.c
 *  LGCP
 *
 *  Created by Timothy Johnson on 3/11/10.
 *  Copyright 2010 University of Michigan. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LGCP.h"
#include "/usr/local/include/fftw3.h"
#include "randgen.h"


extern int *nidx,*midx;
extern fftw_complex *XYZ;

/*float likelihood_tst2(SUB *sub,GP *gp,float cell_vol,int ntot)
{
	int i,j;
	double yprop=0,expyl=0,y;
	
	for (i=0;i<ntot;i++) {
		j = *(midx+i);
		y = (double)(*(sub->mu) + gp[(sub->subtype)+1].Y.realp[j]);
		yprop += y*(sub->ncell[j]);
		expyl += exp(y); 		
		
	}
	
	return (float)(yprop - expyl*cell_vol);	
}
*/
double like_pre(GP *gp,int ntot)
{
	int i,j;
	double yprop=0,expyl=0,y;
	
	for (i=0;i<ntot;i++) {
		j = *(midx+i);
		y = (double)(gp->Y.realp[j]);
		expyl += exp(y); 		
		
	}
	
	return expyl;	

}

float likelihood_tst2(SUB *sub,GP *gp,double expyl,float cell_vol,int ntot)
{
	int i,j;
	double yprop=0,y;
	
	for (i=0;i<ntot;i++) {
		j = *(midx+i);
		y = (double)(*(sub->mu) + gp->Y.realp[j]);
		yprop += y*(sub->ncell[j]);		
	}
	
	return (float)(yprop - exp(*sub->mu)*expyl*cell_vol);	
}



void gamma_gradient(GP *gp,int nsub,float *prop,float mu,int mtot,int *n,int *m,int *g,float cell_volume,float *prop_gr,
				   float *B,float hdiv2,int subtype,int nsubtypes,int dim)
{
	int i,j,ii,isubtype;
	
	void Cx(float *,int);
	for (i=0;i<mtot;i++) 
		XYZ[i][0] = XYZ[i][1] = 0;

	int ntot=1;
	for (i=0;i<dim;i++)
		ntot *= n[i];

	for (i=0;i<ntot;i++) {
		j = *(midx+i);
		XYZ[j][0] = gp[subtype].nn[j] - mu*expf(gp[subtype].Y.realp[j])*cell_volume;			
	}

	Cx(gp[subtype].lambda,mtot);	

	for (i=0;i<mtot;i++)
		XYZ[i][0] *= *(gp[subtype].sigma);
	for (i=0;i<mtot;i++)
		prop_gr[i] = (XYZ[i][0] - prop[i])*hdiv2;
	
}



void MH_gamma(GP *gp,SUB *sub,int nsub,float pgprec,float *prop,float cell_volume,float *sd,int *n,int *m,int *g,int ntot,int mtot,
			  float *B,float *accept,float *attempt,float **yy,int nsubtypes,unsigned long *seed,int iter,int dim)
{
	int i,j,isubtype;
	double prop_dotprd,prop_dotprd2, prop_dotprd3,tmp;
	float new_like;
	float like_ratio,prior_ratio,prop_ratio,*mu;
	float h;// = sd*sd;
	float hdiv2;// = 0.5f*h;
	float *old;
	float *prop_gr,*newlike;
	double expyl=0;
	void Cx(float *lambda,int mtot);
	float likelihood_tst3(GP *gp,float cell_vol,int ntot,int subtype,int nsub);

	old = (float *)calloc(mtot,sizeof(float));
	newlike = (float *)calloc(nsub,sizeof(float));
	mu = (float *)calloc(nsubtypes+1,sizeof(float));
	for (isubtype=1;isubtype<nsubtypes+1;isubtype++) {
		for (i=0;i<nsub;i++) {
			if (sub[i].subtype == isubtype-1)
				mu[isubtype] = expf(*(gp[isubtype].mu));
		}
	}	
	prop_gr = (float *)calloc(mtot,sizeof(float));
	
	for (isubtype=1;isubtype<nsubtypes+1;isubtype++) {
		h = sd[isubtype]*sd[isubtype];
		hdiv2 = 0.5f*h;
		gamma_gradient(gp,*(gp[isubtype].nsubs),gp[isubtype].gamma,mu[isubtype],mtot,n,m,g,cell_volume,gp[isubtype].gamma_gr,B,hdiv2,isubtype,nsubtypes,dim);	
		for (i=0;i<mtot;i++) {
			prop[i] = (float)rnorm((double)(gp[isubtype].gamma[i]+gp[isubtype].gamma_gr[i]),(double)sd[isubtype],seed);// u^(m_+1)
		}
		for (i=0;i<mtot;i++) {
			XYZ[i][0] = prop[i];
			XYZ[i][1] = 0;
		}
		Cx(gp[isubtype].lambda,mtot); // returns y = C^(1/2) lambda in U.realp (y is the proposed gaussian process) 
		
		for (i=0;i<mtot;i++) {
			old[i] = gp[isubtype].Y.realp[i];
			gp[isubtype].Y.realp[i] = *(gp[isubtype].sigma)*XYZ[i][0];
		}
		gamma_gradient(gp,*(gp[isubtype].nsubs),prop,mu[isubtype],mtot,n,m,g,cell_volume,prop_gr,B,hdiv2,isubtype,nsubtypes,dim);	
		prop_dotprd=0.0;prop_dotprd2=0.0;prop_dotprd3=0.0;
		
		for (i=0;i<mtot;i++) {
			tmp = (double)(gp[isubtype].gamma[i]-prop[i]-prop_gr[i]);
			prop_dotprd2 += tmp*tmp;			// ||r^m - e(u^(m+1)) ||^2  
			tmp = (double)(prop[i]-gp[isubtype].gamma[i]-gp[isubtype].gamma_gr[i]);
			prop_dotprd3 += tmp*tmp;		  // ||u^(m+1) - e(r^m) ||^2 
			prop_dotprd  += (double)((prop[i])*(prop[i]));												  		// (u^(m+1))^2
		}
		
		prior_ratio = - 0.5*((float)prop_dotprd - *(gp[isubtype].gamma_dotprd));
		//		prop_ratio = gamma_prop_ratio(prop,prop_gr,gp[isubtype].gamma,gp[isubtype].gamma_gr,h,mtot);
		prop_ratio = -(1.0/(2.0*h))*(float)(prop_dotprd2 - prop_dotprd3);
		
		new_like = 0;
		for (i=0;i<nsub;i++) {
			newlike[i] = 0;
			if (isubtype > 0) {
				if (sub[i].subtype == isubtype-1) {
					newlike[i] = likelihood_tst2(&(sub[i]),&gp[isubtype],expyl,cell_volume,ntot);
					new_like += newlike[i];
				}
			}
			else {
				newlike[i] = likelihood_tst2(&(sub[i]),&gp[isubtype],expyl,cell_volume,ntot);
				new_like += newlike[i];
			}
		}
		
//		new_like = likelihood_tst3(gp,cell_volume,ntot,isubtype,*(gp[isubtype].nsubs));

		like_ratio = new_like - *(gp[isubtype].likehd); // likelihood ratio 
		
	//	printf("%lf %lf \t %lf \t %lf\n",new_like,*(gp[isubtype].likehd),prior_ratio,prop_ratio);
		if ((float)log(kiss(seed)) < prior_ratio + like_ratio + prop_ratio) {
			for (i=0;i<mtot;i++) {
				gp[isubtype].gamma[i] = prop[i];
				gp[isubtype].gamma_gr[i] = prop_gr[i];
			}
			*(gp[isubtype].gamma_dotprd) = prop_dotprd;              // copy proposed gamma dot product to gamma dot product
			if (isubtype>0) {
				*(gp[0].likehd) += like_ratio;
			}
			*(gp[isubtype].likehd) = new_like;                 // update current likelihood
			for (i=0;i<nsub;i++) {
				if (isubtype > 0) {
					if (sub[i].subtype == isubtype-1)
						*(sub[i].likehd) = newlike[i];
				}
				else {
					*(gp[sub[i].subtype+1].likehd) += newlike[i] - *(sub[i].likehd);
					*(sub[i].likehd) = newlike[i];
					
				}
			}
			accept[isubtype]++;
			for (i=0;i<ntot;i++) {
				j = *(midx+i);
				yy[isubtype][i] = *(gp[isubtype].mu) + gp[isubtype].Y.realp[j];
			}
/*			if (dim == 2) {
				int h = 0;
				for(j=0;j<n[1];j++) {
					for(i=0;i<n[0];i++) {
						if (isubtype > 0) 
							yy[isubtype][h] = *(gp[isubtype].mu) + gp[isubtype].Y.realp[i + j*m[0]];
					//	else
					//		yy[isubtype][h] = gp[0].Y.realp[i + j*m[0]];
						
						
						h++;
					}
					
				}
			}
			if (dim == 3) {
				int h = 0;
				for (int k=0;k<n[2];k++) {
					for(j=0;j<n[1];j++) {
						for(i=0;i<n[0];i++) {
							if (isubtype > 0) 
								yy[isubtype][h] = *(gp[isubtype].mu) + gp[isubtype].Y.realp[i + j*m[0] + k*m[0]*m[1]];
						//	else
						//		yy[isubtype][h] = gp[0].Y.realp[i + j*m[0] + k*m[0]*m[1]];
							
							
							h++;
						}
					}	
				}
			}*/
			//		if (iter > 16434 & isubtype == 1)	printf("\t \t accept %.8f %.8f\n",new_like, *(gp[isubtype].likehd));
		} 
		else {
			for (i=0;i<mtot;i++) 
				gp[isubtype].Y.realp[i] = old[i];
			for (i=0;i<ntot;i++) {
				j = *(midx+i);
				yy[isubtype][i] = *(gp[isubtype].mu) + gp[isubtype].Y.realp[j];
			}
	/*		if (dim == 2) {
				int h = 0;
				for(j=0;j<n[1];j++) {
					for(i=0;i<n[0];i++) {
						if (isubtype > 0) 
							yy[isubtype][h] = *(gp[isubtype].mu) + gp[isubtype].Y.realp[i + j*m[0]];
						//	else
						//		yy[isubtype][h] = gp[0].Y.realp[i + j*m[0] + k*m[0]*m[1]];
						
						
						h++;
					}
				}	
			}
			if (dim == 3) {
				int h = 0;
				for (int k=0;k<n[2];k++) {
					for(j=0;j<n[1];j++) {
						for(i=0;i<n[0];i++) {
							if (isubtype > 0) 
								yy[isubtype][h] = *(gp[isubtype].mu) + gp[isubtype].Y.realp[i + j*m[0] + k*m[0]*m[1]];
							//	else
							//		yy[isubtype][h] = gp[0].Y.realp[i + j*m[0] + k*m[0]*m[1]];
							
							
							h++;
						}
					}	
				}
			}*/
			
			//			if (iter > 16434 & isubtype == 1)	printf("\t \t reject %.8f %.8f\n",new_like, *(gp[isubtype].likehd));
		}
		attempt[isubtype]++;
	}
	//	if (iter > 16434 & isubtype == 1)
	//		printf("\n");
	free(prop_gr);

	free(mu);
	
	free(newlike);
	free(old);
}


void MH_mu(GP *gp,SUB *sub,int nsub,float cell_volume,int *n,int *m,int mtot,float *accept,float *attempt,
		   float *sd,int nsubtypes,int dim,unsigned long *seed)
{
	int isub,isubtype,ntot;
	float prop,like_ratio,prior_ratio,new_like,newlike;
	float likelihood_tst3(GP *gp,float cell_vol,int ntot,int subtype,int nsub);
	
	ntot = 1;
	for (int i=0;i<dim;i++)
		ntot *= n[i];
	
	*(gp[0].likehd) = 0;
	for (isubtype=1;isubtype<nsubtypes+1;isubtype++) {
		//	for (isub=0;isub<nsub;isub++) {
		//		if (sub[isub].subtype == isubtype-1) {
		prop = *(gp[isubtype].mu);
		*(gp[isubtype].mu) = (float)rnorm((double)prop,(double)(sd[isubtype]),seed);
		new_like = likelihood_tst3(gp,cell_volume,ntot,isubtype,*(gp[isubtype].nsubs));
		
		like_ratio = new_like - *(gp[isubtype].likehd);
		
		prior_ratio = 0;
//		printf("%f %f\n",*gp[isubtype].mu,like_ratio+prior_ratio);
		if ((float)log(kiss(seed)) < like_ratio + prior_ratio) {
			if (isubtype>0) {
				*(gp[0].likehd) += like_ratio;
			}
			*(gp[isubtype].likehd) = new_like;
			accept[isubtype]++;
		}
		else {
			*(gp[isubtype].mu) = prop;
		}
		
		attempt[isubtype]++;
		//		}
	//	}
	}
}


void MH_beta(GP *gp,SUB *sub,int nsub,float cell_volume,int *n,int *m,int mtot,float *accept,float *attempt,float sd,int *mask,unsigned long *seed)
{
	int i;
	float prop,like_ratio,prior_ratio,new_like,*newlike,*Utmp;
	
//	int ntot = n[0]*n[1];
	
	Utmp = (float *)calloc(mtot,sizeof(float));
	newlike = (float *)calloc(nsub,sizeof(float));

	prop = (float)rnorm((double)(*(gp->beta)),(double)sd,seed);
	for (i=0;i<mtot;i++) {
		gp->Y.imagp[i] = gp->Y.realp[i] + prop*gp->prior_mean[i];
		Utmp[i] = expf(gp->Y.imagp[i]);
	}
	
	new_like = 0;
	for (i=0;i<nsub;i++) {
//		vDSP_vsadd(gp->Y.imagp,1,sub[i].mu,Utmp,1,mtot);
//		newlike[i] =  likelihood_tst(sub[i].ncell,cell_volume,gp->Y.imagp,Utmp,n,m,*(sub[i].mu),expf(*(sub[i].mu)),ntot);
	}
	new_like = 0;
	for (i=0;i<nsub;i++)
		new_like += newlike[i];
		
	like_ratio = new_like - *(gp->likehd);
	prior_ratio = 0;//-0.5*.0001*((prop-0)*(prop-0) - (*(sub->truebeta)-0)*(*(sub->truebeta)-0));
	if ((float)log(kiss(seed)) < like_ratio + prior_ratio) {
		*(gp->beta) = prop;
		*(gp->likehd) = new_like;
		for (i=0;i<nsub;i++)
			*(sub[i].likehd) = newlike[i];
		*accept += 1;
	}
	*attempt += 1;
	free(Utmp);
	free(newlike);
}

void MH_sigma(GP *gp,SUB *sub,int nsub,float psigma,float psigmaprec,float *B,float cell_volume,int *n,int *m,int mtot,
			  float *sd,float *accept,float *attempt,int nsubtypes,int dim,unsigned long *seed)
{
	int i,isubtype,ntot;
	float prop,yn,like_ratio,prior_ratio,new_like,*newlike;
	float *old;
	float likelihood_tst3(GP *gp,float cell_vol,int ntot,int subtype,int nsub);
//	float ALPHA = 3.0, BETA = 2.0;

	if (dim == 2)
		ntot = n[0]*n[1];
	if (dim == 3)
		ntot = n[0]*n[1]*n[2];
		
	old = (float *)calloc(mtot,sizeof(float));
	
	for (isubtype=1;isubtype<nsubtypes+1;isubtype++) {
		newlike = (float *)calloc(nsub,sizeof(float));
		prop = *(gp[isubtype].sigma);
		*(gp[isubtype].sigma) = (float)rnorm((double)(*(gp[isubtype].sigma)),(double)sd[isubtype],seed);
//		if (*(gp[isubtype].sigma) > 0 ) {
			for (i=0;i<mtot;i++) 
				old[i] = gp[isubtype].Y.realp[i];

			yn = (*(gp[isubtype].sigma))/prop;
			for (i=0;i<mtot;i++)
				gp[isubtype].Y.realp[i] = yn*gp[isubtype].Y.realp[i];
			
/*			for (i=0;i<nsub;i++) {
				newlike[i] = 0;
				if (isubtype > 0) {
					if (sub[i].subtype == isubtype-1) 
						newlike[i] = likelihood_tst2(&(sub[i]),gp,cell_volume,ntot);
				}
				else
					newlike[i] = likelihood_tst2(&(sub[i]),gp,cell_volume,ntot);
			}
			new_like = 0;
			for (i=0;i<nsub;i++) {
				if (isubtype > 0) {
					if (sub[i].subtype == isubtype-1)
						new_like += newlike[i];
				}
				else 
					new_like += newlike[i];
			}*/
		
			new_like = likelihood_tst3(gp,cell_volume,ntot,isubtype,*(gp[isubtype].nsubs));
			
			like_ratio = new_like - *(gp[isubtype].likehd);
			
			prior_ratio = 0;// (ALPHA-1)*log(prop/(*(gp->sigma))) + BETA*(*(gp->sigma) - prop);
			
			//		 printf("isubtype = %d %f %f %f %f\n",isubtype,like_ratio,new_like,*(gp[isubtype].likehd),alpha);
			//		printf("\t %f %f\n",prop,*(gp[isubtype].sigma));
			if ((float)log(kiss(seed)) < like_ratio + prior_ratio) {
				if (isubtype>0) {
					*(gp[0].likehd) += like_ratio;
				}
				*(gp[isubtype].likehd) = new_like;                 // update current likelihood
		/*		for (i=0;i<nsub;i++) {
					if (isubtype > 0) {
						if (sub[i].subtype == isubtype-1)
							*(sub[i].likehd) = newlike[i];
					}
					else {
						*(gp[sub[i].subtype+1].likehd) += newlike[i] - *(sub[i].likehd);
						*(sub[i].likehd) = newlike[i];
						
					}
				}*/
				accept[isubtype] += 1;
			}
			else {
				*(gp[isubtype].sigma) = prop;
				for (i=0;i<mtot;i++) 
					gp[isubtype].Y.realp[i] = old[i];
				
				//			yn = 1.0/yn;
				//			vDSP_vsmul(gp[isubtype].Y.realp,1,&yn,gp[isubtype].Y.realp,1,mtot);
			}
			attempt[isubtype] += 1;
//		}
//		else {
//			*(gp[isubtype].sigma) = prop;
//			attempt[isubtype] +=1;
//		}
		free(newlike);
	}
	free(old);
}


void MH_rho(GP *gp,SUB *sub,int nsub,float *B,float cell_volume,int *n,int *m,int mtot,int *g,float *lattice,float *expo,
			float *sd,float *accept,float *attempt,int *count,int nsubtypes,int dim,unsigned long *seed)
{
	int i,flag,isubtype,ntot;
	float *lambdanew,prop,like_ratio,prior_ratio,new_like,*newlike,prop_ratio;
	double expyl;
//	float ALPHA = 1,BETA=1;
//  prior is a G(0,0), i.e. \pi(rho) \propto 1/rho;
//  With this prior and a log-normal proposal, the proposal ratio cancels the prior ratio.
	
	void block_circ_eigenvalues(float *lambda,int mtot,int *flag,int *count);	
	void Cx(float *lambda,int mtot);
	float likelihood_tst3(GP *gp,float cell_vol,int ntot,int subtype,int nsub);

	ntot = 1;
	for (i=0;i<dim;i++)
		ntot *= n[i];
	
	for (isubtype=1;isubtype<nsubtypes+1;isubtype++) {
		prop = *gp[isubtype].rho;
		*gp[isubtype].rho = expf((float)rnorm((double)(logf(prop)),(double)sd[isubtype],seed));
//		*(gp[isubtype].rho) = (float)rnorm((double)prop,(double)sd[isubtype],seed);

		if ((*(gp[isubtype].rho) < 1000.0) && (*(gp[isubtype].rho) > 0.001) ){
			lambdanew = (float *)calloc(mtot,sizeof(float));
			newlike = (float *)calloc(nsub,sizeof(float));

			for (i=0;i<mtot;i++) {
				XYZ[i][0] = expf(-(*gp[isubtype].rho)*gp[isubtype].lattice[i]);
				XYZ[i][1] = 0;
			}
			
			block_circ_eigenvalues(lambdanew,mtot,&flag,count);
//			block_circ_eigenvalues(*(gp[isubtype].rho),lambdanew,lattice,mtot,&flag,count);
//			printf("flag = %d\n",flag);//flag = 1;
			if (flag) {
				
				//				vDSP_vindex(gp[isubtype].gamma,B,1,U.realp,1,mtot);   // copy gamma to U.realp
				//				vDSP_vclr(U.imagp,1,mtot);
				for (i=0;i<mtot;i++) {
					XYZ[i][0] = gp[isubtype].gamma[i];
					XYZ[i][1] = 0;
				}
				Cx(lambdanew,mtot);  // returns y = C^(1/2) lambda in U.realp (y is the proposed gaussian process) 
				
				//			for (i=0;i<10;i++)
				//				printf("%f %f\n",lambdanew[i],gp[isubtype].lambda[i]);
				//				vDSP_vindex(gp[isubtype].Y.realp,B,1,U.imagp,1,mtot);   // copy Y to U.imagp
				//				vDSP_vsmul(U.realp,1,gp[isubtype].sigma,gp[isubtype].Y.realp,1,mtot);   // mult process by sigma
				for (i=0;i<mtot;i++) {
					XYZ[i][1] = gp[isubtype].Y.realp[i];
					gp[isubtype].Y.realp[i] = *(gp[isubtype].sigma)*XYZ[i][0];
				}
				
				new_like = 0;
				expyl = like_pre(&gp[isubtype],ntot);
				for (i=0;i<nsub;i++) {
					newlike[i] = 0;
					if (isubtype > 0) {
						if (sub[i].subtype == isubtype-1) {
							newlike[i] = likelihood_tst2(&(sub[i]),&gp[isubtype],expyl,cell_volume,ntot);
							new_like += newlike[i];
						}
					}
					else {
						newlike[i] = likelihood_tst2(&(sub[i]),&gp[isubtype],expyl,cell_volume,ntot);
						new_like += newlike[i];
					}
				}
				
			//	new_like = likelihood_tst3(gp,cell_volume,ntot,isubtype,*(gp[isubtype].nsubs));

				
				like_ratio = new_like - *(gp[isubtype].likehd);
				
//				prior_ratio = (ALPHA-1)*logf((*(gp[isubtype].rho))/prop) + BETA*(prop - *(gp[isubtype].rho));
				prior_ratio = 0;
				prop_ratio = logf(*(gp[isubtype].rho)/prop);  /* note: old value of rho is in prop */
//				prop_ratio = 0;
				//			printf("%f %f %f %f\n",new_like,like_ratio,prior_ratio,prop_ratio);
//				printf("%f %f \t %f %f\n",*gp[isubtype].rho, like_ratio + prior_ratio + prop_ratio,new_like,*(gp[isubtype].likehd));
				if ((float)log(kiss(seed)) < like_ratio + prior_ratio + prop_ratio) {
//				if ((float)log(kiss(seed)) < -100000000.) {
					if (isubtype>0) {
						*(gp[0].likehd) += like_ratio;
					}
					*(gp[isubtype].likehd) = new_like;                 // update current likelihood
					for (i=0;i<nsub;i++) {
						if (isubtype > 0) {
							if (sub[i].subtype == isubtype-1)
								*(sub[i].likehd) = newlike[i];
						}
						else {
							*(gp[sub[i].subtype+1].likehd) += newlike[i] - *(sub[i].likehd);
							*(sub[i].likehd) = newlike[i];
						}
					}
					//					vDSP_vindex(lambdanew,B,1,gp[isubtype].lambda,1,mtot); 
					for (i=0;i<mtot;i++) {
						gp[isubtype].lambda[i] = lambdanew[i];
					}
					accept[isubtype]++;
				}
				else {
					*(gp[isubtype].rho) = prop;
					//					vDSP_vindex(U.imagp,B,1,gp[isubtype].Y.realp,1,mtot); 
					for (i=0;i<mtot;i++) {
						gp[isubtype].Y.realp[i] = XYZ[i][1];
						XYZ[i][1] = 0;
					}
				}
//				*(gp[isubtype].rho) = prop;
//				printf("gp %f\n",*(gp[isubtype].rho));
			}
			else {
				gp[isubtype].rho[0] = prop;
			}
//			printf("OK\n");fflush(NULL);
			free(lambdanew);
			//printf("OK\");fflush(NULL);
			free(newlike);
		}
		else {
			gp[isubtype].rho[0] = prop;
		}
		attempt[isubtype]++;
	}
}


void MH_expo(GP *gp,SUB *sub,int nsub,float *B,float cell_volume,int *n,int *m,int mtot,int *g,float *lattice,float *expo,
			float *sd,float *accept,float *attempt,int *count,int nsubtypes,int dim,unsigned long *seed)
{
	int i,flag,isubtype,ntot;
	float *lambdanew,prop,like_ratio,prior_ratio,new_like,*newlike,prop_ratio;
	double expyl;
//	float ALPHA = 1,BETA=1;
//  prior is a G(0,0), i.e. \pi(rho) \propto 1/rho;
//  With this prior and a log-normal proposal, the proposal ratio cancels the prior ratio.
	
	void block_circ_eigenvalues(float *lambda,int mtot,int *flag,int *count);	
	void Cx(float *lambda,int mtot);
	float likelihood_tst3(GP *gp,float cell_vol,int ntot,int subtype,int nsub);
	float *new_lattice,prop_expo;
	double *dist;	
	void compute_dist_mat(double *dist,int dim,int *m,int *n);


	ntot = 1;
	for (i=0;i<dim;i++)
		ntot *= n[i];

	for (isubtype=1;isubtype<nsubtypes+1;isubtype++) {
		prop_expo = (float)rnorm((double)*gp[isubtype].delta,(double)sd[isubtype],seed);
		if ((prop_expo >= 0.5f) && (prop_expo < 1.99f)) {

			dist = (double *)calloc(mtot,sizeof(double));
			new_lattice = (float *)calloc(mtot,sizeof(float));
			compute_dist_mat(dist,dim,m,n);
			for (i=0;i<mtot;i++) 
				new_lattice[i] = (float)pow(dist[i],prop_expo);
			free(dist);

			lambdanew = (float *)calloc(mtot,sizeof(float));
			newlike = (float *)calloc(nsub,sizeof(float));

			for (i=0;i<mtot;i++) {
				XYZ[i][0] = expf(-(*gp[isubtype].rho)*new_lattice[i]);
				XYZ[i][1] = 0;
			}
			
			block_circ_eigenvalues(lambdanew,mtot,&flag,count);

			if (flag) {
			
				for (i=0;i<mtot;i++) {
					XYZ[i][0] = gp[isubtype].gamma[i];
					XYZ[i][1] = 0;
				}
				Cx(lambdanew,mtot);  // returns y = C^(1/2) lambda in U.realp (y is the proposed gaussian process) 
		
				for (i=0;i<mtot;i++) {
					XYZ[i][1] = gp[isubtype].Y.realp[i];
					gp[isubtype].Y.realp[i] = *(gp[isubtype].sigma)*XYZ[i][0];
				}
				
				new_like = 0;
				expyl = like_pre(&gp[isubtype],ntot);
				for (i=0;i<nsub;i++) {
					newlike[i] = 0;
					if (isubtype > 0) {
						if (sub[i].subtype == isubtype-1) {
							newlike[i] = likelihood_tst2(&(sub[i]),&gp[isubtype],expyl,cell_volume,ntot);
							new_like += newlike[i];						}
					}
					else {
						newlike[i] = likelihood_tst2(&(sub[i]),&gp[isubtype],expyl,cell_volume,ntot);
						new_like += newlike[i];
					}
				}	

				like_ratio = new_like - *(gp[isubtype].likehd);
					
				prior_ratio = 0;
				prop_ratio = 0;
//				prop_ratio = logf(*(gp[isubtype].rho)/prop);  /* note: old value of rho is in prop */

				if ((float)log(kiss(seed)) < like_ratio + prior_ratio + prop_ratio) {
					if (isubtype>0) {
						*(gp[0].likehd) += like_ratio;
					}
					*(gp[isubtype].likehd) = new_like;                 // update current likelihood
					for (i=0;i<nsub;i++) {
						if (isubtype > 0) {
							if (sub[i].subtype == isubtype-1)
								*(sub[i].likehd) = newlike[i];
						}
						else {
							*(gp[sub[i].subtype+1].likehd) += newlike[i] - *(sub[i].likehd);
							*(sub[i].likehd) = newlike[i];
						}
					}
	
					for (i=0;i<mtot;i++) {
						gp[isubtype].lambda[i] = lambdanew[i];
						gp[isubtype].lattice[i] = new_lattice[i];
					}
					accept[isubtype]++;
					*gp[isubtype].delta = prop_expo;
				}
				else {
					for (i=0;i<mtot;i++) {
						gp[isubtype].Y.realp[i] = XYZ[i][1];
						XYZ[i][1] = 0;
					}
				}
			}
			free(lambdanew);
			free(newlike);
			free(new_lattice);
		}
		attempt[isubtype]++;
	}
}
				

/*
void MH_sigma(GP *gp,SUB *sub,int nsub,float psigma,float psigmaprec,float *B,float cell_volume,int *n,int *m,int mtot,
			  float *sd,float *accept,float *attempt,int nsubtypes,int dim,unsigned long *seed)
{
	int i,isubtype,ntot;
	float prop,yn,like_ratio,prior_ratio,new_like,*newlike;
	float *old;
	//	float ALPHA = 3.0, BETA = 2.0;
	
	if (dim == 2)
		ntot = n[0]*n[1];
	if (dim == 3)
		ntot = n[0]*n[1]*n[2];
	
	old = (float *)calloc(mtot,sizeof(float));
	
	for (isubtype=1;isubtype<nsubtypes+1;isubtype++) {
		newlike = (float *)calloc(nsub,sizeof(float));
		prop = *(gp[isubtype].sigma);
		*(gp[isubtype].sigma) = (float)rnorm((double)(*(gp[isubtype].sigma)),(double)sd[isubtype],seed);
		if (*(gp[isubtype].sigma) > 0 ) {
			for (i=0;i<mtot;i++) 
				old[i] = gp[isubtype].Y.realp[i];
			
			yn = (*(gp[isubtype].sigma))/prop;
			//		vDSP_vsmul(gp[isubtype].Y.realp,1,&yn,gp[isubtype].Y.realp,1,mtot);
			for (i=0;i<mtot;i++)
				gp[isubtype].Y.realp[i] = yn*gp[isubtype].Y.realp[i];
			
			for (i=0;i<nsub;i++) {
				newlike[i] = 0;
				if (isubtype > 0) {
					if (sub[i].subtype == isubtype-1) 
						newlike[i] = likelihood_tst2(&(sub[i]),gp,cell_volume,ntot);
				}
				else
					newlike[i] = likelihood_tst2(&(sub[i]),gp,cell_volume,ntot);
			}
			new_like = 0;
			for (i=0;i<nsub;i++) {
				if (isubtype > 0) {
					if (sub[i].subtype == isubtype-1)
						new_like += newlike[i];
				}
				else 
					new_like += newlike[i];
			}
			
			like_ratio = new_like - *(gp[isubtype].likehd);
			
			prior_ratio = 0;// (ALPHA-1)*log(prop/(*(gp->sigma))) + BETA*(*(gp->sigma) - prop);
			
			//		 printf("isubtype = %d %f %f %f %f\n",isubtype,like_ratio,new_like,*(gp[isubtype].likehd),alpha);
			//		printf("\t %f %f\n",prop,*(gp[isubtype].sigma));
			if ((float)log(kiss(seed)) < like_ratio + prior_ratio) {
				if (isubtype>0) {
					*(gp[0].likehd) += like_ratio;
				}
				*(gp[isubtype].likehd) = new_like;                 // update current likelihood
				for (i=0;i<nsub;i++) {
					if (isubtype > 0) {
						if (sub[i].subtype == isubtype-1)
							*(sub[i].likehd) = newlike[i];
					}
					else {
						*(gp[sub[i].subtype+1].likehd) += newlike[i] - *(sub[i].likehd);
						*(sub[i].likehd) = newlike[i];
						
					}
				}
				accept[isubtype] += 1;
			}
			else {
				*(gp[isubtype].sigma) = prop;
				for (i=0;i<mtot;i++) 
					gp[isubtype].Y.realp[i] = old[i];
				
				//			yn = 1.0/yn;
				//			vDSP_vsmul(gp[isubtype].Y.realp,1,&yn,gp[isubtype].Y.realp,1,mtot);
			}
			attempt[isubtype] += 1;
		}
		else {
			*(gp[isubtype].sigma) = prop;
			attempt[isubtype] +=1;
		}
		free(newlike);
	}
	free(old);
}
*/

