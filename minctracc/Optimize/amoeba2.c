#include "volume_io.h"

#define NMAX 5000

#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0

#define SMALLEST  0.000001

#define GET_PSUM for (j=1;j<=ndim;j++) { for (i=1,sum=0.0;i<=mpts;i++)\
						sum += p[i][j]; psum[j]=sum;}


float amotry2(float **p,float *y,float *psum,int ndim,float (*funk)(),int ihi,int *nfunk,float fac);

int amoeba2(float **p, float y[], int ndim, float ftol, float (*funk)(), int *nfunk)
{
  int i,j,ilo,ihi,inhi,mpts=ndim+1;
  float ytry,ysave,sum,rtol,*psum;
  
  ALLOC(psum,ndim+1);
  *nfunk=0;
  GET_PSUM
    for (;;) {
      ilo=1;
      ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
      for (i=1;i<=mpts;i++) {
	if (y[i] < y[ilo]) ilo=i;
	if (y[i] > y[ihi]) {
	  inhi=ihi;
	  ihi=i;
	} else if (y[i] > y[inhi])
	  if (i != ihi) inhi=i;
      }
      rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
      
      if (rtol < ftol) break;
      
      /*lc*/		if (fabs(y[ilo]) < SMALLEST) break;
      
      if (*nfunk >= NMAX) return(FALSE);
      ytry=amotry2(p,y,psum,ndim,funk,ihi,nfunk,-ALPHA);
      if (ytry <= y[ilo])
	ytry=amotry2(p,y,psum,ndim,funk,ihi,nfunk,GAMMA);
      else if (ytry >= y[inhi]) {
	ysave=y[ihi];
	ytry=amotry2(p,y,psum,ndim,funk,ihi,nfunk,BETA);
	if (ytry >= ysave) {
	  for (i=1;i<=mpts;i++) {
	    if (i != ilo) {
	      for (j=1;j<=ndim;j++) {
		psum[j]=0.5*(p[i][j]+p[ilo][j]);
		p[i][j]=psum[j];
	      }
	      y[i]=(*funk)(psum);
	      ++(*nfunk);
	    }
	  }
	  *nfunk += ndim;
	  GET_PSUM
	  }
      }
    }
  FREE(psum);
  return(TRUE);
}

float amotry2(float **p,float *y,float *psum,int ndim,float (*funk)(),int ihi,int *nfunk,float fac)
{
  int j;
  float fac1,fac2,ytry,*ptry;
  
  ALLOC(ptry,ndim+1);
  fac1=(1.0-fac)/ndim;
  fac2=fac1-fac;
  for (j=1;j<=ndim;j++) 
    ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;

  ytry=(*funk)(ptry);
  ++(*nfunk);

  if (ytry < y[ihi]) {
    y[ihi]=ytry;
    for (j=1;j<=ndim;j++) {
      psum[j] += ptry[j]-p[ihi][j];
      p[ihi][j]=ptry[j];
    }
  }
  FREE(ptry);
  return ytry;
}

#undef ALPHA
#undef BETA
#undef GAMMA
#undef NMAX
