/* note #undef's at end of file */
#define NRANSI
#include "nrutil.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

extern int spropen,Chbcn;
extern double conditioncut;
extern int svdcmp(double **a, int m, int n, double w[], double **v);
extern void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);

void mrqmin(double x[], double y[], double sig[], int ndata, double a[], int ia[],
	    int ma, double **covar, double **alpha, double *chisq,
	    void (*funcs)(double, double [], double *, double [], int), double *alamda)
{
  int inverseSVD(double **,long, double*,double**);
  void covsrt(double **covar, int ma, int ia[], int mfit);
  int gaussj(double **a,int nzero, int n, double **b, int m);
  void mrqcof(double x[], double y[], double sig[], int ndata, double a[],
	      int ia[], int ma, double **alpha, double beta[], double *chisq,
	      void (*funcs)(double, double [], double *, double [], int));
  int j,k,l;
  static int mfit,mfitalloc;
  static int allocate=0;
  static double ochisq,*atry,*beta,*da,**oneda;
    static double wmin,wmax,w[7],**v,B[7],BB[7];
    v=new double*[7];
    assert(v);
    for(j=0;j<7;j++)
      {
    	v[j]=new double[7];
        assert(v[j]);
      }
    if (*alamda < 0.0) {
      if(allocate==1)
	{
	  free_dmatrix(oneda,1,mfitalloc,1,1);
	  free_dvector(da,1,ma);
	  free_dvector(beta,1,ma);
	  free_dvector(atry,1,ma);
	  allocate=0;
	}
      atry=dvector(1,ma);
      beta=dvector(1,ma);
      da=dvector(1,ma);
      for (mfit=0,j=1;j<=ma;j++)
	if (ia[j]) mfit++;
      mfitalloc=mfit;
      oneda=dmatrix(1,mfitalloc,1,1);
      *alamda=0.001;
      mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs);
      ochisq=(*chisq);
      for (j=1;j<=ma;j++) atry[j]=a[j];
      allocate=1;
    }
    for (j=1;j<=mfit;j++) {
      for (k=1;k<=mfit;k++) covar[j][k]=alpha[j][k];
      covar[j][j]=alpha[j][j]*(1.0+(*alamda));
      oneda[j][1]=beta[j];
    }
    if(Chbcn==1)
      gaussj(covar,1,mfit,oneda,1);
    else
      {
	for(j=1;j<=mfit;j++)
	  B[j]=oneda[j][1];
	
	if(0==svdcmp(covar,mfit,mfit,w,v))
	  {
	    wmax=0;
	    wmin=w[1];
	    for(j=1;j<=mfit;j++) if(fabs(w[j]) > wmax&&w[j]!=0) wmax=fabs(w[j]);
	    
	    wmin=wmax*conditioncut;
	    
	    for(j=1;j<=mfit;j++)
	      if (fabs(w[j]) < wmin)
		w[j]=0.0;
	    
	    svbksb(covar,w,v,mfit,mfit,B,BB);
	    inverseSVD(covar,mfit,w,v);
	    
	    for(j=1;j<=mfit;j++)
	      {
		oneda[j][1]=BB[j];
	      }
	  }else
	    if(1==gaussj(covar,1,mfit,oneda,1))
	      {
		for(j=0;j<7;j++)
		  delete [] v[j];
		
		delete [] v;
		return;
	      }
      }
    //koniec dodatku
    
    for (j=1;j<=mfit;j++)
      {
	da[j]=oneda[j][1];
      }
    if (*alamda == 0.0) {
      covsrt(covar,ma,ia,mfit);
      covsrt(alpha,ma,ia,mfit);
      free_dmatrix(oneda,1,mfitalloc,1,1);
      free_dvector(da,1,ma);
      free_dvector(beta,1,ma);
      free_dvector(atry,1,ma);
      for(j=0;j<7;j++)
	delete [] v[j];
      delete [] v;
      allocate=0;
      return;
    }
    for (j=0,l=1;l<=ma;l++)
      if (ia[l])
        {
          atry[l]=a[l]+da[++j];
          if(atry[l]<0)
            atry[l]=a[l];
        }
    if(1<=x[ndata]*a[4])
      atry[4]=a[4];
    
    
    mrqcof(x,y,sig,ndata,atry,ia,ma,alpha,da,chisq,funcs);
    if (*chisq < ochisq) {
      if(*alamda>1e-10)
	*alamda *= 0.1;
      ochisq=(*chisq);
      for (j=1;j<=mfit;j++) {
	for (k=1;k<=mfit;k++) alpha[j][k]=covar[j][k];
	beta[j]=da[j];
      }
      for (l=1;l<=ma;l++)
	if(atry[l]>0&&ia[l]>0)
	  a[l]=atry[l];
    } else {
      *alamda *= 10;
      if(*alamda>1e5&&(*chisq)<(1.001*ochisq))
        {
	  ochisq=(*chisq);
	  for (j=1;j<=mfit;j++) {
	    for (k=1;k<=mfit;k++) alpha[j][k]=covar[j][k];
	    beta[j]=da[j];
	  }
	  for (l=1;l<=ma;l++)
	    if(atry[l]>0&&ia[l]>0)
	      a[l]=atry[l];
	  
        } else
	  *chisq=ochisq;
      
    }
    
    for(j=0;j<7;j++)
      delete [] v[j];
    
    delete [] v;
}
#undef NRANSI
