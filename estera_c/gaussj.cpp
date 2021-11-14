#include <math.h>
#define NRANSI
#include "nrutil.h"
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

int gaussj(double **a,int nzero, int n, double **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,temp;

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
                if(i<nzero)
                {
                   big=2.0;
                   irow=i;
                   icol=i;
                } else
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1)  {
                                                 free_ivector(ipiv,1,n);
        	                                 free_ivector(indxr,1,n);
                                        	 free_ivector(indxc,1,n);
                                                 return 1;
                                                }
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=nzero;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) {
                 free_ivector(ipiv,1,n);
        	 free_ivector(indxr,1,n);
        	 free_ivector(indxc,1,n);
                 return 1;
                }
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
       //         if(icol<nzero)
       //         {
       //         a[icol][icol] *= pivinv;
              //  b[icol][icol] *= pivinv;
       //         }
		for (l=1;l<=n;l++) if(a[icol][l]) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) if(b[icol][l]) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
                           //     if(icol<nzero)
                           //      {
                           //       a[icol][icol] -= a[icol][l]*dum;
                //                  b[icol][icol] -= b[icol][l]*dum;
                           //      }
                                if(dum)
                                {
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
                                }
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
        return 0;
}
#undef SWAP
#undef NRANSI
