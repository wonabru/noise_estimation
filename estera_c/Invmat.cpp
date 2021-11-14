#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

int multiply(double **MM,double ** M1,int mm1x,int mm1y,double ** M2, int mm2x)
{
	static int i,j,k;
for(i=0;i<mm1y;i++)
  for(j=0;j<mm2x;j++)
    MM[i][j]=0.0;

  for(i=0;i<mm1y;i++)
    for(j=0;j<mm2x;j++)
      for(k=0;k<mm1x;k++)
			MM[i][j]+=M1[i][k]*M2[k][j];

return 0;
}
int multiply1(double **MM,double ** M1,int mm1x,int mm1y,double ** M2, int mm2x)
{
	static int i,j,k;
for(i=1;i<=mm1y;i++)
  for(j=1;j<=mm2x;j++)
    MM[i][j]=0.0;

  for(i=1;i<=mm1y;i++)
    for(j=1;j<=mm2x;j++)
      for(k=1;k<=mm1x;k++)
			MM[i][j]+=M1[i][k]*M2[k][j];

return 0;
}



int inverseSVD(double **u,long n,double *w, double **v)
{
	static long i,j;
    static double **INV;
    INV=new double*[n+1];
    assert(INV);
    for(i=0;i<n+1;i++)
    	{
        	INV[i]=new double[n+1];
            assert(INV[i]);
        }

    for(j=1;j<=n;j++)
    {
    	for(i=1;i<=n;i++)
        	{
            	if(w[i])
        	        v[j][i]=v[j][i]/w[i];
                else v[j][i]=0.0;
            }
    }
    for(i=1;i<=n;i++)
    	for(j=1;j<=n;j++)
        	INV[i][j]=u[j][i];
    for(i=1;i<=n;i++)
    	for(j=1;j<=n;j++)
        	u[i][j]=INV[i][j];




    multiply1(INV,v,n,n,u,n);

    for(j=1;j<=n;j++)
    for(i=1;i<=n;i++)
    {
    	u[i][j]=INV[i][j];
    }

 for(i=0;i<n+1;i++)
 {
 	delete [] INV[i];
 }
 delete [] INV;
 return 0;
}
int addmatrix(float **S1,float **S2,long n, long m,float **S3)
{
	static long i,j;
	for(i=1;i<=n;i++)
    	for(j=1;j<=m;j++)
        {
        	S1[i][j]=S2[i][j]+S3[i][j];
        }

	return 0;
}
int minusMAT(float **S1,long n,long m)
{
	static long i,j;
	for(i=1;i<=n;i++)
    	for(j=1;j<=m;j++)
        	S1[i][j]=-S1[i][j];
    return 0;
}

int transpose(float **a,long n, long m,float **atr)
{
     static long i,j;

     for(i=1;i<=n;i++)
     	for(j=1;j<=m;j++)
        	atr[j][i]=a[i][j];

	return 0;
}


