#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include "nrutil.h"
#include <fstream>

int delinear,ia[7],spropen=0,poczatek=0,Nofiles=1,fauto;
long Ndata,i,j,ii,ik,Nwin,shift,Nstart,Ndatagauss;
long double gz;
long double e2s,erf2;
double conditioncut,errorSIG;
int tau=1,idiffmax,d,tau2=3,kk,fittingiterations=10;
double *X,*XOLD,maxx,minx,epsmin=0,epsmax2,deps,beta1=0.564,mindS,maxdS,dSmin=0.0015;
double dS[200],K[200],E[200],dS1[200],powereps=1,sig[200],**alfa,**covar,mean,stddev,dziel=1;
double a[7],omega,chisq,chisq1,alamda,chisqold;
double powerepsmin=0.5,powerepsmax=7,power1,dpower1;
double omegapow=0.55,stddevofcur,daa[7];
char sinputfile[1000],soutputfile[1000];
int  Automatic=1,shufflepar=0,Chbmultifitting=1,Chbcontrol=1,powernumber=1;
int epsnumber;
int Chbpoprawka=1,Chbcn=0,Chbsig=1,ChbK=1,Chba=1,Chbb=1,Chbc=1,Chbd=1,gausssurr=1;
const float Pi=3.14159;

using namespace std;

void mrqmin(double x[], double y[], double sig[], int ndata, double a[], int ia[],
	int ma, double **covar, double **alpha, double *chisq,
	void (*funcs)(double, double [], double *, double [], int), double *alamda);
void four1(double data[], unsigned long nn, int isign);


//------------------------------------------------------

long random01(long numkil)
{
    return (long)(numkil*(rand()-3)*1.0/RAND_MAX);
}
double getmaxdS()
{
    double dSS[200],dSSS[200],pom3;
//        for(i=1;i<=kk;i++)
//                dS[i]/=pow(E[i],powereps);

    dSSS[1]=dS[1];
    for(i=2;i<kk;i++)
    {
        dSSS[i]=(dS[i-1]+dS[i]+dS[i+1])/3.0;
    }
    dSSS[kk]=dS[kk];

    dSS[1]=dSSS[1];
    for(i=2;i<kk;i++)
    {
        dSS[i]=(dSSS[i-1]+dSSS[i]+dSSS[i+1])/3.0;
    }
    dSS[kk]=dSSS[kk];
    pom3=0;
    int pom4=0;
    double maxdss=0;
    for(i=2;i<=0.33*kk;i++)
    {
        if(dSS[i]<dSS[i+1])
                pom4=1;
//        if(dSS[i]>dSS[i+3])
//        if(pom4==1&&maxdss<dSS[i])
        if(dSS[i]>dSS[i+1]&&pom4==1&&maxdss<dSS[i])
        {
              maxdss=dSS[i];
              pom3=E[i];
              pom4=0;
        }
    }
   return pom3;
}

//------------------------------------------------

int calcdS1()
{
   if(powernumber)
     dpower1=(powerepsmax-powerepsmin)/powernumber;

   for(i=1;i<=kk;i++)
   {
        dS1[i]=0;
   for(power1=powerepsmin;power1<=powerepsmax;power1+=dpower1)
   {
        omega=pow(omegapow,-power1);
       dS1[i]+=omega*K[i]*pow((double)E[i],(double)power1);
   }
   }
   return 0;
}

//--------------------------------------------------------------

double getmaxdS1(double *max1,double *max2)
{
    double dSS[200],dSSS[200],pom3;
//        for(i=1;i<=kk;i++)
//                dS[i]/=pow(E[i],powereps);

    calcdS1();
    dSSS[1]=dS1[1];
    for(i=2;i<kk;i++)
    {
        dSSS[i]=(dS1[i-1]+dS1[i]+dS1[i+1])/3.0;
    }
    dSSS[kk]=dS1[kk];

    dSS[1]=dSSS[1];
    for(i=2;i<kk;i++)
    {
        dSS[i]=(dSSS[i-1]+dSSS[i]+dSSS[i+1])/3.0;
    }
    dSS[kk]=dSSS[kk];

    pom3=0;
    *max1=1;
    for(i=1;i<=0.43*kk;i++)
    {

        if(pom3<dSS[i])
        {
              pom3=dSS[i];
              *max1=dSS[i];
        }
    }
    pom3=0;
    *max2=1;
    for(i=(int)(kk/2.0);i<=kk-1;i++)
    {
        if(dSS[i]>pom3)
        {
              pom3=dSS[i];
              *max2=dSS[i];

        }
    }
   return 0;
}
//-----------------------------------------------------

void searchendline(ifstream *filein)
{
   char c[1000];
   filein->getline(c,1000);
}

//-------------------------------------------------------

double gauss(double sd)
{
   double x;
   x=0;
   for(j=0;j<12;j++)
       x+=rand()*1.0/RAND_MAX;
   x-=6;
   x*=sd;
   return x;
}

//---------------------------------------------------

double thetaS(double *X,double eps,double *S,double beta1,int norm)
{
double dS;
double theta,srtheta3,theta3;
double pmax;
long iik,jk,p,iii;

double mnoznik;
mnoznik=(1.0-beta1)*tau*eps;
theta3=0;

if(tau2<=0)
        tau2=1;

pmax=0;
iii=0;
for(iik=0;iik<Ndata;iik++)
   {
   theta=0;
   p=0;
for(jk=fauto;jk+iik<Ndata;jk++)
 {
        theta=0;
      	for(p=0;jk+iik+p*tau<Ndata;p++)
        {
            theta+=fabsl(X[iik+p*tau]-X[iik+jk+p*tau]);
            if(theta>=mnoznik*(p+1)) break;
        }

        if(jk+iik+p*tau>=Ndata)
        {
                theta=0;
                p=0;
        }
        if(p>=tau2)
        {
           theta3+=p;
           iii++;
        }
     }
   }
   if(iii>0)
        srtheta3=theta3*1.0/iii;
   else srtheta3=0;

   if(srtheta3>tau2&&(srtheta3-tau2))
   {
     dS=logl((srtheta3-tau2+1)/(srtheta3-tau2));
    } else dS=0;

  return dS;
}

//--------------------------------------------------

void Calculateentropy()
{
   static double pom,epsmax;
   spropen=2;
   epsmax=stddev*epsmax2;
  /* if(Automatic==1)
     {
       epsmax=stddev*2.4;
       if(epsnumber)
	 deps=(epsmax-epsmin)/epsnumber;
       else deps=1;
       tau2=4;
       tau=1;
       epsmin=0;
       beta1=0.564;
       dSmin=0.015;
     }
  */
    deps=2*(epsmax-epsmin)/epsnumber;

   double sdszumu=0;
   kk=1;
   for(double eps=epsmax*0.5;eps<=epsmax;eps+=deps)
   {
 //    cout <<"Norm:"<<norm<<" eps="<<eps<<endl;
     K[kk]=thetaS(X,eps,&pom,beta1,1);

     if(K[kk]==0) break;
     if(kk>=4)
     {
       if(pow(eps,0.9)*K[kk]>pow((eps-2*deps),0.9)*K[kk-2])
       {
          epsmax=eps;
          break;
       }
     }
     kk++;
   }
    if(epsnumber)
      deps=(epsmax-epsmin)/epsnumber;
    else deps=1;
   int bylsd=0;
   kk=1;
   for(double eps=epsmin;eps<=epsmax;eps+=deps)
   {
     E[kk]=eps;
     K[kk]=thetaS(X,eps,&pom,beta1,1);

     if(K[kk]>0.9) sdszumu=E[kk];
        if(bylsd<1&&sdszumu!=eps&&kk>=10)
        {
            beta1=0.59*sdszumu/stddev*sdszumu/stddev+0.20;
            if(beta1<0.44) beta1=0.44;
            if(beta1>0.65) beta1=0.65;
            eps=epsmin;
            kk=0;
            bylsd++;
       cout <<"beta="<<beta1<<endl;
     }
    if(kk>=epsnumber*0.5)
     {
       if(pow(eps,1.1)*K[kk]>pow(eps-5*deps,1.1)*K[kk-5])
       {
          kk-=5;
          break;
       }
     }
     if(K[kk]<dSmin&&K[kk]!=0)
     {
        kk++;
    	break;
     }
     if(K[kk]!=0&&K[kk]<3.1||kk!=1)
    	 kk++;
 }
 kk--;

 cout <<"Forecast SD of noise = "<< sdszumu<<endl;
 dziel=E[kk]/0.7;
 if(dziel)
   for(i=1;i<=kk+1;i++)
     if(i<epsnumber)
     	 E[i]/=dziel;

 for(i=1;i<=kk;i++)
    dS[i]=K[i]*pow(E[i],powereps);

 maxdS=dS[1];
 mindS=dS[1];
 for(i=1;i<=kk;i++)
 {
    if(maxdS<dS[i]) maxdS=dS[i];
    if(mindS>dS[i]) mindS=dS[i];
 }
 poczatek=1;
}

//--------------------------------------------------

int repowereps(double pom2)
{
      powereps=pom2;
      for(i=1;i<=kk;i++)
          dS[i]=K[i]*pow(E[i],powereps);

      maxdS=dS[1];
      mindS=dS[1];
      for(i=1;i<=kk;i++)
      {
	  if(maxdS<dS[i]) maxdS=dS[i];
          if(mindS>dS[i]) mindS=dS[i];
      }
      return 0;
}

//----------------------------------------------------

double erf(double z)
{
   static double F,t1,Fcal;

   if(z<=0)
   	return 0;
   static double dt,t;
   if(z>10)
   	  z=10;
   dt=1.0e-3;
   Fcal=0;
   for(t=0;t<=z;t+=dt)
   {
   F=exp(-t*t);
   t1=t+0.5*F*dt;
   F+=2*exp(-t1*t1);
   t1=t+0.5*exp(-t1*t1)*dt;
   F+=2*exp(-t1*t1);
   t1=t+exp(-t1*t1)*dt;
   F+=exp(-t1*t1);
   Fcal+=dt/6.0*F;
   }
   return Fcal;
}

//--------------------------------------------------------------------------

void entrope(double eps,double a[], double *yfit, double dyda[], int ma)
{
    if(a[1]<0.001)
        a[1]=0.001;
    if(a[1]>stddev/dziel)
        a[1]=stddev/dziel;

//    for(j=2;j<=6;j++)
//      if(a[j]>100)
 //        a[j]=100;
  //    else if(a[j]<-100)
   //      a[j]=-100;

    if(a[1])
    e2s=eps/(2.0*a[1]);
    else e2s=0;

    erf2=erf(e2s);
 //   if(erf2==0)
 //   	erf2=1;
   if(a[4]*eps>=1&&eps)
        a[4]=0.96/eps;

    long double aa,poweps;
    aa=a[1]*a[1];
    if(erf2)
      gz=2.0/sqrt(Pi)*e2s*exp(-e2s*e2s)/erf2;
    else gz=0;
   long double exp2aa;
   if(aa&&(eps*eps/(2.*aa))>1e-6)
     exp2aa=expl(-eps*eps/(2.*aa));
   else exp2aa=0;
   long double log2eps;
   if(eps>0)
    log2eps=logl(eps);
   else log2eps=0;

   *yfit=0;
   for(i=1;i<=6;i++)
        dyda[i]=0;
  if(Chbpoprawka==0)
   {
   if(Chbmultifitting==1)
   {
    if(powernumber)
   dpower1=(powerepsmax-powerepsmin)/powernumber;
   else dpower1=0;
   for(power1=powerepsmin;power1<=powerepsmax;power1+=dpower1)
   {
    poweps=pow(eps,power1);
    omega=pow(omegapow,-power1);
    if(eps>0&&(1-a[4]*eps)>0)
      *yfit+=omega*((a[2]+a[3]*logl(1-a[4]*eps))*poweps*(1+a[5]*a[1]/eps)-a[6]*gz*poweps*logl(eps));

   if(*yfit>1e2*powernumber)
        *yfit=1e2*powernumber;
   if(*yfit<-10)
        *yfit=-10;
   if((1-eps*a[4])>0&&eps&&aa*aa*erf2*erf2&&exp2aa>=0)
    dyda[1]+=omega*(-((2.0*poweps*eps*eps*a[1]*a[6]*log2eps)*exp2aa+(eps*poweps*sqrt(Pi)*(eps*eps - 2*aa)*a[6]*erf2*log2eps)*sqrtl(exp2aa))/(2.*Pi*aa*aa*erf2*erf2)+poweps/eps*a[5]*(a[2]+a[3]*logl(1-eps*a[4])));
   if(eps)
   dyda[2]+=omega*(poweps*(1+a[1]*a[5]/eps));
   if(eps&&(1.0-eps*a[4])>0)
   dyda[3]+=omega*(poweps*(1+a[1]*a[5]/eps)*logl(1.0-eps*a[4]));
   if(eps&&(-1 + eps*a[4]))
   dyda[4]+=omega*((eps*a[3]*poweps*(1.0 + a[1]*a[5]/eps))/(-1 + eps*a[4]));
   if(eps&&(1.0-eps*a[4])>0)
   dyda[5]+=omega*(poweps/eps*a[1]*(a[2] + a[3]*logl(1.0-eps*a[4])));
   if(eps>0)
   dyda[6]+=-omega*gz*poweps*logl(eps);

   for(i=1;i<=6;i++)
   if(dyda[i]>1e8)
        dyda[i]=1e8;
   for(i=1;i<=6;i++)
   if(dyda[i]<-1e8)
        dyda[i]=-1e8;

   }
   }else
   {
    poweps=pow((double)eps,(double)powereps);
  if(eps>0&&(1-a[4]*eps)>0)
    *yfit=(a[2]+a[3]*logl(1-a[4]*eps))*poweps*(1+a[5]*a[1]/eps)-a[6]*gz*poweps*logl(eps);
   else
      *yfit=0;

   if(*yfit<-10)
        *yfit=-10;
   if(*yfit>10)
        *yfit=10;

   if((1-eps*a[4])>0&&eps&&aa*aa*erf2*erf2&&exp2aa>=0)
   dyda[1]=-((2.0*poweps*eps*eps*a[1]*a[6]*log2eps)*exp2aa+(eps*poweps*sqrt(Pi)*(eps*eps - 2*aa)*a[6]*erf2*log2eps)*sqrtl(exp2aa))/(2.*Pi*aa*aa*erf2*erf2)+poweps/eps*a[5]*(a[2]+a[3]*logl(1-eps*a[4]));
   if(eps)
   dyda[2]=poweps*(1+a[1]*a[5]/eps);
   if(eps&&(1.0-eps*a[4])>0)
   dyda[3]=poweps*(1+a[1]*a[5]/eps)*logl(1.0-eps*a[4]);
   if(eps&&(-1 + eps*a[4]))
   dyda[4]=(eps*a[3]*poweps*(1.0 + a[1]*a[5]/eps))/(-1 + eps*a[4]);
   if(eps&&(1.0-eps*a[4])>0)
   dyda[5]=poweps/eps*a[1]*(a[2] + a[3]*logl(1.0-eps*a[4]));
   if(eps>0)
   dyda[6]=-gz*poweps*logl(eps);
   }
   }else
   {
   double stddevd=stddev/dziel;
   if(stddevd==0) return;
   if(Chbmultifitting==1)
   {
    if(powernumber)
   dpower1=(powerepsmax-powerepsmin)/powernumber;
   else dpower1=0;
   for(power1=powerepsmin;power1<=powerepsmax;power1+=dpower1)
   {
    poweps=pow(eps,power1);
    omega=pow(omegapow,-power1);
    if(eps>0&&(1-a[4]*eps)>0)
      *yfit+=omega*((a[2]+a[3]*logl(1-a[4]*eps))*poweps*(1+a[5]*a[1]/eps)-(1-a[1]/stddevd)*a[6]*gz*poweps*logl(eps));

   if(*yfit>1e2*powernumber)
        *yfit=1e2*powernumber;
   if(*yfit<-10)
        *yfit=-10;
   if((1-eps*a[4])>0&&eps&&aa*aa*erf2*erf2&&exp2aa>=0)
    dyda[1]+=omega*(1.0/stddevd*a[6]*gz*poweps*logl(eps)-((2.0*poweps*eps*eps*a[1]*a[6]*log2eps)*exp2aa+(eps*poweps*sqrt(Pi)*(eps*eps - 2*aa)*a[6]*erf2*log2eps)*sqrtl(exp2aa))/(2.*Pi*aa*aa*erf2*erf2)*(1-a[1]/stddevd)+poweps/eps*a[5]*(a[2]+a[3]*logl(1-eps*a[4])));
   if(eps)
   dyda[2]+=omega*(poweps*(1+a[1]*a[5]/eps));
   if(eps&&(1.0-eps*a[4])>0)
   dyda[3]+=omega*(poweps*(1+a[1]*a[5]/eps)*logl(1.0-eps*a[4]));
   if(eps&&(-1 + eps*a[4]))
   dyda[4]+=omega*((eps*a[3]*poweps*(1.0 + a[1]*a[5]/eps))/(-1 + eps*a[4]));
   if(eps&&(1.0-eps*a[4])>0)
   dyda[5]+=omega*(poweps/eps*a[1]*(a[2] + a[3]*logl(1.0-eps*a[4])));
   if(eps>0)
   dyda[6]+=-omega*(1-a[1]/stddevd)*gz*poweps*logl(eps);

   for(i=1;i<=6;i++)
   if(dyda[i]>1e8)
        dyda[i]=1e8;
   for(i=1;i<=6;i++)
   if(dyda[i]<-1e8)
        dyda[i]=-1e8;

   }
   }else
   {
    poweps=pow((double)eps,(double)powereps);
  if(eps>0&&(1-a[4]*eps)>0)
    *yfit=(a[2]+a[3]*logl(1-a[4]*eps))*poweps*(1+a[5]*a[1]/eps)-(1-a[1]/stddevd)*a[6]*gz*poweps*logl(eps);
   else
      *yfit=0;

   if(*yfit<-10)
        *yfit=-10;
   if(*yfit>10)
        *yfit=10;

   if((1-eps*a[4])>0&&eps&&aa*aa*erf2*erf2&&exp2aa>=0)
   dyda[1]=-1.0/stddevd*a[6]*gz*poweps*logl(eps)-((2.0*poweps*eps*eps*a[1]*a[6]*log2eps)*exp2aa+(eps*poweps*sqrt(Pi)*(eps*eps - 2*aa)*a[6]*erf2*log2eps)*sqrtl(exp2aa))/(2.*Pi*aa*aa*erf2*erf2)*(1-a[1]/stddevd)+poweps/eps*a[5]*(a[2]+a[3]*logl(1-eps*a[4]));
   if(eps)
   dyda[2]=poweps*(1+a[1]*a[5]/eps);
   if(eps&&(1.0-eps*a[4])>0)
   dyda[3]=poweps*(1+a[1]*a[5]/eps)*logl(1.0-eps*a[4]);
   if(eps&&(-1 + eps*a[4]))
   dyda[4]=(eps*a[3]*poweps*(1.0 + a[1]*a[5]/eps))/(-1 + eps*a[4]);
   if(eps&&(1.0-eps*a[4])>0)
   dyda[5]=poweps/eps*a[1]*(a[2] + a[3]*logl(1.0-eps*a[4]));
   if(eps>0)
   dyda[6]=-gz*(1-a[1]/stddevd)*poweps*logl(eps);
   }

   for(i=1;i<=6;i++)
   if(dyda[i]>1e8)
        dyda[i]=1e8;
   for(i=1;i<=6;i++)
   if(dyda[i]<-1e8)
        dyda[i]=-1e8;

   }
}

//--------------------------------------------------------------------

void Fitting1()
{
  long double sumcovar;

  for(i=1;i<=kk;i++)
    sig[i]=1.0;
  if(spropen==2)
    {
      a[1]=E[1];
      for(i=1;i<=kk;i++)
        if(K[i]>0.9&&E[i]<stddev/dziel)
	  a[1]=E[i];
      if(a[1]==E[1])
      {
          double eps1;
	  do
	    {
	      eps1=getmaxdS();
	      if(eps1==0)
		     repowereps(powereps*0.96);
	      if(powereps<0.3441717437358231)
		    break;
	    }while(eps1==0);
          a[1]=eps1;
      }

alamda=-1;
chisq=100000;
a[4]=1.3;
a[3]=0.1;
a[6]=0.7;
a[5]=1;
a[2]=0.1;
if(a[1]>0.24)
        ia[4]=0;
ia[1]=0;
ia[2]=0;
ia[3]=0;
ia[4]=0;
ia[5]=0;
ii=0;
    }
  if(spropen==2&&Chbcontrol==1||spropen==2&&Automatic==1)
    {
      int lk;
      double eps1,maxeps;
      lk=0;
      repowereps(1.4);
      do
	{
	  do
	    {
	      eps1=getmaxdS();
	      if(eps1==0)
		repowereps(powereps*0.96);
	      if(powereps<0.3441717437358231)
		break;
	    }while(eps1==0);
            maxeps=eps1;
            if(fabs(a[1]-maxeps)>a[1]*1e-2&&maxeps>0)
             a[1]=maxeps;

          if(a[1]>stddev/dziel)
            a[1]=stddev/dziel;
    if(maxeps>1e-6&&maxeps!=1)
    repowereps(0.3441717437358231-1.0/logl(maxeps));
     lk++;
	}while(lk<20);

    if(a[1]>1e-6&&a[1]!=1)
    {
   repowereps(0.3441717437358231-(double)(1.0/logl(a[1])));
    }
    }
  if(Automatic==1)
    {
      Chbmultifitting=0;
      if(a[1]>0.24)
	{
	  Chbmultifitting=1;
	  omegapow=0.55;
	  powerepsmax=7;
	  powerepsmin=0.5;
	  int pom10=-1;
	  int lk=0;
	  double max1=0,max2=0;
	  do
	    {
	      getmaxdS1(&max1,&max2);
	      if(max1>max2)
		{
		  if(pom10==0) break;
		  pom10=1;
		  omegapow*=0.96;
		}
	      if(max1<max2)
		{
		  if(pom10==1) break;
		  pom10=0;
		  omegapow*=1.03;
		}
	      calcdS1();
	      lk++;
	    }while(lk<20);
	}
    }
  powernumber=1;

  if(Chbmultifitting==1)
    {
      if(powernumber)
	dpower1=(powerepsmax-powerepsmin)/powernumber;
      else dpower1=0;

      for(i=1;i<=kk;i++)
	{
	  dS1[i]=0;
	  for(power1=powerepsmin;power1<=powerepsmax;power1+=dpower1)
	    {
	      omega=pow(omegapow,-power1);
	      dS1[i]+=omega*K[i]*pow((double)E[i],(double)power1);
	    }
	}
    }




  if(spropen==5)
    {
      alamda=0;
    }

  ik=0;

  do
    {
      if(Chbmultifitting==1)
	{
	  if(powernumber)
	    dpower1=(powerepsmax-powerepsmin)/powernumber;
	  else
	    dpower1=0;

	  for(i=1;i<=kk;i++)
	    {
	      dS1[i]=0;
	      for(power1=powerepsmin;power1<=powerepsmax;power1+=dpower1)
		{
		  omega=pow(omegapow,-power1);
		  dS1[i]+=omega*K[i]*pow((double)E[i],(double)power1);
		}
	    }
	}
      if(spropen<5&&Chbcontrol==1)
	{
	  int lk;
	  double eps1,maxeps;
	  lk=0;
      repowereps(1.4);
      do
	{
	  do
	    {
	      eps1=getmaxdS();
	      if(eps1==0)
		repowereps(powereps*0.96);
	      if(powereps<0.3441717437358231)
		break;
	    }while(eps1==0);
    maxeps=eps1;

    if(fabs(a[1]-maxeps)>a[1]*1e-2&&a[1]>0.1&&maxeps>0.1)
            a[1]=maxeps;
     else break;
//    else break;
    if(a[1]>stddev/dziel)
     a[1]=stddev/dziel;
    if(a[1]>1e-6&&a[1]!=1)
    repowereps(0.3441717437358231-1.0/logl(a[1]));
    lk++;
    }while(lk<4);
    if(a[1]>1e-6&&a[1]!=1)
    repowereps(0.3441717437358231-1.0/logl(a[1]));
      }
      if(spropen<5)
      if(Automatic==1)
	{
    if(ii==0)
        {
        alamda=-1;
        }

    if(ii==5)
    {
        alamda=0;
    }
    if(ii==15)
    {
        alamda=0;
    }
    if(ii==20)
    {
        alamda=0;
    }
  if(ii==3)
  {
        alamda=0;
  }
	}else{
	  if(Chbsig==1)
	    {
	      if(ia[1]==0) alamda=0;
	    }
	  else
	    {
	      if(ia[1]==1) alamda=0;
	    }
	  if(ChbK==1)
	    {
	      if(ia[2]==0) alamda=0;
	    }
	  else
	    {
	      if(ia[2]==1) alamda=0;
	    }
	  if(Chba==1)
	    {
	      if(ia[4]==0) alamda=0;
	    }
	  else
	    {
	      if(ia[4]==1) alamda=0;
	    }
	  if(Chbb==1)
	    {
	      if(ia[3]==0) alamda=0;
	    }
	  else
	    {
	      if(ia[3]==1) alamda=0;
	    }
	  if(Chbc==1)
	    {
	      if(ia[6]==0) alamda=0;
	    }
	  else
	    {
	      if(ia[6]==1) alamda=0;
	    }
	  if(Chbd==1)
	    {
	      if(ia[5]==0) alamda=0;
	    }
	  else
	    {
	      if(ia[5]==1) alamda=0;
	    }
	}
    chisqold=chisq;

if(poczatek==1) alamda=-1;

if(alamda==0||ik==1&&spropen==5)
if(Chbmultifitting==0)
	mrqmin(E,dS,sig,kk,a,ia,6,covar,alfa,&chisq,*entrope,&alamda);
else
	mrqmin(E,dS1,sig,kk,a,ia,6,covar,alfa,&chisq,*entrope,&alamda);

if(spropen==5&&alamda==0)
{
        ia[1]=1;
        for(i=2;i<=6;i++)
           ia[i]=1;
}

if(alamda==0) alamda=-1;

if(ik==1&&spropen==5)
    alamda=0;


poczatek=0;
      if(spropen<5)
	{
	  if(Automatic==1)
	    {
    if(ii==0)
        {
        Chbcontrol=1;
        ia[1]=0;
        ia[5]=0;
        ia[6]=1;
        }

    if(ii==5)
    {
	  ia[5]=0;
    	ia[1]=1;
    }
    if(ii==15)
    {
        ia[4]=1;
	  ia[5]=1;
    	ia[1]=0;
    }
    if(ii==20)
    {
        ia[2]=1;
        ia[3]=1;
        ia[4]=1;
        ia[5]=1;
        ia[6]=1;
    	ia[1]=1;
            Chbcontrol=0;

            }
          if(ii==3)
          {
          	ia[5]=1;
                ia[6]=1;
          }
	    }else{
	      if(Chbsig==1)
		{
		  ia[1]=1;
		}
	      else
		{
		  ia[1]=0;
		}
	      if(ChbK==1)
		{
		  ia[2]=1;
		}
	      else
		{
		  ia[2]=0;
		}
	      if(Chba==1)
		{
		  ia[4]=1;
		}
	      else
		{
		  ia[4]=0;
		}
	      if(Chbb==1)
		{
		  ia[3]=1;
		}
	      else
		{
		  ia[3]=0;
		}
	      if(Chbc==1)
		{
		  ia[6]=1;
		}
	      else
		{
		  ia[6]=0;
		}
	      if(Chbd==1)
		{
		  ia[5]=1;
		}
	      else
		{
		  ia[5]=0;
		}
	    }
	}

  if(a[1]>0.24&&Automatic==1)
        ia[4]=0;


if(Chbmultifitting==0)
	mrqmin(E,dS,sig,kk,a,ia,6,covar,alfa,&chisq,*entrope,&alamda);
else
	mrqmin(E,dS1,sig,kk,a,ia,6,covar,alfa,&chisq,*entrope,&alamda);

if(ik==1&&spropen==5)
    poczatek=1;

       mean=0;
       for(i=0;i<kk;i++)
          mean+=E[i];
       if(kk)
       mean/=kk;
       else mean=0;
       if(kk&&dziel&&mean>1e-8&&chisq>=0)
          errorSIG=sqrt(chisq/kk)*0.75/dziel*stddev/mean;

      if(alamda>10e15)
	{

	  alamda/=10.0;
	  if(chisq==chisqold)
	    {
	      break;
	    }

	}
      if(alamda==0) alamda=-1;
      if(kk)
	chisq1=chisq1/kk;
      ii++;
      ik++;
    }while(ik<fittingiterations&&spropen<5||spropen==5&&ik<2);
  spropen=3;
}

//---------------------------------
long randomones(int reset,long num)
{
   static long *liczbyl,byln=0,num2;
   if(reset==0)
   {
      liczbyl=new long[num];
      for(long i=0;i<num;i++)
           liczbyl[i]=0;
       byln=0;
   }
   if(reset==1)
   {
      if(num-byln>=1)
      {
      num2=random01(num-byln);
      for(long i=0;i<num;i++)
      {
         if(liczbyl[i]==0)
         {
           if(num2==0)
           {
               liczbyl[i]=1;
               byln++;
               return i;
           }
           num2--;
         }
      }
      } else return -1;
   }

   if(reset==2)
   {
      delete [] liczbyl;
   }
   return 0;
}

//--------------------------------

int likwidendl(char *fileins)
{
    ifstream filein(fileins);
    filein.clear();
    assert(filein);
    char c,c1;
    int mk=0;
    do
    {
    if(mk!=0)
      filein.unget();
    c=filein.get();
    filein.get();
    mk++;
    }while(filein.eof()==0);
      filein.close();
      char *st;
      st=new char[mk];
      assert(st);
      filein.open(fileins);
      filein.clear();
      for(int klo=0;klo<mk;klo++)
      {
         st[klo]=filein.get();
      }
      filein.close();
      ofstream fileout(fileins);
      assert(fileout);
      if(c=='\n')
      {
         mk--;
      }
      for(int klo=0;klo<mk;klo++)
      {
        fileout << st[klo];
      }
      fileout.close();
      delete [] st;
      return 0;
}
//----------------------------------------------------------------

void Open(char *sin,long *Ndatao)
{
  ifstream filein;

  filein.open((const char *)sin);
  assert(filein);
  filein.width(16);
  filein.clear();

  for(i=0;i<Nstart&&!filein.eof();i++)
      searchendline(&filein);


  for(i=0;i<*Ndatao&&!filein.eof();i++)
    {
      filein >> X[i];
      if(filein.eof())
	break;
    }
  *Ndatao=i;
  filein.close();
}

//------------------------------------------------

int shuffle(int manyl)
{
       double pom;
       long ist,jst;

   for(int fl=0;fl<manyl;fl++)
   {
       randomones(0,Ndata);
       for(i=0;i<Ndata;i++)
       {
           ist=randomones(1,Ndata);
           jst=randomones(1,Ndata);
           if(ist==-1||jst==-1) break;
           pom=X[ist];
           X[ist]=X[jst];
           X[jst]=pom;
       }
       randomones(2,Ndata);
   }
   return 0;
}
//------------------------------------------

int delinearization(long Ndet)
{
   double *XR,sdmin=0;
   long idiffmin=1;
   XR=new double[Ndet];
   assert(XR);
   Ndet--;
   int ddiff;
   for(ddiff=3;ddiff<=idiffmax;ddiff++)
   {
       for(i=0;i<Ndet+1;i++)
          XR[i]=X[i];
   for(j=0;j<ddiff;j++)
   {
    double means=0;
    long meansco=0;
    for(i=j;i<Ndet;i+=ddiff)
    {
        means+=(XR[i]+XR[i+1])/2.0;
        meansco++;
    }
    if(meansco)
      means/=meansco;
    for(i=j;i<Ndet;i+=ddiff)
    {
       XR[i]-=means;
    }
   }
    double mean1=0;
    double sd1=0;
    for(i=0;i<Ndet;i++)
    {
       mean1+=XR[i];
       sd1+=XR[i]*XR[i];
    }
    if(Ndet)
    {
       mean1/=Ndet;
       sd1/=Ndet;
    }
    if(sd1-mean1*mean1>0)
      sd1=sqrt(sd1-mean1*mean1);
    else sd1=0;
    if(sd1<sdmin||sdmin==0)
    {
       idiffmin=ddiff;
       sdmin=sd1;
    }
   }
   for(j=0;j<idiffmin;j++)
   {
    double means=0;
    long meansco=0;
    for(i=j;i<Ndet;i+=idiffmin)
    {
        means+=(X[i]+X[i+1])/2.0;
        meansco++;
    }
    if(meansco)
      means/=meansco;
    for(i=j;i<Ndet;i+=idiffmin)
    {
       X[i]-=means;
    }
   }
   delete [] XR;
   return 0;
}

//------------------------------------------------
void Gaussiansurrogates(long *Ndatalok)
{
  double *xx,*xxold;
  long *nn,*nnold,k,j,s,Ndatag;
  double pom;
  Ndatag=Ndatagauss;
//  Ndatag--;
//  Ndatag=Ndatag-Ndatag%Ndata;


  if(Ndatag<Ndata) {
     Ndatag=Ndata;
  }
  Open(sinputfile,&Ndatag);
  (*Ndatalok)=Ndatag;
  xx=new double[Ndatag];
  assert(xx);
  xxold=new double[Ndatag];
  assert(xxold);
  nn=new long[Ndatag];
  assert(nn);
  nnold=new long[Ndatag];
  assert(nnold);

  if(delinear==1) delinearization(Ndatag);

  mean=0;
  for(i=0;i<Ndatag;i++)
  {
      mean+=X[i];
  }
  if(i)
    mean/=i;
  else mean=0;
  stddev=0;
  for(i=0;i<Ndatag;i++)
  {
     stddev+=(X[i]-mean)*(X[i]-mean);
  }
  if(i)
    stddev/=i;
  else stddev=0;
  if(stddev>=0)
    stddev=sqrt(stddev);


  for(j=0;j<Ndatag;j++)
  {
      xxold[j]=X[j];
      nnold[j]=j;
  }
  for(j=0;j<Ndatag;j++)
    {
      if(gausssurr==2)
      {
        do
        {
          xx[j]=gauss(stddev);
        }while(fabs(xx[j])>6*stddev);
      } else xx[j]=j*3.46*stddev/Ndatag;
      nn[j]=j;
    }

  if(gausssurr==2)
  for(k=0;k<Ndatag;k++)
    for(j=0;j<Ndatag-1;j++)
      if(xx[j]>xx[j+1])
	{
	  pom=xx[j];
	  xx[j]=xx[j+1];
	  xx[j+1]=pom;
	  pom=nn[j];
	  nn[j]=nn[j+1];
	  nn[j+1]=(long)pom;
	}

  for(k=0;k<Ndatag;k++)
    for(j=0;j<Ndatag-1;j++)
      if(xxold[j]>=xxold[j+1])
	{
	  pom=xxold[j];
	  xxold[j]=xxold[j+1];
	  xxold[j+1]=pom;
	  pom=nnold[j];
	  nnold[j]=nnold[j+1];
	  nnold[j+1]=(long)pom;
	}

  for(k=0;k<Ndatag;k++)
    for(j=0;j<Ndatag;j++)
      {
	s=random01(Ndatag);
	if(xxold[s]==xxold[j])
	  {
	   pom=xxold[j];
	   xxold[j]=xxold[s];
	   xxold[s]=pom;
	   pom=nnold[j];
	   nnold[j]=nnold[s];
	   nnold[s]=(long)pom;
	  }
      }

  for(j=0;j<Ndatag;j++)
    X[nnold[j]]=xx[j];

  if(Ndatag<Ndata) Ndata=Ndatag;

  mean=0;
  for(i=0;i<Ndata;i++)
    {
      mean+=X[i];
     }
  if(i)
    mean/=i;
  else mean=0;
   stddev=0;
   for(i=0;i<Ndata;i++)
     {
       stddev+=(X[i]-mean)*(X[i]-mean);
     }
   if(i)
     stddev/=i;
   else stddev=0;
   if(stddev>=0)
     stddev=sqrt(stddev);
   maxx=minx=X[0];
   for(i=0;i<Ndata;i++)
     {
       if(maxx<X[i]) maxx=X[i];
       if(minx>X[i]) minx=X[i];
     }
   delete [] xx;
   delete [] xxold;
   delete [] nn;
   delete [] nnold;
}
//-----------------------------------------------------------

int readparameters()
{
  ifstream filepar;
  filepar.open("Estera_lim.conf");
  filepar >> sinputfile;
  searchendline(&filepar);
  filepar >> Nofiles;
  searchendline(&filepar);
  filepar >> Nstart;
  searchendline(&filepar);
  filepar >> Ndata;
  searchendline(&filepar);
  filepar >> Nwin;
  searchendline(&filepar);
  filepar >> shift;
  searchendline(&filepar);
  filepar >> Automatic;
  searchendline(&filepar);
  filepar >> fauto;
  searchendline(&filepar);
  filepar >> gausssurr;
  searchendline(&filepar);
  filepar >> Ndatagauss;
  searchendline(&filepar);
  filepar >> delinear;
  searchendline(&filepar);
  filepar >> idiffmax;
  searchendline(&filepar);
  filepar >> Chbmultifitting;
  searchendline(&filepar);
  filepar >> Chbcontrol;
  searchendline(&filepar);
  filepar >> fittingiterations;
  searchendline(&filepar);
  filepar >> beta1;
  searchendline(&filepar);
  filepar >> epsmin;
  searchendline(&filepar);
  filepar >> epsmax2;
  searchendline(&filepar);
  filepar >> epsnumber;
  searchendline(&filepar);
  filepar >> dSmin;
  searchendline(&filepar);
  filepar >> tau;
  searchendline(&filepar);
  filepar >> tau2;
  searchendline(&filepar);
  filepar >> powerepsmin;
  searchendline(&filepar);
  filepar >> powerepsmax;
  searchendline(&filepar);
  filepar >> powernumber;
  searchendline(&filepar);
  filepar >> omegapow;
  searchendline(&filepar);
  filepar >> Chbsig;
  searchendline(&filepar);
  filepar >> ChbK;
  searchendline(&filepar);
  filepar >> Chba;
  searchendline(&filepar);
  filepar >> Chbb;
  searchendline(&filepar);
  filepar >> Chbc;
  searchendline(&filepar);
  filepar >> Chbd;
  searchendline(&filepar);
  filepar >> Chbcn;
  searchendline(&filepar);
  filepar >> conditioncut;
  searchendline(&filepar);
  filepar >> Chbpoprawka;
  searchendline(&filepar);
  filepar >>   shufflepar;
  filepar.close();
  return 0;
}

//------------------------------------------------

int Final()
{
  if(spropen>=3)
    {
      spropen=5;
      Fitting1();
      spropen=3;
    }
  return 0;
}

//------------------------------


//-------------------------------------------------

int main()
{
  FILE *filenames;
  ofstream fileout;
  printf("The Program to evaluate the noise level. Created by K. Urbanowicz (September 2005)\n");

  readparameters();
  alfa=new double*[10];
  assert(alfa);
  for(i=0;i<10;i++)
    {
      alfa[i]=new double[10];
      assert(alfa[i]);
    }
  covar=new double*[10];
  assert(covar);
  for(i=0;i<10;i++)
    {
      covar[i]=new double[10];
      assert(covar[i]);
    }
  long Ndataconst=Ndata,Nstartconst=Nstart,Ndatag2;
  if(Nwin>Ndatagauss)
     Ndatag2=Nwin;
  else
     Ndatag2=Ndatagauss;
   X=new double[Ndatag2];
   assert(X);
   XOLD=new double[Ndatag2];
   assert(XOLD);
  filenames=fopen(sinputfile,"r");
  assert(filenames);
  for(int skl=0;skl<Nofiles;skl++)
  {
   Ndatag2=Ndatagauss;
   Nstart=Nstartconst;
   Ndata=Nwin;
   fscanf(filenames,"%s",sinputfile);
   sprintf(soutputfile,"%s.est.txt",sinputfile);
   fileout.open(soutputfile);
   assert(fileout);
   likwidendl(sinputfile);
   fileout<<"Nstart\t Ndata\t NTS\t errorNTS\t sigma\t errorsigma\t SDofData\tNoise\tLinear\tChaos\n";
   for(;Nstart+Nwin<=Ndataconst;Nstart+=shift)
   {
    Open(sinputfile,&Ndata);
    Ndata=Ndata-Ndata%2;
    if(Ndata<Nwin-10) break;
    printf("\nRead %ld-%ld data of %s file\n",Nstart,Nstart+Ndata,sinputfile);
        mean=0;
    for(i=0;i<Ndata;i++)
      mean+=X[i];
    if(Ndata)
      mean/=Ndata;
    stddev=0;
    for(i=0;i<Ndata;i++)
      stddev+=(X[i]-mean)*(X[i]-mean);
    if(Ndata)
      stddev/=Ndata;
    if(stddev>=0)
      stddev=sqrt(stddev);
    double stddevold2=stddev;
    if(gausssurr!=0)
    {
      printf("Generate Gaussian substitution\n");
      Gaussiansurrogates(&Ndatag2);
    }
    if(shufflepar==1)
    {
       printf("Make a surrogate data by shuffling\n");
       shuffle(3);
    }
  double stddev_delinear = stddev;
  stddev=stddevold2;

  printf("Standard deviation of data = %f\nCalculation of the entropy\n",stddev);
  if(stddev <= 0.00001)
  {
    continue;
  }
  Calculateentropy();
  printf("Fitting procedure\n");
  int  kcount=0;
  do
    {
    if(kcount<2)
         spropen=2;
    else
    {
        spropen=3;
        ia[1]=1;
        ia[2]=1;
        ia[3]=1;
        ia[4]=1;
        ia[5]=1;
        ia[6]=1;
    }
      Fitting1();
      for(int pl=0;pl<2;pl++)
	        Final();
      if(stddev>1e-5)
        cout << "\n\t"<<(kcount+1)<< "\tapproximation\nNTS ="<< (a[1]*dziel/stddev)<<"\terrorNTS = "<<(errorSIG*dziel/stddev)<<endl;
      cout << "SD of noise = "<< (a[1]*dziel)<<"\terror = "<<(errorSIG*dziel)<<endl;

    }while(errorSIG*dziel/stddev>=0.06&&++kcount<3);
    if(stddev>1e-5)
    {
        double nts = (a[1]*dziel/stddev);
        double nts2 = nts*nts;
        double linear = 1 - stddev_delinear * stddev_delinear / stddev / stddev;
        if(linear < 0)
            linear = 0;
        
        if(nts2 + linear > 1)
        {
            linear = 1 - nts2;
        }
        double chaos = 1 - nts2 - linear;

     fileout<<Nstart<<"\t"<<Ndata<<"\t"<<nts<<"\t"<<(errorSIG*dziel/stddev)<<"\t"<<(a[1]*dziel)<<"\t"<<(errorSIG*dziel)<<"\t"<<stddev<<"\t"<<nts2<<"\t"<<linear<<"\t"<<chaos<<endl;
    }
  if(stddev>1e-5)
    cout << "\n\tFinal result for "<<Nstart <<" - "<<(Nstart+Ndata) <<" data of "<<sinputfile <<" file\nNTS = "<<(a[1]*dziel/stddev) <<"\terrorNTS = "<<(errorSIG*dziel/stddev)<<endl;
   cout << "SD of noise = "<<(a[1]*dziel) <<"\terror = "<<(errorSIG*dziel)<<endl;

   }
   fileout.close();

  }
  fclose(filenames);
  delete [] X;
  delete [] XOLD;
  for(i=0;i<10;i++)
  {
      delete [] alfa[i];
      delete [] covar[i];
  }
  delete [] alfa;
  delete [] covar;
  return 0;
}

