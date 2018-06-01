#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

#define NBIN 1000

/**LOCAL VARIABLES**/
static int ***bincnt,**cnt;
static double **sig;
static float binwidth=10;

/**LOCAL FUNCTIONS**/
void allocate_arrays(int n);

void allocate_arrays(int n)
{
  int i,j,k;

  fprintf(stderr,"Allocating histogram arrays.\n");

  bincnt=i3tensor(1,n,1,n,-NBIN,NBIN);
  sig=dmatrix(1,n,1,n);
  cnt=imatrix(1,n,1,n);

  for(i=1;i<=n;++i)
    for(j=1;j<=n;++j)
      {
	sig[i][j]=0;
	cnt[i][j]=0;
	for(k=-NBIN;k<=NBIN;++k)
	  bincnt[i][j][k]=0;
      }
}

void output_histograms(char *a,int nbin, double radius[])
{
  float fac,center;
  int i,j,k,n,binfac,sum,outcnt;
  FILE *fp;
  char fname[100],rt[100];

  /* Get the root of the filename from stdin
   */
  fscanf(stdin,"%s",rt);
  fprintf(stderr,"Root of the filename [%s]\n",rt);
  

  j=1;
  for(i=1;i<=nbin;++i)
    {
      sig[i][j]=sqrt(sig[i][j]/cnt[i][j]);
      if(isnan(sig[i][j]))sig[i][j]=0;
      
      sprintf(fname,"htan%02d.%s",i-1,rt);
      fp=fopen(fname,"w");
      
      fac=(log10(cnt[i][j])*2-2);
      if(fac<1.5)fac=1.5;
      fac=sig[i][j]/fac;
      
      binfac=fac/binwidth;
      /*
      fprintf(stderr,"Binwidth for [%s] is %.0f %.1f %d\n",
	      fname,binfac*binwidth,sig[i][j],cnt[i][j]);
      */
      outcnt=sum=n=0;
      center=0;
      for(k=-NBIN;k<=NBIN;++k)
	{
	  /*fprintf(fp,"%f %d\n",binwidth*(k+0.5),bincnt[i][j][k]);*/
	  n++;
	  sum+=bincnt[i][j][k];
	  center+=binwidth*(k+0.5);
	  if(n==binfac)
	    {
	      outcnt++;
	      fprintf(fp,"%f %d\n",center/n,sum);
	      sum=n=center=0;
	    }
	}
      fclose(fp);
    }


  i=1;
  for(j=1;j<=nbin;++j)
    {
      sig[i][j]=sqrt(sig[i][j]/cnt[i][j]);
      if(isnan(sig[i][j]))sig[i][j]=0;
      
      sprintf(fname,"hrad%02d.%s",j-1,rt);
      fp=fopen(fname,"w");
      
      fac=(log10(cnt[i][j])*2-2);
      if(fac<1.5)fac=1.5;
      fac=sig[i][j]/fac;
      
      binfac=fac/binwidth;
      if(binfac<3)binfac=3;
      /*
      fprintf(stderr,"Binwidth for [%s] is %.0f %.1f %d\n",
	      fname,binfac*binwidth,sig[i][j],cnt[i][j]);
      */
      outcnt=sum=n=0;
      center=0;
      for(k=-NBIN;k<=NBIN;++k)
	{
	  /*fprintf(fp,"%f %d\n",binwidth*(k+0.5),bincnt[i][j][k]);*/
	  n++;
	  sum+=bincnt[i][j][k];
	  center+=binwidth*(k+0.5);
	  if(n==binfac)
	    {
	      outcnt++;
	      fprintf(fp,"%f %d\n",center/n,sum);
	      sum=n=center=0;
	    }
	}
      fclose(fp);
    }

}

void tabulate_histograms(float v, int kbin, int rbin, int nbin)
{
  int i,j;
  float fac;
  FILE *fp;
  static int flag=0;
  
  if(kbin>=nbin || rbin>=nbin)return;
  if(isnan(v) || isinf(v))return;
  if(!flag++)allocate_arrays(nbin);

  /* We want our bins to go from 1..nbin rather than zero.
   */
  kbin++;
  rbin++;

  cnt[kbin][rbin]++;
  sig[kbin][rbin]+=v*v;

  fac=v/binwidth;
  j=(int)fac;
  if(fac<0)j--;
  if(abs(j)<=NBIN)
    bincnt[kbin][rbin][j]++;
}

