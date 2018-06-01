#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "nrutil.h"

void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void calc_gamma(char *fname, float g[], float r[]);
void least_squares(float *x, float *y, int n, float *a, float *b);

int NBODY_FORMAT=1,NSIDE=15;
float RATIO=0.5;

int main(int argc, char **argv)
{
  int i,j,k,NRUNS=5,flag=0;
  char fname[100];
  float *gs,*g2,*r,*x,*g,*gt;
  
  
  if(argc<2)
    {
      fprintf(stderr,"> create_linear_plot xi2d.file [NBODY FORMAT (1/0)] [ratio] [NRUNS]> outfile\n");
      exit(0);
    }

  if(argc>2)
    NBODY_FORMAT=atoi(argv[2]);
  if(argc>3)
    RATIO=atof(argv[3]);
  if(argc>4)
    NRUNS=atof(argv[4]);

  NSIDE=0;
  if(NRUNS>1)
    sprintf(fname,"run%d/%s.%d",1,argv[1],1);
  else
    sprintf(fname,"%s",argv[1],1);
  calc_gamma(fname,g,x);


  gs=vector(1,NSIDE);
  g2=vector(1,NSIDE);
  gt=vector(1,NSIDE);
  g=vector(1,NSIDE);
  r=vector(1,NSIDE);
  x=vector(1,NSIDE);
  for(i=1;i<=NSIDE;++i)
    gs[i]=g2[i]=r[i]=0;

  for(i=1;i<=NRUNS;++i)
    {
      for(k=1;k<=NSIDE;++k)
	gt[k]=0;
      for(j=1;j<=3;++j)
	{
	  if(NRUNS>1)
	    sprintf(fname,"run%d/%s.%d",i,argv[1],j);
	  else
	    sprintf(fname,"%s",argv[i],j);
	  calc_gamma(fname,g,x);
	  if(NRUNS==1)
	    {
	      for(k=1;k<=NSIDE;++k)
		printf("%e %e\n",x[k],g[k]);
	      exit(0);
	    }
	  for(k=1;k<=NSIDE;++k)
	    gt[k]+=g[k]/3;
	}
      for(k=1;k<=NSIDE;++k)
	{
	  gs[k]+=gt[k];
	  r[k]+=x[k];
	  g2[k]+=gt[k]*gt[k];
	}
    }

  for(i=1;i<=NSIDE;++i)
    r[i]/=(1.0*NRUNS);
  for(i=1;i<=NSIDE;++i)
    gs[i]/=(1.0*NRUNS);
  for(i=1;i<=NSIDE;++i)
    g2[i]/=(1.0*NRUNS);
  for(i=1;i<=NSIDE;++i)
    g2[i]=sqrt((g2[i]-gs[i]*gs[i])/(1.0*NRUNS-1));

  for(i=1;i<=NSIDE;++i)
    fprintf(stdout,"%e %e %e\n",r[i],gs[i],g2[i]);
}

void calc_gamma(char *fname, float g[], float rr[])
{
  int i,j,nsize,ix,iy,nbins,sm_out=0,niter;
  float x1,x2,x3,x4,x5,err,xhi,xlo,rsigma,b,RATIO2=0.1,xa[3],ya[3];
  float r,rs,rp,dx1,dx2,dx3,dx4,*xi,*xx,*yy,*y2,rmax,a,error=1,tolerance=0.01;
  FILE *fp,*fp2;
  char aa[1000];
  
  j=-1;
  if(!(fp=fopen(fname,"r")))
    {
      fprintf(stderr,"ERROR opening [%s]\n",fname);
      exit(0);
    }
  if(NBODY_FORMAT)
    {
      while(!(feof(fp)))
	{
	  fscanf(fp,"%f %f %f %f",&x1,&x5,&x4,&err);
	  fgets(aa,1000,fp);
	  j++;
	}
    }
  else
    {
      while(!(feof(fp)))
	{
	  fscanf(fp,"%f %f %f %f %f %f %f",&x1,&x5,&x4,&err,&x2,&x3,&x3);
	  j++;
	}
    }    
  rewind(fp);
  
  nsize=sqrt(j*1.0);
  if(!NSIDE)
    {
      NSIDE=nsize;
      return;
    }
  fprintf(stderr,"Read %d lines from file [%s] %d.\n",j,fname,nsize);
  xi=vector(1,nsize);
  xx=vector(1,nsize);
  yy=vector(1,nsize);
  y2=vector(1,nsize);

  for(j=1;j<=nsize;++j)
    {
      for(i=1;i<=nsize;++i)
	{
	  if(NBODY_FORMAT)
	    fscanf(fp,"%f %f %f %f",&x1,&x5,&x4,&err);
	  else
	    fscanf(fp,"%f %f %f %f %f %f %f",&x5,&x4,&x1,&x2,&x3,&x3,&err);
	  fgets(aa,1000,fp);

	  xi[i]=x4;
	  xx[i]=x5;
	  yy[i]=x1;
	}
      x1=xi[1];

      /* Fit a straight line (in log space) to the first three points
       * and interpolate the value at r_pi=0.1 Mpc/h.
       * 
       * NEW-- exclude the innermost point.
       */
      for(i=0;i<3;++i)
	{
	  //xa[i]=log(xx[i+1]);
	  //ya[i]=log(xi[i+1]);
	  xa[i]=log(xx[i+2]);
	  ya[i]=log(xi[i+2]);
	}
      least_squares(xa,ya,3,&a,&b);
      x1=exp(a+b*log(0.1));
      
      // take MEAN
      x1 = (xi[1] + xi[2] + xi[3])/3.0;

      for(i=1;i<=nsize;++i)
	{
	  xi[i]/=x1;
	}
      spline(xx,xi,nsize,1.0E+30,1.0E+30,y2);
      i=1;
      while(xi[i]>RATIO && i<=nsize)i++;

      /* Find the location, r_pi, where the correlation function 
       * drops by a factor (RATIO)
       */
      xhi=xx[i];
      xlo=xx[i-1];
      error=1;
      niter=0;
      while(error>tolerance && niter<200)
	{
	  niter++;
	  r=0.5*(xlo+xhi);
	  splint(xx,xi,y2,nsize,r,&a);
	  if(a<RATIO)xhi=r;
	  else xlo=r;
	  error=fabs(a-RATIO)/RATIO;
	}
      for(i=1,rsigma=0;i<=nsize;++i)
	rsigma+=yy[i]/nsize;
      
      rr[j]=rsigma;
      g[j]=r;

      continue;

      /* Find the ratio of the correlation functions at the
       * radii (0.1,10.0).
       */
      splint(xx,xi,y2,nsize,0.1,&a);
      splint(xx,xi,y2,nsize,20.0,&b);
      x1=b/a;


      /*
      fprintf(stdout,"%e %e %e ",rsigma,r,x1);
      fflush(stdout);
      */

      /* Find the location, r_pi, where the correlation function 
       * drops by a factor of 100
       */
      i=1;
      while(xi[i]>RATIO2 && i<=nsize)i++;

      xhi=xx[i];
      xlo=xx[i-1];
      error=1;
      niter=0;
      while(error>tolerance && niter<100)
	{
	  niter++;
	  r=0.5*(xlo+xhi);
	  splint(xx,xi,y2,nsize,r,&a);
	  if(a<RATIO2)xhi=r;
	  else xlo=r;
	  error=fabs(a-RATIO2)/(RATIO2);
	}
      /*
      fprintf(stdout,"%e\n",r);
      fflush(stdout);
      */
    }
}
