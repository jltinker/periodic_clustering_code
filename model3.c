#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

/** Same as model.c, but this time we're using two-sided exponentials rather than
 ** 1-sigma Gaussians. These have to be rad in from a seperate file.
 */

#define rt2  1.414214
#define NMAX 50
#define N2d  25
#define pi 3.1415926535898

/** NUMERICAL RECIPES **/
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);


/** INTERNAL ROUTINES **/
float exp1(float x, float s);
float gauss1(float x, float s);
float interpolate2d(float x, float y, float **xi, int n, float rmax);

int main(int argc, char **argv)
{
  float *v,*r,*xi,*sig,**xi2d,rsig,rpi,y,realrad,*y2xi,*y2sig,*y2v,x1,x2,x3,x4,x5,x6,
    v1,s1,dv,dy,rmax=0,ymin,ymax,*sign,*sigp,*y2sign,*y2sigp;
  int i=0,j,k,npts,ncol,Ny,i1,i2;
  FILE *fp;

  r=vector(1,NMAX);
  v=vector(1,NMAX);
  sig=vector(1,NMAX);
  sign=vector(1,NMAX);
  sigp=vector(1,NMAX);
  xi=vector(1,NMAX);
  xi2d=matrix(1,N2d,1,N2d);
  y2sig=vector(1,NMAX);
  y2sign=vector(1,NMAX);
  y2sigp=vector(1,NMAX);
  y2xi=vector(1,NMAX);
  y2v=vector(1,NMAX);

  /* Input the correlation data (+other data).
   */
  if(!(fp=fopen(argv[1],"r")))
    kill("Error opening input file.\n");
  while(!feof(fp))
    {
      i++;
      fscanf(fp,"%f %f %f %f %f %f %f %f",&r[i],&x1,&v[i],&x2,&sig[i],&x3,&xi[i],&x4);
      v[i]*=1;
      if(r[i]>rmax)rmax=r[i];
    }
  npts=i-1;
  fclose(fp);
  rmax/=sqrt(2);


  /* Input exponential fit parameters
   */
  if(!(fp=fopen(argv[3],"r")))
    kill("Error opening stdev file.\n");
  for(i=1;i<=npts;++i)
    fscanf(fp,"%f %f %d %d %f %f %f",&sign[i],&sigp[i],&i1,&i2,&x1,&x2,&x3);

  /* Set up the spline interpolation routine for xi(r), v(r), sig(r).
   */
  spline(r,xi,npts,2.0E+30,2.0E+30,y2xi);
  spline(r,sig,npts,2.0E+30,2.0E+30,y2sig);
  spline(r,sigp,npts,2.0E+30,2.0E+30,y2sigp);
  spline(r,sign,npts,2.0E+30,2.0E+30,y2sign);
  
  /* Integrate in linear space, then convert to log space (later).
   */
  for(i=1;i<=N2d;++i)
    for(j=1;j<=N2d;++j)
      {
	rsig=(i-0.5)/N2d*rmax;
	rpi=(j-0.5)/N2d*rmax;
	ymin=r[1];
	ymax=sqrt(2*rmax*rmax-rsig*rsig);
 	Ny=ymax/ymin*4.5;
	dy=ymax/Ny;

	for(xi2d[i][j]=0,k=-Ny+1;k<=Ny;++k)
	  {
	    y=(k-0.5)/Ny*ymax;
	    realrad=sqrt(y*y+rsig*rsig);	
	    dv=(rpi-y)*100.0;
	    if(y<0)dv=-dv;

	    splint(r,xi,y2xi,npts,realrad,&x1);
	    if(dv>0)
	      splint(r,sigp,y2sigp,npts,realrad,&s1);
	    else
	      splint(r,sign,y2sign,npts,realrad,&s1);

	    xi2d[i][j]+=(x1+1)*exp1(dv,s1)*dy*100;

	    /*
	    splint(r,xi,y2xi,npts,sqrt(rsig*rsig+rpi*rpi),&x2);
	    if(j==9 && i==9)
	      printf("VEL1 %f %f %f %f %f %e %f %f %f %f\n",
		     realrad,rsig,rpi,y,dv,exp1(dv,s1),dy,x1,xi2d[i][j],x2);
	    
	    if(j==15 && i==3)
	      printf("VEL2 %f %f %f %f %f %e %f %f %f %f\n",
		     realrad,rpi,rsig,y,dv,exp1(dv,s1),dy,x1,xi2d[i][j],x2);
	    */
	  }
	/*
	if(i==4)
	  printf("BUH %f %f %e\n",rsig,rpi,xi2d[i][j]);
	*/
      }

  for(i=1;i<=N2d;++i)
    for(j=1;j<=N2d;++j)
      xi2d[i][j]-=1;


  /* Read in the xi(sigma,pi) from the file and 
   * interpolate the model onto the same positions.
   */
  if(!(fp=fopen(argv[2],"r")))
    kill("Error opening 2dcorr file.");

  i=0;
  while(!feof(fp))
    {
      fscanf(fp,"%f %f %f",&x1,&rsig,&rpi);
      i++;
    }
  rewind(fp);

  /* Test the interpolation scheme.
   *
  rsig=(10-0.5)/N2d*rmax;
  i=1;
  printf("%f %f %f\n",(i-0.5)/N2d*rmax,xi2d[10][i],0.0);
  for(i=1+1;i<=N2d-1;++i)
    printf("%f %f %f\n",(i-0.5)/N2d*rmax,xi2d[10][i],
	   interpolate2d(rsig,(i-0.5)/N2d*rmax,xi2d,N2d,rmax));
  exit(0);
  */

  for(j=1;j<i;++j)
    {
      fscanf(fp,"%f %f %f",&x1,&rsig,&rpi);
      x2=interpolate2d(rsig,rpi,xi2d,N2d,rmax);
      splint(r,xi,y2xi,npts,sqrt(rpi*rpi+rsig*rsig),&x3);
      printf("%f %f %f %f %f\n",x1,rsig,rpi,x2,x3);
    }
  exit(0);

  for(i=1;i<=N2d;++i)
    for(j=1;j<=N2d;++j)
      fprintf(stdout,"%e %e %f %f\n",xi2d[i][j],xi2d[j][i],(i-0.5)/N2d*rmax,(j-0.5)/N2d*rmax);

  /* Get things in a format easy to read in 
   * sm columns
   */
  ncol=sqrt(i);
  for(j=1;j<=ncol;++i)
    {
      for(k=1;k<=ncol;++k)
	{
	  fscanf(fp,"%f %f %f",&x1,&rsig,&rpi);
	  x2=interpolate2d(rsig,rpi,xi2d,N2d,rmax);
	  printf("");
	}
    }

}


/* Little ditty to calculate the exponential
 * with x, sig (centered on zero).
 */
float exp1(float x, float s)
{
  return(1.0/(rt2*s)*exp(-rt2*fabs(x)/s));
}
    
/* Little ditty to calculate the Gaussian
 * with x, sig (centered on zero).
 */
float gauss1(float x, float s)
{
  return(1.0/(sqrt(2*pi)*s)*exp(-x*x*0.5/(s*s)));
}
    
/* Returns the interpolated value of xi(sig,pi)
 * from the model (which is on a cartesian grid
 * of spacing rmax/N2d.
 */
float interpolate2d(float x, float y, float **xi, int n, float rmax)
{
#define cnint(xx) ((xx-floor(xx)) < 0.5 ? floor(xx) : ceil(xx))

  int ip1,jp1,ip,jp;
  float x1,t,u,t1,u1,a,b;

  a=x;
  b=y;

  if(x>rmax || y>rmax)
    kill("rmax too small for input xi2d.");

  /* Convert spatial coords into grid coords
   */
  x=x/rmax*n+1;
  y=y/rmax*n+1;
  
  if(x<1 || y<1)
    kill("N2d too small.");
  
  ip1=cnint(x);
  jp1=cnint(y);
  
  ip=ip1-1;
  jp=jp1-1;
  
  t=x-(ip+0.5);
  u=y-(jp+0.5);
  
  t1=1-t;
  u1=1-u;

  x1=t1*u1*xi[ip][jp] + 
    t*u1*xi[ip1][jp] +
    t1*u*xi[ip][jp1] + 
    t*u*xi[ip1][jp1];

  return(x1);
}

