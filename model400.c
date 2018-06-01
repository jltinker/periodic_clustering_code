#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

#define NMAX 50
#define N2d  100
#define pi 3.1415926535898

/** NUMERICAL RECIPES **/
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);


/** INTERNAL ROUTINES **/
float gauss1(float x, float s);
float interpolate2d(float x, float y, float **xi, int n, float rmax);

int main(int argc, char **argv)
{
  float *v,*r,*xi,*sig,**xi2d,rsig,rpi,y,realrad,*y2xi,*y2sig,*y2v,x1,x2,x3,x4,x5,x6,
    v1,s1,dv,dy,rmax=0,ymin,ymax;
  int i=0,j,k,npts,ncol,Ny;
  FILE *fp;

  r=vector(1,NMAX);
  v=vector(1,NMAX);
  sig=vector(1,NMAX);
  xi=vector(1,NMAX);
  xi2d=matrix(1,N2d,1,N2d);
  y2sig=vector(1,NMAX);
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
  rmax=rmax/sqrt(2);

  /* Integrate in linear space, then convert to log space (later).
   */
  for(i=1;i<=N2d;++i)
    for(j=1;j<=N2d;++j)
      {
	rsig=(i-0.5)/N2d*rmax;
	rpi=(j-0.5)/N2d*rmax;
	ymin=r[1];
	ymax=sqrt(rmax*rmax-rsig*rsig);
 	Ny=ymax/ymin*0.5;
	dy=ymax/Ny;

	for(xi2d[i][j]=0,k=-Ny+1;k<=Ny;++k)
	  {
	    y=(k-0.5)/Ny*ymax;
	    realrad=sqrt(y*y+rsig*rsig);
	    
	    v1=0;
	    s1=400*sqrt(1);

	    dv=(rpi-y)*100.0-v1*y/realrad;
	    xi2d[i][j]+=x1*gauss1(dv,s1)*dy*100/2;
	  }
      }

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

  for(j=1;j<i;++j)
    {
      fscanf(fp,"%f %f %f",&x1,&rsig,&rpi);
      x2=interpolate2d(rsig,rpi,xi2d,N2d,rmax);
      printf("%f %f %f %f\n",x1,rsig,rpi,x2);
    }

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

