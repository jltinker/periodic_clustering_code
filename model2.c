#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

/* This is basically the same as the model.c program, only this time we actually
 * have hard-wired distribution functions f(v) along the line of observer's sight.
 * (Histograms of N-body data.)
 *
 * New and improved! Now we will be using 2-d histograms--> this should be EXACTLY 
 * the pdf(v) that is involved.
 */

#define NMAX 50
#define N2d  50
#define pi 3.1415926535898
#define NHMAX 1200

/** NUMERICAL RECIPES **/
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);


/** INTERNAL ROUTINES **/
float gauss1(float x, float s);
float interpolate2d(float x, float y, float **xi, int n, float rmax);
void get_histograms(float ***v, float ***n, int **nbin, int npts, char *root, char *fn2);
float prob2d(float v, float x, float y, float rmax);
float internal_bin_sum(float **xi2d, int n2d, float rmax, float xhi, float yhi);

float ***nhist, ***vhist, *rbin, ***y2hist,**rxbin,**rybin,VSCALE=1;
int **nbins,NPTS;

int main(int argc, char **argv)
{
  float *v,*r,*xi,*sig,**xi2d,rsig,rpi,y,realrad,*y2xi,*y2sig,*y2v,x1,x2,x3,x4,x5,x6,
    v1,s1,dv,dy,rmax=0,ymin,ymax,err,sighi,pihi,xx,yy;
  int i=0,j,k,npts,ncol,Ny,SIGN=-1;
  FILE *fp;

  if(argc<6)
    kill("model2 xifile xi2dfile nhist path rfile [SIGN]\n");

  if(argc>6)
    SIGN=atof(argv[6]);

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
  rmax/=sqrt(2);
  if(rmax<20)rmax=20;
  fprintf(stderr,"RMAX: %f\n",rmax);

  /* Get the f(v) histograms
   */
  NPTS=atoi(argv[3]);
  vhist=f3tensor(1,NPTS,1,NPTS,1,NHMAX);
  nhist=f3tensor(1,NPTS,1,NPTS,1,NHMAX);
  nbins=imatrix(1,NPTS,1,NPTS);
  rbin=vector(1,NPTS);
  y2hist=f3tensor(1,NPTS,1,NPTS,1,NHMAX);
  rxbin=matrix(1,NPTS,1,NPTS);
  rybin=matrix(1,NPTS,1,NPTS);

  get_histograms(vhist,nhist,nbins,NPTS,argv[4],argv[5]);

  /* Set up the spline interpolation routine for xi(r)
   */
  spline(r,xi,npts,2.0E+30,2.0E+30,y2xi);
  
  /* Integrate in linear space, then convert to log space (later).
   */
  for(i=1;i<=N2d;++i)
    for(j=1;j<=N2d;++j)
      { 
	rsig=(i-0.5)/N2d*rmax;
	rpi=(j-0.5)/N2d*rmax;	
	ymin=r[1];
	ymax=sqrt(2*rmax*rmax-rsig*rsig);
 	Ny=ymax/ymin*2.5;
	dy=ymax/(Ny);

	for(xi2d[i][j]=0,k=-Ny+1;k<=Ny;++k)
	  {
	    y=(k-0.5)/Ny*ymax;
	    realrad=sqrt(y*y+rsig*rsig);
	    
	    splint(r,xi,y2xi,npts,realrad,&x1);
	    dv=(rpi-y)*100.0;
	    x3=dv;
	    if(y<0)dv=SIGN*dv;

	    xi2d[i][j]+=(1+x1)*prob2d(dv,rsig,y,rmax)*dy*100;
	    
	    splint(r,xi,y2xi,npts,sqrt(rsig*rsig+rpi*rpi),&x2);
	    if(j==20 && i==32)
	      printf("VEL1 %f %f %f %f %f %e %f %f %f %f\n",
		     realrad,rsig,rpi,y,dv,prob2d(dv,rsig,y,rmax),dy,x1,xi2d[i][j],x2);
	    
	    if(j==2 && i==32)
	      printf("VEL2 %f %f %f %f %f %e %f %f %f %f\n",
		     realrad,rpi,rsig,y,x3,prob2d(dv,rsig,y,rmax),dy,x1,xi2d[i][j],x2);
	    
	  }
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
      fscanf(fp,"%f %f %f %f %f %f",&x1,&rsig,&rpi,&x1,&x1,&x1);
      i++;
    }
  rewind(fp);

  xx=0;
  yy=0;
  for(j=1;j<i;++j)
    {
      fscanf(fp,"%f %f %f %f %f %f",&x1,&rsig,&rpi,&err,&sighi,&pihi);
      x2=internal_bin_sum(xi2d,N2d,rmax,sighi,pihi);
      if(x2==0)x2=interpolate2d(rsig,rpi,xi2d,N2d,rmax);

      /*
      x2=interpolate2d(xx,yy,xi2d,N2d,rmax);
      if(j%10==0)xx=sighi;
      yy=pihi;
      */

      splint(r,xi,y2xi,npts,sqrt(rpi*rpi+rsig*rsig),&x3);
      printf("%f %f %f %f %f %f\n",x1,rsig,rpi,x2,x3,err/sqrt(5.));
    }
  exit(0);
 
  for(i=1;i<=N2d;++i)
    for(j=1;j<=N2d;++j)
      fprintf(stdout,"%e %e %f %f\n",xi2d[i][j],xi2d[j][i],(i-0.5)/N2d*rmax,(j-0.5)/N2d*rmax);
  exit(0);

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
    {
      fprintf(stderr,"rmax too small for input xi2d.\n");
      if(x>rmax)x=rmax;
      if(y>rmax)y=rmax;
    }

  /* Convert spatial coords into grid coords
   */
  x=x/rmax*n+1;
  y=y/rmax*n+1;

  if(x<1 || y<1)
    kill("N2d too small.");
  
  ip1=cnint(x);
  jp1=cnint(y);
  
  if(ip1>=n+1)
    ip1--;
  if(jp1>=n+1)
    jp1--;

  ip=ip1-1;
  jp=jp1-1;
  
  if(ip<1)
    {
      ip++;
      ip1++;
    }
  if(jp<1)
    {
      jp++;
      jp1++;
    }

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

/* Function to use the histograms to return the probability of having the
 * input velocity at input coordinates.
 */
float prob2d(float v, float x, float y,float rmax)
{
  int i,j,k,ix,iy,imin,jmin,n;
  float p,dx,dy,r,rmin,sum=0,wsum=0;
  static int flag=0;

  n=NPTS;
  y=fabs(y);

  /* Initialize the spline interpolation
   */
  if(!flag++)
    for(i=1;i<=NPTS;++i)
      for(j=1;j<=NPTS;++j)
	spline(vhist[i][j],nhist[i][j],nbins[i][j],2.0E+30,2.0E+30,y2hist[i][j]);

  /* For now, we'll just use the NGP method.
   */
  i=1;
  while(rbin[i]<x && i<NPTS)i++;
  j=1;
  while(rbin[j]<y && j<NPTS)j++;
  ix=i;
  iy=j;

  
  rmin=100;
  for(i=ix-1;i<=ix+1;++i)
    for(j=iy-1;j<=iy+1;++j)
      {
	if(i<1 || j<1)continue;
	if(i>NPTS || j>NPTS)continue;
	dx=rxbin[i][j]-x;
	dy=rybin[i][j]-y;
	r=sqrt(dx*dx+dy*dy);
	splint(vhist[i][j],nhist[i][j],y2hist[i][j],nbins[i][j],v,&p);

	sum+=1/r/r;
	wsum+=p/r/r;
	if(r<rmin)
	  {
	    rmin=r;
	    imin=i;
	    jmin=j;
	  }
      }
  i=imin;
  j=jmin;

  /*
  if((rbin[i]-x > x-rbin[i+1]) && i!=1 && i!=NPTS)
    i++;
  if((rbin[j]-x > x-rbin[j+1]) && j!=1 && j!=NPTS)
    j++;
  if(i!=NPTS)
    i++;
  if(j!=NPTS)
    j++;
  if(i!=1)
    i--;
  if(j!=1)
    j--;
  */
  
  splint(vhist[i][j],nhist[i][j],y2hist[i][j],nbins[i][j],v,&p);
  /*printf("H2 %.2f %.2f %.1f %d %d %.2f %.2f %f %f\n",x,y,v,i,j,rxbin[i][j],rybin[i][j],p,rmin);*/
  fflush(stdout);
  if(p<0)return(0);
  return(p);
}




/* File to retrieve the f(v) histograms from
 * some pre-determined location.
 */

void get_histograms(float ***v, float ***n, int **nbin, int npts, char *root, char *fn2)
{
  int i,j,k;
  float binwidth,x,ntot;
  FILE *fp;
  char fn[100];

  for(i=1;i<=npts;++i)
    for(j=1;j<=npts;++j)
      {
	sprintf(fn,"%s/h2d.%02d.%02d",root,i,j);
	if(!(fp=fopen(fn,"r")))
	  {
	    fprintf(stderr,"Error opening [%s]\n",fn);
	    kill(" ");
	  }
	fflush(stderr);
	k=ntot=0;
	while(!feof(fp))
	  {
	    k++;
	    fscanf(fp,"%f %f",&v[i][j][k],&n[i][j][k]);
	    ntot+=n[i][j][k];
	    v[i][j][k]*=VSCALE;
	  }
	nbin[i][j]=k-1;
	binwidth=v[i][j][2]-v[i][j][1];
	for(k=1;k<=nbin[i][j];++k)
	  n[i][j][k]/=(ntot*binwidth);
	fclose(fp);
      }
 
  /* Now get the radii for each radial bin.
   */
  if(!(fp=fopen(fn2,"r")))
    {
      fprintf(stderr,"Error opening [%s]\n",fn);
      kill(" ");
    }
  for(i=1;i<=npts;++i)
    fscanf(fp,"%f %f",&x,&rbin[i]);
  rewind(fp);
    
  for(i=1;i<=npts;++i) 
    for(j=1;j<=npts;++j) 
      fscanf(fp,"%f %f",&rxbin[i][j],&rybin[i][j]);
}

float internal_bin_sum(float **xi2d, int n2d, float rmax, float xhi, float yhi)
{
  float r,dx,dy,sum=0,wsum=0,x,y,x2;
  int i,j,nx,ny;
  static float xlo=0,ylo=0,xprev=0;
  
  /* If the bin is smaller than 1 h^-1 Mpc in both directions,
   * then just interpolate and return 0 for interpolation
   */
  if(xhi-xlo<=1 && yhi-ylo<=1)
    {
      if(ylo>yhi)
	{
	  xlo=xprev;
	  ylo=0;
	}
      xprev=xhi;
      ylo=yhi;
      return(0);
    }
  if(ylo>yhi)
    {
      xlo=xprev;
      ylo=0;
    }

  nx=(xhi-xlo)+4;
  ny=(yhi-ylo)+4;
  dx=(xhi-xlo)/nx;
  dy=(yhi-ylo)/ny;

  for(i=1;i<=nx;++i)
    for(j=1;j<=ny;++j)
      {
	x=xlo+(i-0.5)*dx;
	y=ylo+(j-0.5)*dy;
	x2=interpolate2d(x,y,xi2d,n2d,rmax);
	wsum+=x*(1+x2);
	sum+=x;
      }
  wsum=wsum/sum-1;
  xprev=xhi;
  ylo=yhi;
  return(wsum);

}
