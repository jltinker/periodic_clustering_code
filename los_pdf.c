#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

/** Same as model.c, but this time we're using the asymmetric probability distribution
 *  assuming that the underlying distributions are exponential and that the mode is at
 *  zero velocity (rather than at the mean).
 */

#define rt2  1.414214
#define NMAX 50
#define N2d  25
#define pi 3.1415926535898
#define RT2PI 2.506628


/** NUMERICAL RECIPES **/
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);


/** INTERNAL ROUTINES **/
float exp1(float x, float s);
float gauss1(float x, float s);
float interpolate2d(float x, float y, float **xi, int n, float rmax);
float vel_prob(float x, float sn, float sp, float st, float rs, float rp);
float term2(double v, float st, float sr, float mu, float sintheta, float costheta);
float vel_prob2g(float x, float sigr1, float sigr2, float sigt1, float sigt2, float mu1, float mu2,
		 float ft, float fr, float rs, float rp);

int main(int argc, char **argv)
{
  float *v,*r,*xi,*sig,**xi2d,rsig,rpi,y,realrad,*y2xi,*y2sig,*y2v,x1,x2,x3,x4,x5,x6,
    v1,s1,dv,dy,rmax=0,ymin,ymax,*sign,*sigp,*sigtan,*y2sign,*y2sigp,*y2sigtan,s1tan,s1p,s1n,
    *sigr1,*sigr2,*sigt1,*sigt2,*mu1,*mu2,*y2sigr1,*y2sigr2,*y2sigt1,*y2sigt2,*y2mu1,*y2mu2,
    sr1,sr2,st1,st2,m1,m2,*fact,*facr,xlo,xhi,ylo,yhi,m1a,m2a,sr1a,sr2a,st1a,st2a,
    m1b,m2b,sr1b,sr2b,st1b,st2b;
  int i=0,j,k,npts,ncol,Ny,i1,i2,nhist;
  FILE *fp,*fp1,*fp2;
  char fname[100];

  /* These are read in from the xi-file (avg of covar3 results)
   */
  r=vector(1,NMAX);
  v=vector(1,NMAX);
  sig=vector(1,NMAX);
  xi=vector(1,NMAX);

  /* These are the parameters of the two-gaussian fit,
   * read in from two files.
   */
  mu1=vector(1,NMAX);
  mu2=vector(1,NMAX);
  sigr1=vector(1,NMAX);
  sigr2=vector(1,NMAX);
  sigt1=vector(1,NMAX);
  sigt2=vector(1,NMAX);
  fact=vector(1,NMAX);
  facr=vector(1,NMAX);

  /* The calculated 2d-correlation function
   */
  xi2d=matrix(1,N2d,1,N2d);

  /* Stuff for the splint interpolation routines.
   */
  y2sig=vector(1,NMAX);
  y2sigr1=vector(1,NMAX);
  y2sigr2=vector(1,NMAX);
  y2sigt1=vector(1,NMAX);
  y2sigt2=vector(1,NMAX);
  y2mu1=vector(1,NMAX);
  y2mu2=vector(1,NMAX);
  y2xi=vector(1,NMAX);
  y2v=vector(1,NMAX);

  if(argc<5)
    kill("model4 xifile g2tan g2rad nhist rfile > out");

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


  /* Input 2-gaussian fit parameters
   */
  if(!(fp=fopen(argv[2],"r")))
    kill("Error opening 2g-fit file tan.\n");
  for(i=1;i<=npts;++i)
    fscanf(fp,"%f %f %f %f %f",&fact[i],&x1,&sigt1[i],&x1,&sigt2[i]);
  fclose(fp);

  if(!(fp=fopen(argv[3],"r")))
    kill("Error opening 2g-fit file rad.\n");
  for(i=1;i<=npts;++i)
    fscanf(fp,"%f %f %f %f %f",&facr[i],&mu1[i],&sigr1[i],&mu2[i],&sigr2[i]);
  fclose(fp);

  /* Set up the spline interpolation routine for xi(r), v(r), sig(r).
   */
  spline(r,xi,npts,2.0E+30,2.0E+30,y2xi);
  spline(r,sigt1,npts,2.0E+30,2.0E+30,y2sigt1);
  spline(r,sigt2,npts,2.0E+30,2.0E+30,y2sigt2);
  spline(r,sigr1,npts,2.0E+30,2.0E+30,y2sigr1);
  spline(r,sigr2,npts,2.0E+30,2.0E+30,y2sigr2);
  spline(r,mu1,npts,2.0E+30,2.0E+30,y2mu1);
  spline(r,mu2,npts,2.0E+30,2.0E+30,y2mu2);
  
  fp1=fopen("g2dat","w");
  if(!(fp=fopen(argv[5],"r")))
    kill("Error opening rfile.\n");

  nhist=atoi(argv[4]);

  for(i=1;i<=nhist;++i)
    for(j=1;j<=nhist;++j)
      {
	sprintf(fname,"new-histo/h2d.%02d.%02d",i,j);
	fp2=fopen(fname,"w");
	fscanf(fp,"%f %f %f %f %f %f",&rsig,&rpi,&xlo,&ylo,&xhi,&yhi);
	/*fprintf(stderr,"%f %f %f %f %f %f\n",rsig,rpi,xlo,ylo,xhi,yhi);*/

	realrad=sqrt(rsig*rsig+rpi*rpi);
	splint(r,sigr1,y2sigr1,npts,realrad,&sr1);
	splint(r,sigr2,y2sigr2,npts,realrad,&sr2);
	splint(r,sigt1,y2sigt1,npts,realrad,&st1);
	splint(r,sigt2,y2sigt2,npts,realrad,&st2);
	splint(r,mu1,y2mu1,npts,realrad,&m1);
	splint(r,mu2,y2mu2,npts,realrad,&m2);
	
	realrad=sqrt(xlo*xlo+yhi*yhi);
	splint(r,sigr1,y2sigr1,npts,realrad,&sr1a);
	splint(r,sigr2,y2sigr2,npts,realrad,&sr2a);
	splint(r,sigt1,y2sigt1,npts,realrad,&st1a);
	splint(r,sigt2,y2sigt2,npts,realrad,&st2a);
	splint(r,mu1,y2mu1,npts,realrad,&m1a);
	splint(r,mu2,y2mu2,npts,realrad,&m2a);
	
	realrad=sqrt(xhi*xhi+ylo*ylo);
	splint(r,sigr1,y2sigr1,npts,realrad,&sr1b);
	splint(r,sigr2,y2sigr2,npts,realrad,&sr2b);
	splint(r,sigt1,y2sigt1,npts,realrad,&st1b);
	splint(r,sigt2,y2sigt2,npts,realrad,&st2b);
	splint(r,mu1,y2mu1,npts,realrad,&m1b);
	splint(r,mu2,y2mu2,npts,realrad,&m2b);
	
	fprintf(fp1,"%f %f %f %f %f %f %f %f %.1f %.1f\n",
		rsig,rpi,st1,st2,sr1,sr2,m1,m2,fact[1],facr[1]);
	for(k=-4000;k<=4000;k+=10)
	  {
	    x1=vel_prob2g(k,sr1,sr2,st1,st2,m1,m2,fact[1],facr[1],rsig,rpi);
	    x2=vel_prob2g(k,sr1a,sr2a,st1a,st2a,m1a,m2a,fact[1],facr[1],xlo,yhi);
	    x3=vel_prob2g(k,sr1b,sr2b,st1b,st2b,m1b,m2b,fact[1],facr[1],xhi,ylo);
	    if(isnan(x2))x2=0;
	    printf("%d %e %e %e\n",k,x1,x2,x3);
	  }
	for(k=-4000;k<=4000;k+=10)
	  fprintf(fp2,"%d %e\n",k,vel_prob2g(k,sr1,sr2,st1,st2,m1,m2,fact[1],facr[1],rsig,rpi));
	fclose(fp2);
      }

  exit(0);

}

float vel_prob2g(float x, float sigr1, float sigr2, float sigt1, float sigt2, float mu1, float mu2,
		 float ft, float fr, float rs, float rp)
{
  float x1,x2,x3,x4,costheta,sintheta,r;

  r=sqrt(rs*rs+rp*rp);
  costheta=rs/r;
  sintheta=rp/r;

  x1=ft*fr*term2(x,sigt1,sigr1,mu1,sintheta,costheta);
  x2=ft*(1-fr)*term2(x,sigt1,sigr2,mu2,sintheta,costheta);
  x3=(1-ft)*fr*term2(x,sigt2,sigr1,mu1,sintheta,costheta);
  x4=(1-ft)*(1-fr)*term2(x,sigt2,sigr2,mu2,sintheta,costheta);
  return(x1+x2+x3+x4);
}

/* This is the term that is repeated in vel_prob2g, just
 * with different paramters.
 */
float term2(double x, float st, float sr, float mu, float sintheta, float costheta)
{
  return(exp(-(mu*sintheta-x)*(mu*sintheta-x)/2/
	     (sintheta*sintheta*sr*sr+costheta*costheta*st*st))/
	 RT2PI/costheta/sr/st/
	 sqrt(1/sr/sr+sintheta*sintheta/costheta/costheta/st/st));
}



/* This is the asymmetric probability distribution, assuming that the
 * underlying pdfs are exponential.
 */
float vel_prob(float x, float sn, float sp, float st, float rs, float rp)
{
  double A_p,Ap_p,A_m,Ap_m,B_p,B_m,C_p,C_m,vcrit,r,p;

  rp=fabs(rp);
  r=sqrt(rs*rs+rp*rp);
  vcrit=r/rp*x;
  A_p=r/(rs*2*sp*st)*exp(-rt2*x*r/st/rs);
  Ap_p=r/(rs*2*sp*st)*exp(rt2*x*r/st/rs);
  B_p=1.0/sp + rp/(rs*st);
  C_p=1.0/sp - rp/(rs*st);

  A_m=r/(rs*2*sn*st)*exp(-rt2*x*r/st/rs);
  Ap_m=r/(rs*2*sn*st)*exp(rt2*x*r/st/rs);
  B_m=1.0/sn + rp/(rs*st);
  C_m=1.0/sn - rp/(rs*st);

  if(x<0)
    {
      p = A_m/(rt2*B_m)*exp(rt2*vcrit*B_m) + 
	Ap_m/(rt2*C_m)*(1.0 - exp(rt2*vcrit*C_m)) + 
	Ap_p/(rt2*B_p);
    }
  else
    {
      p = A_m/(rt2*B_m) + 
	A_p/(rt2*C_p)*(1.0 - exp(-rt2*vcrit*C_p)) + 
	Ap_p/(rt2*B_p)*exp(-rt2*vcrit*B_p);
    }
  if(isnan(p) || isinf(p))return(0);
  return(p);

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

