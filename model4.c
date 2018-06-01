#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

/** Same as model.c, but this time we're using the asymmetric probability distribution
 *  assuming that the underlying distributions are exponential and that the mode is at
 *  zero velocity (rather than at the mean).
 *
 * NEW--> Now we have implemented the analytic version of the 2-gaussian fit. We
 * will be using that from now on.
 */

#define rt2  1.414214
#define NMAX 70
#define N2d  25
#define pi 3.1415926535898
#define RT2PI 2.506628

/** NUMERICAL RECIPES **/
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
float qromo(float (*func)(float), float a, float b,
	    float (*choose)(float(*)(float), float, float, int));
float midpnt(float (*func)(float), float a, float b, int n);
float trapzd(float (*func)(float), float a, float b, int n);


/** INTERNAL ROUTINES **/
float exp1(float x, float s);
float gauss1(float x, float s);
float interpolate2d(float x, float y, float **xi, int n, float rmax);
float ngp2d(float x, float y, float **xi, int n, float rmax);
float vel_prob(float x, float sn, float sp, float st, float rs, float rp);
float term2(float v, float st, float sr, float mu, float sintheta, float costheta);
float vel_prob2g(float x, float sigr1, float sigr2, float sigt1, float sigt2, float mu1, float mu2,
		 float ft, float fr, float rs, float rp);
float internal_bin_sum(float **xi2d, int n2d, float rmax, float xhi, float yhi);

float integrate_v_los_prob(float x, float sigr1, float sigr2, float sigt1, float sigt2, 
			   float mu1, float mu2, float ft, float fr, float rs, float rp);
float func2(float v);
float func2_exp(float v);

/** GLOBALS FOR USE WITH QROMO **/
float g_logr,
  g_sintheta,
  g_costheta,
  g_vlos1,
  g_sigt1,
  g_sigt2,
  g_sigr1,
  g_sigr2,
  g_mu1,
  g_mu2,
  g_fr,
  g_ft;

/** COSMO FACTORS **/
float OMEGA_M=0.3,
  SIGMA8=0.8;

int main(int argc, char **argv)
{
  float *v,*r,*xi,*sig,**xi2d,rsig,rpi,y,realrad,*y2xi,*y2sig,*y2v,x1,x2,x3,x4,x5,x6,
    v1,s1,dv,dy,rmax=0,ymin,ymax,*sign,*sigp,*sigtan,*y2sign,*y2sigp,*y2sigtan,s1tan,s1p,s1n,
    *sigr1,*sigr2,*sigt1,*sigt2,*mu1,*mu2,*y2sigr1,*y2sigr2,*y2sigt1,*y2sigt2,*y2mu1,*y2mu2,
    sr1,sr2,st1,st2,m1,m2,*fact,*facr,err,xup,yup,chi2,xx,yy,*y2f1,*y2f2,f1,f2;
  int i=0,j,k,npts,ncol,Ny,i1,i2,gridsize;
  FILE *fp,*fp2;

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
  y2f1=vector(1,NMAX);
  y2f2=vector(1,NMAX);
  y2xi=vector(1,NMAX);
  y2v=vector(1,NMAX);

  if(argc<5)
    kill("model4 xifile xi2dfile g2tan g2rad > out");

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
  fprintf(stderr,"RMAX %f\n",rmax);



  /* Input 2-gaussian fit parameters
   */
  if(!(fp=fopen(argv[3],"r")))
    kill("Error opening 2g-fit file tan.\n");
  for(i=1;i<=npts;++i)
    fscanf(fp,"%f %f %f %f %f",&fact[i],&x1,&sigt1[i],&x1,&sigt2[i]);
  fclose(fp);

  if(!(fp=fopen(argv[4],"r")))
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
  spline(r,fact,npts,2.0E+30,2.0E+30,y2f1);
  spline(r,facr,npts,2.0E+30,2.0E+30,y2f2);
  

  /* Integrate in linear space, then convert to log space (later).
   */
  fprintf(stderr,"  ");
  for(i=1;i<=N2d;++i)
    for(j=1;j<=N2d;++j)
      {
	if(j==1)fprintf(stderr,"\b\b%02d",i);	
	rsig=(i-0.5)/N2d*rmax;
	rpi=(j-0.5)/N2d*rmax;
	ymin=r[1];
	ymax=sqrt(2*rmax*rmax-rsig*rsig);
 	Ny=ymax/ymin*2.5;
	dy=ymax/(Ny);

	rsig=0.668478;
	rpi=0.463007;

	for(xi2d[i][j]=0,k=-Ny+1;k<=Ny;++k)
	  {
	    y=(k-0.5)/Ny*ymax;
	    y=-5;
	    realrad=sqrt(y*y+rsig*rsig);	
	    dv=(rpi-y)*100.0;

	    splint(r,xi,y2xi,npts,realrad,&x1);
	    splint(r,sigr1,y2sigr1,npts,realrad,&sr1);
	    splint(r,sigr2,y2sigr2,npts,realrad,&sr2);
	    splint(r,sigt1,y2sigt1,npts,realrad,&st1);
	    splint(r,sigt2,y2sigt2,npts,realrad,&st2);
	    splint(r,mu1,y2mu1,npts,realrad,&m1);
	    splint(r,mu2,y2mu2,npts,realrad,&m2);
	    splint(r,fact,y2f1,npts,realrad,&f1);
	    splint(r,facr,y2f2,npts,realrad,&f2);

	    x2=integrate_v_los_prob(dv,sr1,sr2,st1,st2,m1,m2,f1,f2,rsig,y);
	    exit(0);
	    /*
	    x2=vel_prob2g(dv,sr1,sr2,st1,st2,m1,m2,f1,f2,rsig,y);
	    */
	    if(x2<0)
	      {
		fprintf(stderr,"BOOO %f %f %f %f\n",x2,y,rsig,rpi);
		x2=0;
	      }
	    xi2d[i][j]+=(1+x1)*x2*dy*100;

	    printf("%f %f %f %f %f %e %f %f %f %f\n",
		   realrad,rsig,rpi,y,dv,x2,dy,x1,xi2d[i][j],x3);

	    /*
	    if(k==1 && i==32)
	      printf("PROB %f %f %f %f %f %e %f %f %f %f\n",
		     realrad,rsig,rpi,y,dv,x2,dy,x1,xi2d[i][j],x3);

	    if(j==2 && i==32)
	      printf("G2 %f %f %f %f %f %f %f %f %f\n",realrad,y,dv,m1,sr1,m2,sr2,st1,st2);

	    splint(r,xi,y2xi,npts,sqrt(rsig*rsig+rpi*rpi),&x3);
	    if(j==20 && i==32)
	      printf("VEL1 %f %f %f %f %f %e %f %f %f %f\n",
		     realrad,rsig,rpi,y,dv,x2,dy,x1,xi2d[i][j],x3);
	    if(j==2 && i==32)
	      printf("VEL2 %f %f %f %f %f %e %f %f %f %f\n",
		     realrad,rsig,rpi,y,dv,x2,dy,x1,xi2d[i][j],x3);
	    */
	    fflush(stdout);
	  }
	fprintf(stderr,"SUM= %e\n",xi2d[i][j]);
	exit(0);
      }
  fprintf(stderr,"\b\b\n");

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
      fscanf(fp,"%f %f %f %f %f %f",&x1,&rsig,&rpi,&err,&xup,&yup);
      i++;
    }
  rewind(fp);

  gridsize=sqrt(i-1);

  chi2=0;
  xx=0;
  yy=0;
  for(j=1;j<i;++j)
    {      
      fscanf(fp,"%f %f %f %f %f %f",&x1,&rsig,&rpi,&err,&xup,&yup);
      err/=sqrt(14.0);
      x2=internal_bin_sum(xi2d,N2d,rmax,xup,yup);      
      
      if(x2==0)x2=interpolate2d(rsig,rpi,xi2d,N2d,rmax);

      /*x2=interpolate2d(xx,yy,xi2d,N2d,rmax);*/
      if(j%10==0)xx=xup;
      yy=yup;

      splint(r,xi,y2xi,npts,sqrt(rpi*rpi+rsig*rsig),&x3);
      printf("%f %f %f %f %f %e\n",x1,rsig,rpi,x2,x3,err);
      chi2+=(x2-x1)*(x2-x1)/err/err;

      fflush(stdout);
    }
  fprintf(stderr,"CHI2: %e\n",chi2);
  

  /* Now print out the perpendicular stuff
   */
  rewind(fp);
  fp2=fopen("out2","w");
  for(i=1;i<=gridsize;++i)
    {
      for(j=1;j<=gridsize*gridsize;++j)
	{
	  fscanf(fp,"%f %f %f %f %f %f",&x1,&rsig,&rpi,&err,&xup,&yup);
	  if(j%10!=i%gridsize)continue;
 
	  x2=internal_bin_sum(xi2d,N2d,rmax,xup,yup);      
	  
	  if(x2==0)x2=interpolate2d(rsig,rpi,xi2d,N2d,rmax);
	  
	  splint(r,xi,y2xi,npts,sqrt(rpi*rpi+rsig*rsig),&x3);
	  fprintf(fp2,"%f %f %f %f %f %e\n",x1,rsig,rpi,x2,x3,err/sqrt(14.));
	  fflush(fp2);
	}
      rewind(fp);
    }

  

  exit(0);

  for(i=1;i<=N2d;++i)
    for(j=1;j<=1;++j)
      fprintf(stdout,"%e %e %f %f\n",xi2d[i][j],xi2d[i][j+1],(i-0.5)/N2d*rmax,(j-0.5)/N2d*rmax);
  exit(0);

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
float term2(float x, float st, float sr, float mu, float sintheta, float costheta)
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
  float A_p,Ap_p,A_m,Ap_m,B_p,B_m,C_p,C_m,vcrit,r,p;

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

/* This is the nearest grid point (ngp) method of getting
 * the model value of xi2d
 */
float ngp2d(float x, float y, float **xi, int n, float rmax)
{
  int i,j;

  i=x/rmax*n+1;
  j=y/rmax*n+1;
  if(i<1)i=1;
  if(j<1)j=1;
  printf("BOO %f %f %d %d %f %d %f\n",x,y,i,j,n,rmax,xi[i][j]);
  return(xi[i][j]);
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
      fprintf(stderr,"%f %f %f\n",a,b,rmax);
      kill("rmax too small for input xi2d.");
    }

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
  if(ip==n)
    {
      ip--;
      ip1--;
    }
  if(jp==n)
    {
      jp--;
      jp1--;
    }


  t=x-(ip+0.5);
  u=y-(jp+0.5);
  
  t1=1-t;
  u1=1-u;
  
  x1=t1*u1*xi[ip][jp] + 
    t*u1*xi[ip1][jp] +
    t1*u*xi[ip][jp1] + 
    t*u*xi[ip1][jp1];
  /*
  printf("INT %f %f %f %f %d %d %f %f %f %f %f\n",x,y,t,u,ip,jp,x1,xi[ip][jp],
	 xi[ip1][jp1],xi[ip][jp1],xi[ip1][jp1]);
  */
  return(x1);
}


float internal_bin_sum(float **xi2d, int n2d, float rmax, float xhi, float yhi)
{
  float dx,dy,sum=0,wsum=0,x,y,x2;
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
	sum+=x;
	wsum+=x*(1+x2);
      }
  xprev=xhi;
  ylo=yhi;
  return(wsum/sum-1);

}

float integrate_v_los_prob(float x, float sigr1, float sigr2, float sigt1, float sigt2, 
			   float mu1, float mu2, float ft, float fr, float rs, float rp)
{
  float s1,r,error=1,tol=1.0E-4,prev;
  int n=1;

  /* Gotta put all this shit into global variables
   */
  r=sqrt(rp*rp+rs*rs);
  g_logr=log10(r);
  g_sintheta=rp/r;
  g_costheta=rs/r;
  g_vlos1=x;
  g_sigt1=sigt1;
  g_sigt2=sigt2;
  g_sigr1=sigr1;
  g_sigr2=sigr2;
  g_mu1=mu1;
  g_mu2=mu2;
  g_fr=fr;
  g_ft=ft;

  /*
  fprintf(stderr,"starting %f %f %f\n",rp,rs,x);
  s1=qromo(func2,-2000,2000,midpnt); 
  fprintf(stderr,"finished %f %e\n",rp,s1/g_costheta);
  return(s1/g_costheta);
  */

  /* Let's try a different integration routine
   */
  prev=trapzd(func2_exp,-2000,2000,n);
  while(error>tol) {
    n++;
    s1=trapzd(func2_exp,-2000,2000,n);
    error=fabs(s1-prev)/prev;
    if(n>5 && s1==0)break;
    if(n>7 && s1<1.0E-10)break;
    /*fprintf(stderr,"%d %e %e %e\n",n,s1,prev,error);*/
    prev=s1;
  }
  /*fprintf(stderr,"finished %d %e %f %e\n",n,error,rp,s1/g_costheta);*/
  return(s1/g_costheta);

}

/* Now I will put the double-gaussian with vrad-dependent sigma_t
 * in this program, but I'd really rather just integrate that numerically.
 */
float func2(float v)
{
  float a1,a2,b1,b2,sigt1,sigt2,slope,vmin,smin,inter;

  vmin=11.3+36.8*g_logr+367.4*g_logr*g_logr;
  /*vmin=-131.4+430.4*g_logr;*/
  slope=-.00157+0.000904*g_logr;
  if(v>vmin)slope*=-1;
  smin=.387+.34*g_logr;
  inter=smin-slope*vmin;
  sigt2=(slope*v+inter)*g_sigt2;

  if(v>vmin)slope*=-1;
  vmin=-51.4+430.4*g_logr;
  if(v>vmin)slope*=-1;
  slope/=2;
  inter=smin-slope*vmin;  
  sigt1=(slope*v+inter)*g_sigt1;

  if(g_logr<0.0)
    {
      sigt1=g_sigt1;
      sigt2=g_sigt2;
    }
  a1=g_ft*exp(-(g_vlos1-g_sintheta*v)*(g_vlos1-g_sintheta*v)/
	    g_costheta/g_costheta/2/sigt1/sigt1)/RT2PI/sigt1;
  a2=(1-g_ft)*exp(-(g_vlos1-g_sintheta*v)*(g_vlos1-g_sintheta*v)/
		g_costheta/g_costheta/2/sigt2/sigt2)/RT2PI/sigt2;
  b1=g_fr*exp(-(v-g_mu1)*(v-g_mu1)/2/g_sigr1/g_sigr1)/RT2PI/g_sigr1;
  b2=(1-g_fr)*exp(-(v-g_mu2)*(v-g_mu2)/2/g_sigr2/g_sigr2)/RT2PI/g_sigr2;
  return((a1+a2)*(b1+b2));
}
  
/* Now I will put the double-gaussian with vrad-dependent sigma_t
 * in this program, but I'd really rather just integrate that numerically.
 * NEW-> Been modified such that the non-linear component of the tan-pdf is
 * an exponential.
 */
float func2_exp(float v)
{
  float a1,a2,b1,b2,sigt1,sigt2,slope,vmin,smin,inter,a[7],vrad1,vt;

  a[1]=0.634475;
  a[2]=0.191846;
  a[3]=-0.348138;
  a[4]=-0.097268;
  a[5]=0.000754;
  a[6]=-0.000409;

  vrad1=v/(pow(OMEGA_M,0.6)*SIGMA8);
  smin=a[1]+a[2]*g_logr;
  vmin=a[3]+a[4]*g_logr;
  slope=a[5]+a[6]*g_logr;
  if(vrad1<vmin)slope*=-1;
  sigt1=g_sigt1*(slope*(vrad1-vmin)+smin);
  sigt2=g_sigt2*(slope*(vrad1-vmin)+smin);

  a1=g_ft*exp(-fabs(g_vlos1-g_sintheta*v)*rt2/
	      g_costheta/sigt1)/(rt2*sigt1);
  a2=(1-g_ft)*exp(-(g_vlos1-g_sintheta*v)*(g_vlos1-g_sintheta*v)/
		  g_costheta/g_costheta/2/sigt2/sigt2)/RT2PI/sigt2;
  b1=g_fr*exp(-(v-g_mu1)*(v-g_mu1)/2/g_sigr1/g_sigr1)/RT2PI/g_sigr1;
  b2=(1-g_fr)*exp(-(v-g_mu2)*(v-g_mu2)/2/g_sigr2/g_sigr2)/RT2PI/g_sigr2;

  vt=(g_vlos1-g_sintheta*v)/g_costheta;
  printf("%f %f %e %e %e %e %f %f\n",v,vt,a1,a2,b1,b2,g_sintheta,g_costheta);
  return((a1+a2)*(b1+b2));
}
  
