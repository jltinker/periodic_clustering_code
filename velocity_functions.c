#include <math.h>

/* This is a set of funtions to get the parameters for the two-gaussian
   fits, both as a function of halo mass and radius. Also, Omega_m and
   Sigma_8 will be input parameters.
   
   Input variables in all functions:

   - r - this is in h^{-1}Mpc.
   - m - this is in log_2(M/M_c) where M_c=706 particles (128 box), 45 particles (320 box)
   - omega - matter density
   - sig8 - sigma 8

   So for sigma2, the normalization is at halo mass log_2(m)=-4
*/

float sigma1rad(float r, float m, float omega, float sig8)
{
  float fac;
  r=log10(r);
  fac=(pow(omega,0.6)*sig8);
  return(fac*(165.0+294*r-174*r*r+143.0*r*r*r));
}

float mu2rad(float r, float m, float omega, float sig8)
{
  float fac;
  r=log10(r);
  fac=(pow(omega,0.6)*sig8);
  return(-fac*(296.0+1030*r-646*r*r));
}

float sigma2rad(float r, float m, float omega, float sig8)
{
  float fac,a,b,a0,b0,f;
  r=log10(r);
  fac=(pow(omega,0.6)*pow(sig8,4./3.));
  f=(630.0+618*r-526*r*r+195*r*r*r);

  a0=232-71.1*(-4);
  b0=398-41.6*(-4);

  a=232-71.1*(m);
  b=398-41.6*(m);

  return(fac*f/(a0+b0*r)*(a+b*r));
}

float mu1rad(float r, float m, float omega, float sig8)
{
  float fac,lgr;
  lgr=log10(r);
  fac=(pow(omega,0.6));
  
  return(-fac*(22.1+5.6*lgr+400*lgr*lgr-697*lgr*lgr*lgr+267*lgr*lgr*lgr*lgr + 
	 pow(1.0/r+1.0,3)*pow(2.0,m)*(0.25)));
}

float sigma1tan(float r, float m, float omega, float sig8)
{
  float fac;
  r=log10(r);
  fac=(pow(omega,0.6)*sig8);
  return(fac*(147.0+230*r-53*r*r+56*r*r*r));
}


float sigma2tan(float r, float m, float omega, float sig8)
{
  float fac,a,b,a0,b0,f;
  r=log10(r);
  fac=(pow(omega,0.6)*pow(sig8,4./3.));
  f=635.0+663*r-641*r*r+243*r*r*r;

  a0=232-71.1*(-4);
  b0=398-41.6*(-4);

  a=232-71.1*(m);
  b=398-41.6*(m);

  return(fac*f/(a0+b0*r)*(a+b*r));
}

void tangential_factor(float r, float omega, float sigma, float v, float *fac1, float *fac2)
{
  float slope, inter, vmin, smin, fac;

  fac=pow(omega,0.6)*sigma;
  slope=(6.5E-4*r - 1.18E-3)/fac;
  vmin=(351*r-42)*fac/sigma;
  if(v>vmin)slope*=-1;
  smin=0.5+0.3*r;
  
  inter=smin-slope*vmin;
  *fac1=(slope*v+inter);

  slope/=2;
  inter=smin-slope*vmin;  
  *fac2=(slope*v+inter);

}
