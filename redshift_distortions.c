#include <stdio.h>
#include <math.h>
#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define pi 3.1415926535898

/*
 * This function applies redshift-space distortions to the dimension passed.
 * The velocity accross the box is rcube*100 [km/s], so the distortion in a
 * particle's position is:
 *
 *     dz = vz/vcube*rcube = vz/100
 *
 * We're adjusting the routine to take out the redshift dependence of the 
 * velocities v(proper)=v(z)/(1+z)
 */

/** NUMERICAL RECIPES ROUTINE **/
float gasdev(long *idum);

void redshift_distortions(float *z, float *v, float rcube, int np, int imode)
{
  int i,idum=55;
  float dv,zp1,fac;

  switch(imode) {
    case 1: zp1=1.0; fac=pow(0.1,-0.6); break;
    case 2: zp1=1.2; fac=pow(0.161,-0.6); break;
    case 3: zp1=1.55; fac=pow(0.295,-0.6); break;
    case 4: zp1=2.0; fac=pow(0.471,-0.6); break;
    case 5: zp1=2.47; fac=pow(0.636,-0.6); break;
  default: zp1=1.0; fac=1;
  }
  
  /* TEMP--> normalizing everything to Omega=1
   */
  /*zp1/=fac;*/
  fprintf(stderr,"Dividing each velocity by: %f\n",zp1);
  
  for(i=0;i<np;++i)
    {
      z[i]+=v[i]/100.0/zp1;
      v[i]/=zp1;
      if(z[i]<0.0)z[i]+=rcube;
      if(z[i]>=rcube)z[i]-=rcube;
    }
  return;

  /* Temporary distraction. Testing the model.
   */
  for(i=0;i<np;++i)
    {
      dv=400.0*gasdev(&idum);
      z[i]+=dv/100.0;
      if(z[i]<0.0)z[i]+=rcube;
      if(z[i]>=rcube)z[i]-=rcube;
    }
  return;

}
