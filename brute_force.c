#include <stdio.h>
#include <math.h>
#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define pi 3.1415926535898


void brute_force_method(int p, float *x, float *y, float *z, int npair[101][101], int binlookup[], 
			float binfac, int np, int nrbin, float zmax)
{
  int i,j,k,kbin,rbin;
  float dx,dy,dz,r,rcube,halfcube;
  static int ntot=0,im,pm;
  static float rmax=0,mdx,mdy;

  rcube=128;
  halfcube=rcube*0.5;

  for(i=p+1;i<np;++i)
    {
      dx=fabs(x[i]-x[p]);
      if(dx>halfcube)dx=rcube-dx;
      
      /*
      dy=fabs(y[i]-y[p]);
      if(dy>halfcube)dy=rcube-dy;
      r=sqrt(dx*dx+dy*dy) ;
      */
      r=dx;

      if(r>=zmax)continue;
      kbin=binlookup[(int)(binfac*r)] ;

      dz=mabs(z[p]-z[i]);
      if(dz>halfcube)dz=rcube-dz;
      if(dz>=zmax)continue;

      rbin=binlookup[(int)(binfac*dz)];

      npair[kbin][rbin]++;
    }
  /*
  if(p==np-1)
    printf("NTOT %d %f %.1f %.1f %.1f %.1f %.1f %.1f %d %d %f\n",ntot,rmax,x[im],y[im],x[pm],y[pm],
	   mdx,mdy,im,pm,halfcube);
  */
}
