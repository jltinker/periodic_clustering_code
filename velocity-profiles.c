#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

#define NBIN 1000
#define pi 3.1415926535898


/**LOCAL VARIABLES**/
static int ***bincnt,**cnt;
static double **sig;
static float binwidth=10;

char OUTFILENAME[100];

/**LOCAL FUNCTIONS**/
void allocate_arrays(int n);
void output_angle_histograms(int ***n, int nr, int nphi, int nbin, float binwidth, char *root);
void output_radial_histograms(int **n, int nr, int nbin, float binwidth_a, char *root);
float velocity_components(float x1, float y1, float z1, 
			  float vx1, float vy1, float vz1, float x2, float y2, float z2, 
			  float vx2, float vy2, float vz2, float *vrad, float *vtan);

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
  char fname[100];

  fprintf(stderr,"HERE [%s] %d %f %d\n",a,nbin,radius[0],cnt[1][1]);
  fflush(stderr);


  for(i=1;i<=nbin;++i)
    for(j=1;j<=nbin;++j)
      {
	sig[i][j]=sqrt(sig[i][j]/cnt[i][j]);
	if(isnan(sig[i][j]))sig[i][j]=0;

	sprintf(fname,"h2d.%02d.%02d",i,j);
	fp=fopen(fname,"w");

	fac=(log10(cnt[i][j])*2-2);
	if(fac<1.5)fac=1.5;
	fac=sig[i][j]/fac;

	binfac=fac/binwidth;
	fprintf(stderr,"Binwidth for [%s] is %.0f %.1f %d\n",
		fname,binfac*binwidth,sig[i][j],cnt[i][j]);

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
  
  return;
 
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

/* This is for testing purposes. We will check how the l-o-s velocity disribution
 * changeswith angle for a constant radial seperation, plus to see if there is any
 * change in the radial and tangential distributions.
 */

void angular_dependence(float rsig, float rpi, float dvz, float x1, float y1, float z1, 
			float vx1, float vy1, float vz1, float x2, float y2, float z2, 
			float vx2, float vy2, float vz2, int rbin, char fname[])
{
  static double *rbar;
  static float **phibin,phiwidth;
  static int ***vzbin, **vrbin, **vtbin, flag=0,nr;
  int cnt, nbin=200,nphi=9,i,j,k,ibin,jbin,sum;
  float binwidth_a=20,r,rlo=21,rhi=40,phi,vrad,vtan,rcheck;
  float dx,dy,dz,rcube=128,rhalf=64;

  if(rbin<0)
    {
      nr=-rbin;
      fprintf(stderr,"Setting number of radial bins to [%d]\n",nr);
      return;
    }

  if(rsig<-100)
    {
      sprintf(OUTFILENAME,"%s",fname);
      for(k=1;k<=nr;++k)
	{
	  sum=0;
	  for(i=1;i<=nphi;++i)
	    for(j=-nbin;j<=nbin;++j)
	      sum+=vzbin[k][i][j];
	  rbar[k]/=(double)sum;
	  if(isnan(rbar[k]))rbar[k]=0;
	  printf("%d %f\n",k,rbar[k]);
	  fflush(stdout);
	}
    }

  /*
  if(rsig<-100)
    output_angle_histograms(vzbin,nr,nphi,nbin,binwidth_a,"zhist");
  */
  if(rsig<-100)
    output_radial_histograms(vrbin,nr,nbin,binwidth_a,"rhist");
  if(rsig<-100)
    output_radial_histograms(vtbin,nr,nbin,binwidth_a,"thist");

  if(!flag)
    {
      flag=1;
      phiwidth=90.0/nphi;
      rbar=dvector(1,nr);
      phibin=matrix(1,nr,1,nphi);
      vzbin=i3tensor(1,nr,1,nphi,-nbin,nbin);
      vrbin=imatrix(1,nr,-nbin,nbin);
      vtbin=imatrix(1,nr,-nbin,nbin);
      for(i=1;i<=nr;++i)
	rbar[i]=0;
      for(i=1;i<=nphi;++i)
	phibin[i]=0;
      for(i=1;i<=nr;++i)
	for(j=-nbin;j<=nbin;++j) 
	  vrbin[i][j]=vtbin[i][j]=0;
      for(i=1;i<=nr;++i)
	for(j=1;j<=nphi;++j)
	  for(k=-nbin;k<=nbin;++k)
	    vzbin[i][j][k]=0;
    }
  r=sqrt(rsig*rsig+rpi*rpi);
  rbar[rbin]+=r;


  /*if(r<rlo || r>rhi)return;*/
  rcheck=velocity_components(x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2,&vrad,&vtan);

  /*
  if(rbin==20 && vrad<-500)
    {
      printf("BOO %f %f %f %f %f %f %f\n",x1,y1,z1,x2,y2,z2,vrad);
      fflush(stdout);
    }
  */
  phi=fabs(atan(rpi/rsig))*180/pi;
  ibin=phi/phiwidth+1;
  if(ibin>nphi)ibin=nphi;

  jbin=dvz/binwidth_a;
  if(dvz<0)jbin--;
  if(abs(jbin)<=nbin)
    {
      vzbin[rbin][ibin][jbin]++;
    }

  jbin=vrad/binwidth_a;
  if(vrad<0)jbin--;
  if(abs(jbin)<=nbin)
    vrbin[rbin][jbin]++;

  jbin=vtan/binwidth_a;
  if(vtan<0)jbin--;
  if(abs(jbin)<=nbin)
    vtbin[rbin][jbin]++;
}

float velocity_components(float x1, float y1, float z1, 
			  float vx1, float vy1, float vz1, float x2, float y2, float z2, 
			  float vx2, float vy2, float vz2, float *vrad, float *vtan)
{
  float rcube=324,rhalf=162,dx,dy,dz,v12,r,dvx,dvy,dvz,ex,ey,ez,rcheck,vtot2,vt2,vtan2;
  static int flag=0;

  if(!flag++)
    srand48(-444432232);

  dx=x2-x1 ;
  dy=y2-y1 ;
  dz=z2-z1 ;
  if (dx>rhalf)  dx -= rcube ;
  if (dx<-rhalf) dx += rcube ;
  if (dy>rhalf)  dy -= rcube ;
  if (dy<-rhalf) dy += rcube ;
  if (dz>rhalf)  dz -= rcube ;
  if (dz<-rhalf) dz += rcube ;		
  r=sqrt(dx*dx+dy*dy+dz*dz);
  rcheck=r;

  dvx=vx2-vx1 ;
  dvy=vy2-vy1 ;
  dvz=vz2-vz1 ;
  vtot2=dvx*dvx+dvy*dvy+dvz*dvz;

  *vrad=v12=(dvx*dx+dvy*dy+dvz*dz)/r ;
  vtan2=dvx*dvx+dvy*dvy+dvz*dvz-v12*v12;

  dx/=r;
  dy/=r;
  dz/=r;
  if(dz!=0)
    {
      ex=drand48()-0.5;
      ey=drand48()-0.5;
      ez=-(dx*ex + dy*ey)/dz;
      r=sqrt(ex*ex+ey*ey+ez*ez);
    }
  else
    {
      if(dy!=0)
	{
	  ex=drand48()-0.5;
	  ez=drand48()-0.5;
	  ey=-(dx*ex + dz*ez)/dy;
	  r=sqrt(ex*ex+ey*ey+ez*ez);
	}
      else
	{
	  ez=drand48()-0.5;
	  ey=drand48()-0.5;
	  ex=-(dz*ez + dy*ey)/dx;
	  r=sqrt(ex*ex+ey*ey+ez*ez);
	}
    }
  ex/=r;
  ey/=r;
  ez/=r;

  /* Currently I think this line is wrong.
   */
  *vtan=(dvx-dvx*dx)*ex + (dvy-dvy*dy)*ey + (dvz-dvz*dz)*ez;

  /* Evidence suggests this line is correct
   */
  *vtan=dvx*ex+dvy*ey+dvz*ez;

  *vrad=v12;
  return(rcheck);

}

void output_angle_histograms(int ***n, int nr, int nphi, int nbin, float binwidth_a, char *root)
{
  int i,j,k,cnt;
  char fname[100],root2[100],root3[100];
  FILE *fp;

  /* If the output filename has not been passed via the command
   * line, use the name in the text file "filename".
   */
  if(strncmp(OUTFILENAME,"999",3)==0)
    {
      if(fp=fopen("filename","r"))
	{
	  fscanf(fp,"%s",root2);
	  sprintf(root3,"%s.%s",root2,root);
	  fclose(fp);
	}
      else
	sprintf(root3,"%s",root);
    }
  else
    {
      sprintf(root3,"%s.%s",OUTFILENAME,root);
    }
  for(i=1;i<=nr;++i)
    for(j=1;j<=nphi;++j)
      {
	sprintf(fname,"%s%d.%d",root3,i,j);
	fp=fopen(fname,"w");
	for(cnt=0,k=-nbin;k<=nbin;++k)
	  cnt+=n[i][j][k];
	
	fprintf(stderr,"opening [%s] %d %.1f\n",fname,cnt,binwidth_a);
	if(!cnt)cnt=1;
	for(k=-nbin;k<=nbin;++k)
	  fprintf(fp,"%f %e\n",(k+0.5)*binwidth_a,(float)n[i][j][k]/binwidth_a/cnt);
	/*
	  for(j=-nbin;j<=nbin;++j)
	  fprintf(fp,"%f %d\n",(j+0.5)*binwidth_a,n[i][j]);
	*/
	fclose(fp);
      }	  

}

void output_radial_histograms(int **n, int nr, int nbin, float binwidth_a, char *root)
{
  int i,j,k,cnt;
  char fname[100],root2[100],root3[100];
  FILE *fp;

  /* If the output filename has not been passed via the command
   * line, use the name in the text file "filename".
   */
  if(strncmp(OUTFILENAME,"999",3)==0)
    {
      if(fp=fopen("filename","r"))
	{
	  fscanf(fp,"%s",root2);
	  sprintf(root3,"%s.%s",root2,root);
	  fclose(fp);
	}
      else
	sprintf(root3,"%s",root);
    }
  else
    {
      sprintf(root3,"%s.%s",OUTFILENAME,root);
    }
  for(i=1;i<=nr;++i)
    {
      sprintf(fname,"%s.%d",root3,i);
      fp=fopen(fname,"w");
      for(cnt=0,k=-nbin;k<=nbin;++k)
	cnt+=n[i][k];

      fprintf(stderr,"opening [%s] %d %.1f\n",fname,cnt,binwidth_a);
      if(cnt==0)cnt=1;
      for(k=-nbin;k<=nbin;++k)
	fprintf(fp,"%f %e\n",(k+0.5)*binwidth_a,(float)n[i][k]/binwidth_a/cnt);
      fclose(fp);
    }	  

}

