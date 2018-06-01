#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "nrutil.h"

void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
float qtrap(float (*func)(float), float a, float b);

void average_files(int ipow, int fnumber);
float avg_xi_m(float r);
float ***xi_m,***rad_m,***y2_m;
int K,J,LINES,NRUNS=5;
char PATH[100],OUTPATH[100];

int main(int argc, char **argv)
{
  int i,fn,plo,phi,j,outlo,outhi;
  
  if(argc<4)
    {
      fprintf(stderr,"xi-multipoles outlo outhi pow2low pow2hi [NRUNS] [INPATH] [OUTPATH]\n\n");
      exit(0);
    }

  outlo=atoi(argv[1]);
  outhi=atoi(argv[2]);
  plo=atoi(argv[3]);
  phi=atoi(argv[4]);
  if(argc>5)
    NRUNS=atoi(argv[5]);
  if(argc>6)
    sprintf(PATH,"%s",argv[6]);
  else
    sprintf(PATH,"run");
  if(argc>7)
    sprintf(OUTPATH,"%s",argv[7]);
  else
    sprintf(OUTPATH,".");
  
  for(i=outlo;i<=outhi;++i)
    for(j=plo;j<=phi;++j)
      average_files(j,i);

  return(0);
}

void average_files(int ipow, int fnumber)
{
  int n,j,i,ix,lines=0,k,nr,power,nsig,npi,nxi,r_lines=0;
  FILE *fp[5][3],*outf,*fpr[5];
  char filename[200],mstar[10];
  float x1,x2,x3,x4,x5,rsig,rpi,esig,epi,rsig2,rpi2,xi,xi2,exi,xup,yup,
    xm,xm2,xq,xq2,em,eq,r,xx[15],**xi_r,**rad,**xi_y2,xm_hi,xq_hi,xm_lo,xq_lo,
    r_temp,xm_temp,qp_temp;

  /* THIS IS THE OLD FILE NOTATION
   *
   *
  power=pow(2.0,abs(ipow));
  if(ipow<0)
    sprintf(mstar,"mstar%d",power);
  if(ipow>=0)
    sprintf(mstar,"%dmstar",power);
  */

  /* New file notation:
   */
  sprintf(mstar,"g%d",ipow);

  for(i=0;i<NRUNS;++i)
    for(j=0;j<3;++j)
      {
	if(NRUNS>1)
	  sprintf(filename,"%s%d/xi_ang.%02d.%s.%d",PATH,i+1,fnumber,mstar,j+1);
	else
	  sprintf(filename,"%s/xi_ang.%02d.%s.%d",PATH,fnumber,mstar,j+1);
	if(!(fp[i][j]=fopen(filename,"r")))
	  {
	    fprintf(stderr,"ERROR with file %d %d [%s]\n",i+1,j+1,filename);
	    return;
	  }
      }
  
  for(i=0;i<NRUNS;++i)
    {
      if(NRUNS>1)
	sprintf(filename,"%s%d/xi.%02d.%s",PATH,i+1,fnumber,mstar);
      else
	sprintf(filename,"%s/xi.%02d.%s",PATH,fnumber,mstar);	
      if(!(fpr[i]=fopen(filename,"r")))
	{
	  fprintf(stderr,"ERROR with file %d [%s]\n",i+1,filename);
	  return;
	}
    }
  
  sprintf(filename,"%s/xi_m.%02d.%s.avg",OUTPATH,fnumber,mstar);
  if(!(outf=fopen(filename,"w")))
    {
      fprintf(stderr,"ERROR opening outfile: [%s]\n",filename);
      exit(0);
    }

  while(!feof(fp[0][0]))
    {
      lines++;
      fgets(filename,100,fp[0][0]);
    }
  lines--;
  LINES=lines;
  fprintf(stderr,"z_lines: %d\n",lines);

  rewind(fp[0][0]);    

  while(!feof(fpr[0]))
    {
      r_lines++;
      fgets(filename,200,fpr[0]);
    }
  r_lines-=2;

  fprintf(stderr,"r_lines: %d\n",r_lines);

  rewind(fpr[0]);    

  /* Before reading in the multipole information, read in the real-space xi(r)
   * to interpolate the values of xi(r) at the centers of the multipole bins.
   */
  xi_r=matrix(0,NRUNS-1,1,r_lines);
  rad=matrix(0,NRUNS-1,1,r_lines);
  xi_y2=matrix(0,NRUNS-1,1,r_lines);

  for(j=0;j<NRUNS;++j)
    fgets(filename,200,fpr[j]);

  for(i=1;i<=r_lines;++i)
    for(j=0;j<NRUNS;++j)
      {
	for(k=0;k<13;++k)
	  fscanf(fpr[j],"%f",&xx[k]);
	fscanf(fpr[j],"%d",&ix);
	xi_r[j][i]=xx[3];
	rad[j][i]=xx[2];
      }

  for(j=0;j<NRUNS;++j)
    spline(rad[j],xi_r[j],r_lines,1.0E+30,1.0E+30,xi_y2[j]);

  /* For convenience in calculating the spherically averaged monopole,
   * tabulate the monopole values.
   */
  xi_m=f3tensor(0,NRUNS-1,0,2,1,lines);
  rad_m=f3tensor(0,NRUNS-1,0,2,1,lines);
  y2_m=f3tensor(0,NRUNS-1,0,2,1,lines);
 
  for(i=1;i<=lines;++i)
    for(j=0;j<NRUNS;++j)
      for(k=0;k<3;++k)
	{
	  fscanf(fp[j][k],"%f %f %f",&x1,&x2,&x3);
	  xi_m[j][k][i]=(x2);
	  rad_m[j][k][i]=(x1);
	}
  for(j=0;j<NRUNS;++j)
    for(k=0;k<3;++k)
      spline(rad_m[j][k],xi_m[j][k],lines,1.0E+30,1.0E+30,y2_m[j][k]);

  for(j=0;j<NRUNS;++j)
    for(k=0;k<3;++k)
      rewind(fp[j][k]);

  /* Now actually read in the values and caluclate:
   * 1) Monopole to real-space.
   * 2) Quadrupole to (avg. monopole - monopole)
   */

  for(i=1;i<=lines;++i)
    {
      r=0;
      xm=xm2=0;
      xq=xq2=0;
      nr=0;
      xm_hi=xq_hi=-100;
      xm_lo=xq_lo=100;
      for(j=0;j<NRUNS;++j)
	{
	  r_temp=xm_temp=qp_temp=0;
	  for(k=0;k<3;++k)
	    {
	      J=j;
	      K=k;
	      fscanf(fp[j][k],"%f %f %f",&x1,&x2,&x3);
	      splint(rad[j],xi_r[j],xi_y2[j],r_lines,x1,&x4);
	      if(x1>rad_m[j][k][1])
		x5=qtrap(avg_xi_m,(rad_m[j][k][1]),(x1))*3.0/(x1*x1*x1);
	      else
		x5=0;
	      
	      r_temp+=x1/3;
	      xm_temp+=x2/x4/3;
	      qp_temp+=x3/(x2-x5)/3;
	    }
	  nr++;
	  r+=r_temp;
	  xm+=xm_temp;
	  xm2+=xm_temp*xm_temp;
	  xq+=qp_temp;
	  xq2+=qp_temp*qp_temp;
	}
      r/=nr;

      xm/=nr;
      xm2/=nr;
      xq/=nr;
      xq2/=nr;

      if(nr>1)
	{
	  em=sqrt((xm2-xm*xm)/(nr-1.0));
	  eq=sqrt((xq2-xq*xq)/(nr-1.0));
	}
      else
	em=eq=0;

      fprintf(outf,"%f %e %e %e %e\n",r,xm,em,xq,eq);
    }
  fclose(outf);
  for(i=0;i<NRUNS;++i)
    for(j=0;j<3;++j)
      fclose(fp[i][j]);
  for(i=0;i<NRUNS;++i)
    fclose(fpr[i]);
}

float avg_xi_m(float r)
{
  float a;
  splint(rad_m[J][K],xi_m[J][K],y2_m[J][K],LINES,r,&a);
  return(r*r*a);
}
