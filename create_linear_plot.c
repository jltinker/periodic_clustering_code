#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "nrutil.h"

void splie2(float x1a[], float x2a[], float **ya, int m, int n, float **y2a);
void splin2(float x1a[], float x2a[], float **ya, float **y2a, int m, int n,
	    float x1, float x2, float *y);


int main(int argc, char **argv)
{
  int i,j,nsize,nsize2,ix,iy,nbins,sm_out=0,NBODY_FORMAT=0,i1,i2,ng=20,j1;
  float x1,x2,x3,x4,x5,err;
  float r,rs,rp,dx1,dx2,dx3,dx4,**xi2d,*xx2d,*yy2d,**xi2d_data,**xi2d_kaiser,
    xi0_m,xi2_m,xi0_k,xi2_k,xi0_d,xi2_d,xi_r,**y2model,**y2data,rmax,
    rs1,rs2,rp1,rp2,xsum,dx,dy;
  FILE *fp,*fp2;
  char fname[100];

  
  if(argc<4)
    {
      fprintf(stderr,"> create_linear_plot xi2d.file rmax nbins [SM (1/0)] [NBODY FORMAT (1/0)> outfile\n");
      exit(0);
    }

  rmax=atof(argv[2]);
  nbins=atoi(argv[3]);
  if(argc>4)
    sm_out=atoi(argv[4]);
  if(argc>5)
    NBODY_FORMAT=atoi(argv[5]);

  j=-1;
  fp=fopen(argv[1],"r");
  j = filesize(fp);
  
  nsize=sqrt(j*1.0);
  nsize = atoi(argv[6]);
  nsize2 = atoi(argv[7]);
  fprintf(stderr,"Read %d x %d lines from file [%s] %d.\n",nsize,nsize2,argv[1],j);
  xi2d=matrix(1,nsize,1,nsize2);
  xx2d=vector(1,nsize);
  yy2d=vector(1,nsize2);
  xi2d_data=matrix(1,nsize,1,nsize2);
  xi2d_kaiser=matrix(1,nsize,1,nsize2);
  for(i=1;i<=nsize;++i)
    xx2d[i]=yy2d[i]=0;
  
  ix=iy=0;
  for(i=1;i<=j;++i)
    {
      iy++;
      if(i%nsize2==1)
	{
	  ix++;
	  iy=1;
	}
      switch(NBODY_FORMAT){
      case 1:
	fscanf(fp,"%f %f %f %f %f %f",&x2,&x1,&x3,&err,&x4,&x5);
	//printf("%d %d %f %f %f\n",ix,iy,x2,x1,x3);
	//fscanf(fp,"%f %f %f %f %f %f",&x1,&x5,&x4,&err,&x2,&x3);
	break;
      case 0:
	fscanf(fp,"%f %f %f %f %f %f %f",&x5,&x4,&x1,&x2,&x3,&err,&x5);
	break;
      case 2: //idit
	fscanf(fp,"%f %f %f %f %f %f",&x2,&x1,&x1,&x3,&x4,&x4);
	x1 = 2*(iy-0.5);
	x2 = pow(10.0,-1.09+0.2*(ix-0.5));
	break;
        case 3:
	fscanf(fp,"%f %f %f %f %f",&x3,&x2,&x1,&x4,&x4);
	break;
      }    
      /*
	x2 = rsig
	x1 = rpi
	x3 = xi2d

      if(NBODY_FORMAT)
	fscanf(fp,"%f %f %f %f %f %f",&x1,&x5,&x4,&err,&x2,&x3);
      else
	fscanf(fp,"%f %f %f %f %f %f %f",&x5,&x4,&x1,&x2,&x3,&err,&x5);
      */
      xi2d[ix][iy]=x3;
      xi2d_data[ix][iy]=x3;
      xx2d[ix]+=x2;
      yy2d[iy]+=x1;
      //printf("BOO %d %d %f %f %f\n",ix,iy,xx2d[ix],xx2d[iy],xi2d[ix][iy]);
    }      
  fclose(fp);
  for(i=1;i<=nsize;++i)
    {
      xx2d[i]/=(float)nsize2;
    }
  for(i=1;i<=nsize2;++i)
    {
      yy2d[i]/=(float)nsize;
    }

  for(i=1;i<=-nsize;++i)
    printf("%d %f\n",i,xx2d[i]);
  for(j=1;j<=-nsize2;++j)
    printf("%d %f\n",j,yy2d[j]);

  y2model=matrix(1,nsize,1,nsize2);
  y2data=matrix(1,nsize,1,nsize2);
  splie2(xx2d,yy2d,xi2d,nsize,nsize2,y2model);
  splie2(xx2d,yy2d,xi2d_data,nsize,nsize2,y2data);

  for(i=1;i<=nbins;++i)
    for(j=1;j<=nbins;++j)
      {
	if(sm_out)
	  {
	    rp2=rmax/(nbins*1.0)*i;
	    rs2=rmax/(nbins*1.0)*j;
	    rp1=rmax/(nbins*1.0)*(i-1);
	    rs1=rmax/(nbins*1.0)*(j-1);
	  }
	else
	  {
	    rs=rmax/(nbins*1.0)*i;
	    rp=rmax/(nbins*1.0)*j;
	  }
	dx=(rs2-rs1)/ng;
	dy=(rp2-rp1)/ng;
	xsum=0;
	for(i1=1;i1<=ng;++i1)
	  for(j1=1;j1<=ng;++j1)
	    {
	      rs=rs1+(i1-0.5)*dx;
	      rp=rp1+(j1-0.5)*dy;
	      splin2(xx2d,yy2d,xi2d,y2model,nsize,nsize2,rs,rp,&dx1);
	      xsum+=2*dy*rs*dx*dx1;
	    }
	xsum/=((rs2*rs2-rs1*rs1)*(rp2-rp1));
	
	/*
	splin2(xx2d,yy2d,xi2d,y2model,nsize,nsize,rs,rp,&dx1);
	splin2(xx2d,yy2d,xi2d_data,y2data,nsize,nsize,rs,rp,&dx2);
	*/
	printf("%e %e %e\n",0.5*(rs2+rs1),0.5*(rp2+rp1),xsum);
	fflush(stdout);
      }

}
