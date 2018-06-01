#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void average_files(int ipow, int fnumber);

int main(int argc, char **argv)
{
  int i,fn,plo,phi,j,outlo,outhi;
  
  if(argc<4)
    {
      fprintf(stderr,"avg-files outlo outhi pow2low pow2hi\n\n");
      exit(0);
    }

  outlo=atoi(argv[1]);
  outhi=atoi(argv[2]);
  plo=atoi(argv[3]);
  phi=atoi(argv[4]);

  for(i=outlo;i<=outhi;++i)
    for(j=plo;j<=phi;++j)
      average_files(j,i);

  return(0);
}

void average_files(int ipow, int fnumber)
{
  int n,j,i,ix,lines=0,k,nr,power,nsig,npi,nxi;
  FILE *fp[5][3],*outf;
  char filename[100],mstar[10];
  float x1,x2,x3,rsig,rpi,esig,epi,rsig2,rpi2,xi,xi2,exi,xup,yup,
    xm,xm2,xq,xq2,em,eq,r;

  /*
  power=pow(2.0,abs(ipow));
  if(ipow<0)
    sprintf(mstar,"mstar%d",power);
  if(ipow>=0)
    sprintf(mstar,"%dmstar",power);
  */

  /* New file notation:
   */
  sprintf(mstar,"g%d",ipow);


  for(i=0;i<5;++i)
    for(j=0;j<3;++j)
      {
	sprintf(filename,"run%d/xi_ang.%02d.%s.%d",i+1,fnumber,mstar,j+1);
	if(!(fp[i][j]=fopen(filename,"r")))
	  {
	    fprintf(stderr,"ERROR with file %d %d [%s]\n",i+1,j+1,filename);
	    return;
	  }
      }
  
  sprintf(filename,"xi_ang.%02d.%s.avg",fnumber,mstar);
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

  rewind(fp[0][0]);    

  for(i=1;i<=lines;++i)
    {
      r=0;
      xm=xm2=0;
      xq=xq2=0;
      nr=0;
      for(j=0;j<5;++j)
	for(k=0;k<3;++k)
	  {
	    fscanf(fp[j][k],"%f %f %f",&x1,&x2,&x3);
	    nr++;
	    r+=x1;
	    xm+=x2;
	    xm2+=x2*x2;
	    xq+=x3;
	    xq2+=x3*x3;
	  }
      r/=nr;
      xm/=nr;
      xm2/=nr;
      xq/=nr;
      xq2/=nr;

      em=sqrt((xm2-xm*xm)/(nr-1.0));
      eq=sqrt((xq2-xq*xq)/(nr-1.0));
      
      fprintf(outf,"%f %e %e %e %e\n",r,xm,em,xq,eq);
    }
  fclose(outf);
  for(i=0;i<5;++i)
    for(j=0;j<3;++j)
      fclose(fp[i][j]);
}
