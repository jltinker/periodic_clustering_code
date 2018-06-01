#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void average_files(int ipow, int fnumber);

int NRUNS=5;
char PREFIX[100]="xi2d";

int main(int argc, char **argv)
{
  int i,fn,plo,phi,j,outlo,outhi;
  
  if(argc<4)
    {
      fprintf(stderr,"avg-files outlo outhi pow2low pow2hi [NRUNS] [PREFIX]\n\n");
      exit(0);
    }

  outlo=atoi(argv[1]);
  outhi=atoi(argv[2]);
  plo=atoi(argv[3]);
  phi=atoi(argv[4]);
  if(argc>5)
    NRUNS=atoi(argv[5]);
  if(argc>6)
    sprintf(PREFIX,"%s",argv[6]);

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
  float x1,x2,x3,rsig,rpi,esig,epi,rsig2,rpi2,xi,xi2,exi,xup,yup;

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
	if(NRUNS==1)
	  sprintf(filename,"%s.%d.%s.%d",PREFIX,fnumber,mstar,j+1);
	else
	  sprintf(filename,"run%d/%s.%02d.%s.%d",i+1,PREFIX,fnumber,mstar,j+1);

	if(!(fp[i][j]=fopen(filename,"r")))
	  {
	    fprintf(stderr,"ERROR with file %d %d [%s]\n",i+1,j+1,filename);
	    return;
	  }
      }
  
  sprintf(filename,"%s.%02d.%s.avg",PREFIX,fnumber,mstar);
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
      rsig=rpi=xi=0;
      esig=epi=exi=0;
      rsig2=rpi2=xi2=0; 
      nr=0;
      for(j=0;j<NRUNS;++j)
	for(k=0;k<3;++k)
	  {
	    fscanf(fp[j][k],"%f %f %f %f %f",&x1,&x2,&x3,&xup,&yup);
	    if(x1>-1.99)
	      {
		nr++;
		xi+=x1;
		xi2+=x1*x1;
		rsig+=x2;
		rsig2+=x2*x2;
		rpi+=x3;
		rpi2+=x3*x3;
	      }
	  }
      xi/=nr; xi2/=nr;
      rsig/=nr;
      rpi/=nr;
      exi=sqrt(xi2-xi*xi)/sqrt(nr-1);
      if(!nr)
	{
	  xi=-1;
	  rsig=x2;
	  rpi=x3;
	}
      
      /* I think for now all I will do is put xi, rsig, rpi
       */
      /*
      if(nr>1)er=sqrt((r2-r*r)/(nr-1.0));
      */
      fprintf(outf,"%e %f %f %e %f %f\n",xi,rsig,rpi,exi,xup,yup);
    }
  fclose(outf);
  for(i=0;i<NRUNS;++i)
    for(j=0;j<3;++j)
      fclose(fp[i][j]);
}
