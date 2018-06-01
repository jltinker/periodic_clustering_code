#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void average_files(int ipow, int ipow2, int fnumber);

int main(int argc, char **argv)
{
  int i,fn,plo,phi,j,k;
  
  if(argc<4)
    {
      fprintf(stderr,"avg-2dx-corr pow2low pow2hi no_of_outputs\n\n");
      exit(0);
    }

  plo=atoi(argv[1]);
  phi=atoi(argv[2]);
  fn=atoi(argv[3]);

  for(i=1;i<=fn;++i)
    for(j=plo;j<=phi;++j)
      for(k=j+1;k<=phi;++k)
	average_files(j,k,i);

  return(0);
}

void average_files(int ipow, int ipow2, int fnumber)
{
  int n,j,i,ix,lines=0,k,nr,power,nsig,npi,nxi;
  FILE *fp[5][3],*outf;
  char filename[100],mstar[10],mstar2[10];
  float x1,x2,x3,rsig,rpi,esig,epi,rsig2,rpi2,xi,xi2,exi,xup,yup;

  power=pow(2.0,abs(ipow));
  if(ipow<0)
    sprintf(mstar,"mstar%d",power);
  if(ipow>=0)
    sprintf(mstar,"%dmstar",power);

  power=pow(2.0,abs(ipow2));
  if(ipow2<0)
    sprintf(mstar2,"mstar%d",power);
  if(ipow2>=0)
    sprintf(mstar2,"%dmstar",power);

  for(i=0;i<5;++i)
    for(j=0;j<3;++j)
      {
	sprintf(filename,"run%d/xi2d.%02d.%s.%s.%d",i+1,fnumber,mstar,mstar2,j+1);
	if(!(fp[i][j]=fopen(filename,"r")))
	  {
	    fprintf(stderr,"ERROR with file %d %d [%s]\n",i+1,j+1,filename);
	    return;
	  }
      }
  
  sprintf(filename,"xi2d.%02d.%s.%s.avg",fnumber,mstar,mstar2);
  if(!(outf=fopen(filename,"w")))
    {
      fprintf(stderr,"ERROR opening outfile: [%s]\n",filename);
      kill("Exiting...\n");
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
      for(j=0;j<5;++j)
	for(k=0;k<3;++k)
	  {
	    fscanf(fp[j][k],"%f %f %f %f %f",&x1,&x2,&x3,&xup,&yup);
	    if(x1>-0.99)
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
      exi=sqrt(xi2-xi*xi);
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
  for(i=0;i<5;++i)
    for(j=0;j<3;++j)
      fclose(fp[i][j]);
}
