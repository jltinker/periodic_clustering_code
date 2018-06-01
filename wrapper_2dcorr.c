#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

void help_message(void);

int main(int argc, char **argv)
{
  int power,i,j,k,nstar,n1,n2,nlo,nhi,n,nbin,num,outlo,outhi;
  FILE *fp;
  char command[200],outfile[100],xifile[100],anum[2];
  float rlo,rhi,rcube;

  if(argc<7)
    help_message();

  outlo=1;
  outhi=5;
  if(argc>=9)
    {
      outlo=atoi(argv[7]);
      outhi=atoi(argv[8]);
    }

  n1=atoi(argv[1]);
  n2=atoi(argv[2]);

  rlo=atof(argv[3]);
  rhi=atof(argv[4]);
  nbin=atoi(argv[5]);
  rcube=atof(argv[6]);

  for(j=outlo;j<=outhi;++j)
    {
      for(i=n1;i<=n2;++i)
	{
	  power=pow(2.0,abs(i));
	  n=nstar*pow(2.0,i);
	  nlo=n/sqrt(2);
	  nhi=n*sqrt(2);
	  if(i<0)
	    sprintf(outfile,"halo.%02d.mstar%d",j,power);
	  if(i>=0)
	    sprintf(outfile,"halo.%02d.%dmstar",j,power);
	  if(i<0)
	    sprintf(xifile,"xi2d.%02d.mstar%d",j,power);
	  if(i>=0)
	    sprintf(xifile,"xi2d.%02d.%dmstar",j,power);
	  for(k=1;k<=3;++k)
	    {
	      sprintf(command,"2dcorr %.2f %.2f %d %.1f 0 %.1f 1 %s a 0 1 %d %d > %s.%d",
		      rlo,rhi,nbin,rcube,rcube,outfile,j,k,xifile,k);
	      fprintf(stderr,">%s\n",command);
	      system(command);
	    }
	  /*
	  sprintf(command,"awk -f avg.awk out.1 > %s",xifile);
	  fprintf(stderr,">%s\n",command);
	  system(command);
	  */
	}  
    }

}

void help_message()
{
  fprintf(stderr,"\nwrapper_2dcorr pow2low pow2hi rmin rmax nbin rcube\n\n");
  fprintf(stderr," > This wrapper calls 2dcorr to get halos in mass bins:\n");
  fprintf(stderr," > N_star*pow(2,lo) -> N_star*pow(2,hi)\n\n");
  fprintf(stderr," > Then the wrapper calls 2dcorr assuming autocorrelation, no\n");
  fprintf(stderr," > and the parameters above.\n\n");
  kill("Exiting...");
}
