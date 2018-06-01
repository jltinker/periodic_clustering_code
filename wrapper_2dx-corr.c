#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

void help_message(void);

int main(int argc, char **argv)
{
  int power,i,j,k,kk,nstar,n1,n2,nlo,nhi,n,nbin,num,outlo,outhi,pow2lo,pow2hi;
  FILE *fp;
  char command[400],outfile[100],xifile[100],anum[2],file1[100],file2[100],mstar1[10],mstar2[10];
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

  pow2lo=atoi(argv[1]);
  pow2hi=atoi(argv[2]);

  rlo=atof(argv[3]);
  rhi=atof(argv[4]);
  nbin=atoi(argv[5]);
  rcube=atof(argv[6]);

  for(i=outlo;i<=outhi;++i)
    for(j=pow2lo;j<=pow2hi;++j)
      for(k=j+1;k<=pow2hi;++k)
	{
	  power=pow(2.0,abs(j));
	  if(j<0)
	    sprintf(mstar1,"mstar%d",power);
	  if(j>=0)
	    sprintf(mstar1,"%dmstar",power);
	  sprintf(file1,"halo.%02d.%s",i,mstar1);

	  power=pow(2.0,abs(k));
	  if(k<0)
	    sprintf(mstar2,"mstar%d",power);
	  if(k>=0)
	    sprintf(mstar2,"%dmstar",power);
	  sprintf(file2,"halo.%02d.%s",i,mstar2);
	  
	  sprintf(xifile,"xi2d.%02d.%s.%s",i,mstar1,mstar2);

	  for(kk=1;kk<=3;++kk)
	    {
	      sprintf(command,
		      "2dx-corr %.2f %.2f %d %.1f 0 %.1f 1 %s a 0 1 %s a 0 1 %d %d > %s.%d",
		      rlo,rhi,nbin,rcube,rcube,file1,file2,i,kk,xifile,kk);
	      fprintf(stderr,">%s\n",command);
	      system(command);
	    }
	}  
}

void help_message()
{
  fprintf(stderr,"\nwrapper_2dcorr pow2low pow2hi rmin rmax nbin rcube [outlo] [outhi]\n\n");
  fprintf(stderr," > This wrapper calls 2dcorr to get halos in mass bins:\n");
  fprintf(stderr," > N_star*pow(2,lo) -> N_star*pow(2,hi)\n\n");
  fprintf(stderr," > Then the wrapper calls 2dcorr assuming autocorrelation, no\n");
  fprintf(stderr," > and the parameters above.\n\n");
  kill("Exiting...");
}
