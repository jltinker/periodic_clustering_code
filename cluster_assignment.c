#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

/*
 * implementing equation B4 from Rykoff et al 2012 redmapper paper
 *  ln(M200b/10^14*h70) = 1.72 + 1.08*ln(L/60)
 * 
 */
float gasdev(long *idum);


int main(int argc, char **argv)
{
  int i,j,k,nhalo, i1, imass;
  FILE *fp;
  float mass, r1, vc, x[6];
  char aa[1000];
  long IDUM=-555;

  float *hmass, *xxh, *yxh, *zxh, *vxh, *vyh, *vzh;
  float lnL60, L, sig = 0.22/1.08, Lobs, fac, pmass, x1, omega, resolution;

  fp = openfile(argv[1]);

  // check to see if this is an FOF halo file or a QPM halo file
  if(argc>2) 
    {
      nhalo = filesize(fp);
      hmass = vector(1,nhalo);
      xxh = vector(1,nhalo);
      yxh = vector(1,nhalo);
      zxh = vector(1,nhalo);
      vxh = vector(1,nhalo);
      vyh = vector(1,nhalo);
      vzh = vector(1,nhalo);

      omega = atof(argv[2]);
      resolution = atof(argv[3]);
      pmass = 2.775e11*omega/(resolution*resolution*resolution);
      fprintf(stderr,"Particle mass: %e\n",pmass);
      for(i=1;i<=nhalo;++i)
	{
	  fscanf(fp,"%d %d %f %f %f %f %f %f %f",&i1,&imass,&xxh[i],&yxh[i],&zxh[i],&x1,&vxh[i],&vyh[i],&vzh[i]);
	  hmass[i] = imass*pmass;
	  fgets(aa,1000,fp);
	}
      fclose(fp);
    }
  else 
    {
      fread(&nhalo,sizeof(int),1,fp);
      
      // get total number of halos and allocate memory
      hmass = vector(1,nhalo);
      xxh = vector(1,nhalo);
      yxh = vector(1,nhalo);
      zxh = vector(1,nhalo);
      vxh = vector(1,nhalo);
      vyh = vector(1,nhalo);
      vzh = vector(1,nhalo);
      
      
      // for the QPM halo files.
      fprintf(stderr,"Reading in the PB-halos [%d]...\n",nhalo);
      fread(&hmass[1],sizeof(float),nhalo,fp);
      fread(&xxh[1],sizeof(float),nhalo,fp);
      fread(&yxh[1],sizeof(float),nhalo,fp);
      fread(&zxh[1],sizeof(float),nhalo,fp);
      fread(&vxh[1],sizeof(float),nhalo,fp);
      fread(&vyh[1],sizeof(float),nhalo,fp);
      fread(&vzh[1],sizeof(float),nhalo,fp);
      fclose(fp);
    }

  //sig = 0.25;
  for(i=1;i<=nhalo;++i)
    {
      //fscanf("%d %f %f %f",&j,&mass, &r1, &vc);
      //for(j=0;j<6;++j)fscanf("%f",&x[j]);
      //fgets(aa,1000,fp);

      mass = hmass[i];

      // get function for sig^2 factor (approximates dln(n)/dln(m) of HMF)
      fac = 1 + mass/0.8E14;

      // get mean lam for this mass from (B4)
      lnL60 = (log(mass/1.4E14) - 1.72)/1.08;
      L = exp(lnL60 - fac*sig*sig/2)*60;

      // gaussian deviate with scatter 0.25*1.08
      Lobs = exp(gasdev(&IDUM)*sig)*L ;

      // print out the result:
      if(Lobs>10)
	printf("%e %e %e 0.0 0.0 0.0 %e %f %f\n",xxh[i], yxh[i], zxh[i], mass, Lobs, L);
      continue;

      if(mass>1.0E14)
	printf("%e %e %e 0.0 0.0 0.0 %e %f %f\n",xxh[i], yxh[i], zxh[i], mass, Lobs, L);
      continue;

    }


}
