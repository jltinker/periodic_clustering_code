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

  int ibin, nbin=6, bincnt[10];
  double avgmass[10], avgsat[10], avgnum[10];
  float M1, M_cut, sigma_logM, alpha,  ncen, nsat; 

  fp = openfile(argv[1]);

  for(i=1;i<=nbin;++i)
    avgmass[i] = bincnt[i] = avgsat[i] = avgnum[i] = 0;

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
  
  // read in the HOD from a file
  fp = fopen(argv[4],"r");
  fscanf(fp,"%f %f %f %f %f",&x1,&M1, &alpha, &M_cut, &sigma_logM);
  
  
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

      // if above 10, then get the M/N for this
      if(Lobs<10)
	continue;

      // get the number of satellites
      //ncen = 0.5*(1+erf((log10(mass) - log10(M_min))/sigma_logM));
      ncen = 1;
      nsat = exp(-M_cut/mass)*pow(mass/M1,alpha)*ncen;
      //printf("%e %e %e\n",mass, nsat, Lobs);
	   

      // put this in the proper bin
      if(Lobs>=10 && Lobs<11)ibin = 1;
      if(Lobs>=11 && Lobs<13)ibin = 2;
      if(Lobs>=13 && Lobs<17)ibin = 3;
      if(Lobs>=17 && Lobs<25)ibin = 4;
      if(Lobs>=25 && Lobs<41)ibin = 5;
      if(Lobs>=41)ibin = 6;
      
      avgmass[ibin] += mass;
      avgnum[ibin] += (ncen+nsat);
      avgsat[ibin] += nsat;
      bincnt[ibin]++;
    }


  // rho_crit div nbar for this sample
  fac = 2.775e11/0.0003;
  for(i=1;i<=nbin;++i)
    printf("%e %e %e %e %e %d\n",avgmass[i]/bincnt[i],avgmass[i]/avgnum[i]/fac, avgmass[i]/avgsat[i]/fac,
	   avgnum[i]/bincnt[i], avgsat[i]/bincnt[i], bincnt[i]);

}
