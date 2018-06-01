#include <stdio.h>
#include <math.h>
#include "nrutil.h"

#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define pi 3.1415926535898

#define muh(x) fprintf(stderr,"muh %d\n",x)
#define muf(x) fprintf(stderr,"muh %f\n",x)

#define HANDHOLD 0
#define NBINLOOKUP 50000     		/* entries in radial bin lookup table */
#define NBINMAX1 101			/* max. number of xi bins, +1 */
#define RBINMAX1 101

void die(char *errmsg) ;
void getposvel_fastfood(FILE *fp, float zstep, 
			float **x, float **y, float **z,
			float **vx, float **vy, float **vz, int *np) ;
void getposvel_tipsy(FILE *fp, double time, char type, 
                 float **x, float **y, float **z,
		 float **vx, float **vy, float **vz,
                 int *nsph, int *ndark, int *nstar, int *np)  ;
void getposvel_ascii(FILE *fp, 
		     float **x, float **y, float **z,
		     float **vx, float **vy, float **vz, int *np) ;

/** TINKER **/
void brute_force_method(int p, float *x, float *y, float *z, int npair[101][101], int binlookup[], 
			float binfac, int np, int nrbin, float zmax);
void redshift_distortions(float *z, float *vz, float rcube, int np, int imode);

main(int argc,char **argv)
{
  int np1,np2,          /* particle numbers in 2 samples */
    nlattice,		/* chaining mesh lattice size */
    i, j, k,		/* assorted integers */
    nrbin,		/* number of radial bins */
    nbin,               /* number of projected seperation bins */
    ibin,kbin,rbin,		/* bin counters */
    nsph,ndark,nstar,	/* for use in tipsy binary reads */
    iauto,		/* 1->autocorrelation, 0->cross-correlation */
    ioct,ioct0,		/* octant counter */
    nreset,		/* number of particles with reset coordinates */
    nrmax ;		/* number of particles within rmax */
  
  double weight0,	/* sum of weight over octants */
    fac;		/* 1/(summed weight in bin) */

  float rmin,rmax,rmax2,	/* radial range, and square of rmax */
    zmin,zmax,          /* range in line-of-site distance */
    rcube,		/* size of cube in Mpc/h */
    rhalf,		/* half of rcube */
    xmin,xmax,		/* range of particle coordinates */
    vscale,		/* velocity scaling factor */
    time1,time2,        /* times to read from input files */
    rscale,		/* scaling factor for positions */
    avgweight1,		/* mean weight of particles in set 1 */
    avgweight2,		/* mean weight of particles in set 2 */
    cmin,cmax,		/* used to reset particles on boundary */
    lrstep,		/* logarithmic radial step */
    binfac,		/* used in radial bin calculation */
    xp,yp,zp,		/* position of reference particle */
    vxp,vyp,vzp,	/* velocity of reference particle */
    dx,dy,dz,		/* difference in position */
    dvx,dvy,dvz,	/* difference in velocity */
    r,   		/* radial separation */
    v12,		/* line-of-sight velocity */
    rlow,		/* lower radius of bin */
    density,		/* mean weighted particle density in cube */
    vstream0,		/* sum of vstream over octants */
    vpar0,		/* sum of vpar2 over octants */
    vtan0,		/* sum of vtan2 over octants */
    vxi0,		/* sum of vxi2 over octants */
    vol,		/* volume of radial bin */
    weightrandom,	/* weight expected for random distribution */
    xi,			/* correlation function */
    ovstream,		/* output value of streaming velocity */
    ovpar,		/* output value of parallel dispersion */
    ovtan,		/* output value of tangential dispersion */
    ovxi,		/* output value of velocity correlation */
    errxi,		/* error in correlation function */
    errvstream,		/* error in streaming velocity */
    errvpar,		/* error in parallel dispersion */
    errvtan,		/* error in tangential dispersion */
    errvxi,		/* error in velocity correlation */
    jackweight,		/* weight in jackknife subsample */
    diff,		/* difference of jackknife value from mean */
    vstream1,		/* streaming velocity in jackknife sample */
    vxi1 ,		/* velocity correlation in jackknife sample */
    v;                  /* velocity magnitude of a particle */
  
  double weight ;	/* product of weights of two particles */
  
  int *listbods, **lattice, 	/* linklist arrays (2-D)*/
    *indx,			/* index array for particles within rmax */
    *binlookup ;		/* radial bin given separation */
  
  float *rupp,		/* upper radius of bin */
    *rsqr,		/* squared distance for particles within rmax */
    *x1, *y1, *z1, 	/* positions for first particle set */
    *vx1, *vy1, *vz1,   /* velocities for first particle set */
    *weight1,		/* weights for first particle set */
    *x2, *y2, *z2, 	/* positions for second particle set */
    *vx2, *vy2, *vz2,	/* velocities for second particle set */
    *weight2 ;		/* weights for second particle set */
  
  /* note that second subscript in subsequent arrays is for values in 
   * individual octants 
   */
  int npair[NBINMAX1],          /* number of pairs in radial bin */
    npair2d[NBINMAX1][RBINMAX1];/* number of pairs in 2d grid */
  double rbar[NBINMAX1],	/* mean separation in radial bin */
    weightsum[NBINMAX1][8],     /* summed weight product */
    vstream[NBINMAX1][8],	/* mean streaming velocity */
    vpar2[NBINMAX1][8],		/* parallel pairwise dispersion */
    vtan2[NBINMAX1][8],		/* tangential pairwise dispersion */
    vxi2[NBINMAX1][8] ;		/* velocity correlation function */

  double
    weightsum2d[NBINMAX1][RBINMAX1][8] ; /* number of pairs in [sigma,pi] bin */

  char *calloc() ;
  FILE *fp ;

  /** TINKER's workspace **/
  float rpi,rsigma,             /* coordinates in (pi,sigma) */
    density2d,                  /* surface density */
    density1d,                  /* linear density */
    los_low,                    /* lower bound of l-o-s bin */
    phi,                        /* angle between sigma,pi */
    *temp_arr;
 
  double  sigbar[NBINMAX1][RBINMAX1], /* average value of projected seperation */
    pibar[NBINMAX1][RBINMAX1],  /* average value of radial seperation */
    mbar[NBINMAX1][NBINMAX1];
  int iredshift=0;              /* flag for distorting z-positions */
  float DELTA_PHI,phi_low,xi_array[NBINMAX1][RBINMAX1];
  int redshift_flag,            /* Output file number */
    velocity_axis;              /* which direction is line-f-sight */

  FILE *jackfile;
  char jackname[100];

  // Stuff for MR mass stuff
  float halomass,**xtemp;
  int ngals,ihalo,itemp,**ixtemp,total_sats;
  float x[18];
  int ix[3];

  if (argc < 13 )  {
    fprintf(stderr,"usage: 2dcorr rmin rmax nbin rcube xmin xmax vscale ") ;
    fprintf(stderr,
	    "file1 format1 time1 weight1 [output #] [los axis->1,2,3]\n");
    die(" ") ;
  }


  if(argc>=15)
    sprintf(jackname,"%s",argv[14]);
  else
    sprintf(jackname,"xi2d.jackfile");

  /* read control parameters */
  sscanf(argv[1],"%f",&rmin) ;		
  sscanf(argv[2],"%f",&rmax) ;
  sscanf(argv[3],"%d",&nbin) ;
  sscanf(argv[4],"%f",&rcube) ;
  sscanf(argv[5],"%f",&xmin) ;
  sscanf(argv[6],"%f",&xmax) ;
  sscanf(argv[7],"%f",&vscale) ;
  if (nbin>NBINMAX1-1)
    die("requested number of bins exceeds maximum allowed") ;


  /* set the line-of-sight range equal to the input range */
  zmin=rmin;
  zmax=rmax;
  nrbin=nbin;
  
  avgweight1 = avgweight2 = 1;

  /* read or assign weights for first object list */
  weight1 = (float *) calloc(np1,sizeof(float)) ;	
  if (strcmp(argv[11],"1")!=0)  {
    fp=fopen(argv[11],"r") ;
    if (fp==NULL)
      die("couldn't open first weight file") ;
    for (i=0;i<np1;i++)  {
      if (fscanf(fp,"%f",&weight1[i])==EOF)  {
	die("premature end of first weight file") ;
      }
    }
    fclose(fp) ;
  }
  else  {
    for (i=0;i<np1;i++)  {
      weight1[i]=1.0 ;
    }
  }


  /*
   * currently only doing the auto-correlation function.
   */

  iauto=1 ;
  x2=x1 ;
  y2=y1 ;
  z2=z1 ;
  vx2=vx1 ;
  vy2=vy1 ;
  vz2=vz1 ;
  weight2=weight1 ;
  np2=np1 ;

  /* rescale positions (to 0--1) and velocities, culling weight=0 objects along the 
   * way for efficiency 
   */
  rhalf=rcube/2.0 ;

  /*********************** 
   *initializing the logarithmic bins 
   */  
  rupp = (float *) calloc(nbin+1,sizeof(float)) ;
  lrstep=log(rmax/rmin)/(float)(nbin-1) ;
  binlookup=(int *)calloc(NBINLOOKUP+2,sizeof(int)) ;
  ibin=0 ;
  for (i=0;i<=NBINLOOKUP;i++)  {
    r=rmax*i/NBINLOOKUP ;
    if (r>0)  {
      kbin=(int)floor(log(r/rmin)/lrstep+1.0) ;

      /* if bins are linear use this 
       *
       * kbin=(int)(r/rmax*nbin);
       */
    }
    else {
      kbin=0 ;
    }
    if (kbin<0) kbin=0 ;
    if (kbin>ibin)  {
      rupp[ibin]=r ;
      ibin=kbin ;
    }
    binlookup[i]=kbin ;
  }
  binlookup[NBINLOOKUP+1]=nbin ;
  rupp[nbin-1]=rmax ;
  rupp[nbin]=rmax ;
  binfac=NBINLOOKUP/rmax ;


  /* ensure that sums start from zero */
  for (i=0;i<NBINMAX1;i++)  {
    npair[i]=0 ;
    rbar[i]=0.0 ;
    for(j=0;j<RBINMAX1;j++)
      {
	sigbar[i][j]=0;
	pibar[i][j]=0;
	npair2d[i][j]=0;
	mbar[i][j]=0;
      }
    for (j=0;j<8;j++)  {
      weightsum[i][j]=0.0 ;
      vstream[i][j]=0.0 ;
      vpar2[i][j]=0.0 ;
      vtan2[i][j]=0.0 ;
      vxi2[i][j]=0.0 ;
    }
  }

  /* final preparations */
  indx = (int *) calloc(np2,sizeof(int)) ;	
  rsqr = (float *) calloc(np2,sizeof(float)) ;	
  if (indx==NULL || rsqr==NULL)  {
    die("not enough memory for indx and rsqr arrays") ;
  }
  rmax2=rmax*rmax ;
  

  /***********************************************
   * Do the actual correlation function.
   */

  xp=0;

  fp = fopen("allfixed","r");
  //  fread(&ngals,1,sizeof(int),fp);
  ngals = 794634;

  
  xtemp = matrix(0,ngals-1,0,17);
  ixtemp = imatrix(0,ngals-1,0,2);

  /*
  for(j=0;j<ngals;++j)
    {
      fread(ixtemp[j],3,sizeof(int),fp);
      fread(xtemp[j],18,sizeof(float),fp);
    }
  */
  for(j=0;j<ngals;++j)
    {
      for(i=0;i<3;++i)
	fscanf(fp,"%d",&ixtemp[j][i]);
      for(i=0;i<18;++i)
	fscanf(fp,"%f",&xtemp[j][i]);
      xtemp[j][10] = -22;
    }
  
  fprintf(stderr,"Done reading %d galaxies.\n",ngals);
  
  /* Read in the galaxies for each halo.
   */
  //  fread(ix,3,sizeof(int),fp);
  //fread(x,18,sizeof(float),fp);

  
  for(i=0;i<3;++i)
    ix[i] = ixtemp[0][i];
  for(i=0;i<18;++i)
    x[i] = xtemp[0][i];
  

  /* Create array for gals.
   */
  x1 = malloc(100000*sizeof(float));
  y1 = malloc(100000*sizeof(float));
  z1 = malloc(100000*sizeof(float));

  ihalo = 0;
  itemp = 0;
  total_sats=0;
  for(j=1;j<ngals;++j)
    {
      if(x[10]>-19.5)goto SKIP_GAL;
      
      i = itemp;
      /*
      x1[i] = x[0] + x[6]/100;
      if(x1[i]>rcube)x1[i]-=rcube;
      if(x1[i]<0)x1[i]+=rcube;

      y1[i] = x[1] + x[7]/100;
      if(y1[i]>rcube)y1[i]-=rcube;
      if(y1[i]<0)y1[i]+=rcube;
      */
      x1[i] = x[0];
      y1[i] = x[1];
      z1[i] = x[2] + x[8]/100;
      if(z1[i]>rcube)z1[i]-=rcube;
      if(z1[i]<0)z1[i]+=rcube;
      i++;itemp++;

      if(ix[2]==0 || ix[2]==4)
	halomass = x[15];

    SKIP_GAL:
      // fread(ix,3,sizeof(int),fp);
      //fread(x,18,sizeof(float),fp);
      
      for(i=0;i<3;++i)
	ix[i] = ixtemp[j][i];
      for(i=0;i<18;++i)
	x[i] = xtemp[j][i];
      

      if(ix[2]!=0 && ix[2]!=4)continue;
      ihalo++;

      np1 = itemp;
      itemp = 0;

      //      if(ihalo==100)exit(0);
      if(np1<2)continue;
      if(np1>100 && halomass<14) {
	printf("HALO %d %f %d\n",ihalo,halomass,np1);fflush(stdout);
	continue;
      }
      total_sats += np1-1;


      /* loop over particles in 1st list */
      for (i=0;i<np1;i++)  {		    
	xp=x1[i] ;
	yp=y1[i] ;
	zp=z1[i] ;
	
	/* compute octant index based on particle position */
	ioct0=0 ;
	if (xp>rhalf)  ioct0++ ;
	if (yp>rhalf)  ioct0+=2 ;
	if (zp>rhalf)  ioct0+=4 ;
	
	for (k=0;k<np1;k++)  {
	  if(k==i)continue;

	  dx = mabs(x1[i] - x1[k]);
	  if(dx>rhalf)dx=rcube-dx;
	  dy = mabs(y1[i] - y1[k]);
	  if(dy>rhalf)dy=rcube-dy;
	  dz = mabs(z1[i] - z1[k]);
	  if(dz>rhalf)dz=rcube-dz;
	  rsqr[k] = dx*dx + dy*dy;

	  r=sqrt(rsqr[k]) ;
	  if(r>rmax)continue;

	  kbin=binlookup[(int)(binfac*r)] ;
	  dz=mabs(z1[i]-z1[k]);
	  if(dz>rhalf)dz=rcube-dz;
	  if(dz<=zmax)
	    rbin=binlookup[(int)(binfac*dz)];
	  else
	    rbin=nrbin+1;

	  /* for an autocorrelation, choose the octant based on particle i half the
	   * time (when j is odd) and particle j half the time, to avoid biases
	   * in number of pairs per octant caused by only looking at pairs with
	   * j > i 
	   */
	  if (iauto==1 && j%2==0)  {
	    ioct=0 ;
	    if (x1[k]>rhalf)  ioct++ ;
	    if (y1[k]>rhalf)  ioct+=2 ;
	    if (z1[k]>rhalf)  ioct+=4 ;
	  }
	  else  {
	    ioct=ioct0 ;
	  }

	  /* add weight to total, which will be used for xi(r) and weighting vel stats */
	  weight=1;//weight1[i]*weight2[j] ;
	  weightsum[kbin][ioct]+=weight ;
	  rbar[kbin]+=weight*r ;
	  npair[kbin]++ ;
	  if(rbin<nrbin)
	    {
	      weightsum2d[kbin][rbin][ioct]+=weight;
	      npair2d[kbin][rbin]++;
	      sigbar[kbin][rbin]+=weight*r;
	      pibar[kbin][rbin]+=weight*dz;
	      mbar[kbin][rbin]+=halomass;
	    }
	}
      }
    }

  fprintf(stderr,"sats: %d\n",total_sats);

  /* compute the output quantities and error bars */
  rlow=0.0 ;
  density=np2*avgweight2/(rcube*rcube*rcube) ;
  density2d=density*rcube;      /* surface density of objects */
  density1d=density2d*rcube;

  if (iauto==1)  {		/* pairs not double counted */
    density=density/2. ;
    density2d=density2d/2.;
    density1d=density1d/2;
  }

  for (kbin=0;kbin<nbin;kbin++)  {	/* loop over radial bins */

    weight0=0.0 ;

    for (ioct=0;ioct<8;ioct++)  {		/* add up values from octants */
      weight0+=weightsum[kbin][ioct] ;
    }
    
    if (weight0>0.0)  {
      fac=1./weight0 ;
      rbar[kbin] *= fac ;
    }
    else  {					/* avoid errors in empty bins */
      fac=0.0 ;
      rbar[kbin]=(rupp[kbin]+rlow)/2. ;
    }
    
    /* compute xi, dividing summed weight by that expected for a random set */
    vol=pi*(rupp[kbin]*rupp[kbin]-rlow*rlow) ;
    weightrandom=np1*avgweight1*density2d*vol ;	
    xi=weight0/weightrandom-1 ;

    /* TINKER temp-- volume of wedge of spherical shell 
     *
     vol=4./3.*pi*(rupp[kbin]*rupp[kbin]*rupp[kbin]-rlow*rlow*rlow);
     vol*=(-cos(0.5*pi-phi_low*pi/180.)+cos(0.5*pi-(phi_low+DELTA_PHI)*pi/180.))*2.0;
     weightrandom=np1*avgweight1*density*vol ;	
     xi=weight0/weightrandom-1 ;
    */

    /* compute the jackknife errors by subtracting each octant from the total
     * in turn, computing the difference between this jackknife subsample (seven
     * of the eight octants) and the result for the whole cube, and summing
     * the differences in quadrature.
     */
    weightrandom *= 7.0/8.0 ;	/* each jack sample is 7/8 of volume */
    errxi=0.0 ;
    for (ioct=0;ioct<8;ioct++)  {
      jackweight=weight0-weightsum[kbin][ioct] ;
      if (jackweight>0.0)  {
	fac=1./jackweight ;
      }
      else  {
	fac=0.0 ;
      }
      diff=jackweight/weightrandom-1-xi ;	  /* sample xi minus mean */
      errxi+=diff*diff ;
    }
    
    errxi=sqrt(errxi*8./7.) ;

    /* This is the calculation of the 1-d correlation function.
     * We will not print this out.
     printf("%7.4f %7.4f %7.4f %10.3e %10.3e %10d %10.0f\n",
     rlow,rupp[kbin],rbar[kbin],xi,errxi,npair[kbin],weightrandom*8.0/7.0) ;
     */
    rlow=rupp[kbin] ;
  }					/* next radial bin */


  /******************************************** 
   * now do the same thing for the 2d-array 
   */
  
  rlow=0;

  jackfile = fopen(jackname,"w");

  for (kbin=0;kbin<nbin;kbin++)  {	/* loop over proj-sep bins */

    /* loop over the l-o-s bins */

    los_low=0;
    for (rbin=0;rbin<nrbin;rbin++)
      {
	weight0=0.0 ;

	for (ioct=0;ioct<8;ioct++)  {		/* add up values from octants */
	  weight0+=weightsum2d[kbin][rbin][ioct] ;
	}
	if (weight0>0.0)  {
	  fac=1./weight0 ;
	  sigbar[kbin][rbin] *= fac ;
	  pibar[kbin][rbin] *= fac;
	  mbar[kbin][rbin] *= fac;
	}
	else  {					/* avoid errors in empty bins */
	  fac=0.0 ;
	  sigbar[kbin][rbin]=(rupp[kbin]+rlow)*0.5;
	  pibar[kbin][rbin]=(rupp[rbin]+los_low)*0.5;
	}
    
	/* compute xi, dividing summed weight by that expected for a random set */
	vol=pi*(rupp[kbin]*rupp[kbin]-rlow*rlow)*(rupp[rbin]-los_low)*2.0;
	weightrandom=np1*avgweight1*density*vol ;	
	xi=weight0/weightrandom-1 ;

	/* compute the jackknife errors by subtracting each octant from the total
	 * in turn, computing the difference between this jackknife subsample (seven
	 * of the eight octants) and the result for the whole cube, and summing
	 * the differences in quadrature.
	 */
	weightrandom *= 7.0/8.0 ;	/* each jack sample is 7/8 of volume */
	errxi=0.0 ;
	for (ioct=0;ioct<8;ioct++)  {
	  jackweight=weight0-weightsum2d[kbin][rbin][ioct] ;
	  if (jackweight>0.0)  {
	    fac=1./jackweight ;
	  }
	  else  {
	    fac=0.0 ;
	  }
	  fprintf(jackfile,"%e ",jackweight/weightrandom-1);
	  diff=jackweight/weightrandom-1-xi ;	  /* sample xi minus mean */
	  errxi+=diff*diff ;
	}
	fprintf(jackfile,"\n");

	/* pi = radial sep
	 * sigma = projected sep
	 */
	rpi=(rupp[rbin]+los_low)*0.5;
	rsigma=(rupp[kbin]+rlow)*0.5;

	errxi=sqrt(errxi*7./8.) ;
	xi_array[kbin][rbin]=xi;
	
	printf("%f %f %e %e %d %.0f %e\n",sigbar[kbin][rbin],pibar[kbin][rbin],xi,errxi,npair2d[kbin][rbin],weightrandom*8./7.,mbar[kbin][rbin]);
	/*
	printf("%f %f %f %f %f\n",xi,sigbar[kbin][rbin],pibar[kbin][rbin],rupp[kbin],rupp[rbin]);
	printf("%f\n",xi);
	printf("%d %d %f %.0f\n",kbin,rbin,xi,weight0);
	printf("%7.4f %7.4f %10.3e %10.3e %10d %.0f %f\n",
	       rsigma,rpi,xi,errxi,npair2d[kbin][rbin],weightrandom*8./7.,
	       npair2d[kbin][rbin]/(weightrandom*8./7.)) ;
	*/
	fflush(stdout);
	los_low=rupp[rbin];
	
      }					/* next radial bin */

    rlow=rupp[kbin] ;

  } /* next proj-sep bin */

}

/********************************
 * report fatal error and quit 
 */

void die(char errmsg[])  		
{
    fprintf(stderr,"2dcorr> %s\n",errmsg) ;
    fprintf(stderr,"exiting\n") ;
    exit(-1) ;
}
