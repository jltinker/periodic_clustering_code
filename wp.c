#include <stdio.h>
#include <math.h>
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
    pibar[NBINMAX1][RBINMAX1];  /* average value of radial seperation */
  int iredshift=0;              /* flag for distorting z-positions */
  float DELTA_PHI,phi_low,xi_array[NBINMAX1][RBINMAX1];
  int redshift_flag,            /* Output file number */
    velocity_axis;              /* which direction is line-f-sight */
  

  if (argc < 13 )  {
    fprintf(stderr,"usage: 2dcorr rmin rmax nbin rcube xmin xmax vscale ") ;
    fprintf(stderr,
	    "file1 format1 time1 weight1 [output #] [los axis->1,2,3]\n");
    die(" ") ;
  }



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
  zmax=4;  /* <--- SLOAN SPECIFIC NUMBER!!! */
  nrbin=nbin;
  
  /* read position and velocity data for first object list */
  if(!(fp=fopen(argv[8],"r")))
    die("couldn't open first input file") ;
  sscanf(argv[10],"%f",&time1) ;
  if (strcmp(argv[9],"a")==0)  { 
    /*Don't even ask...*/
    fprintf(stderr,"");
    getposvel_ascii(fp,&x1,&y1,&z1,&vx1,&vy1,&vz1,&np1) ;
  }
  else if (strcmp(argv[9],"f")==0)  { 
    getposvel_fastfood(fp,time1,&x1,&y1,&z1,&vx1,&vy1,&vz1,&np1) ;
  }
  else if (strcmp(argv[9],"t")==0)  { 
    getposvel_tipsy(fp,(double)time1,'a',&x1,&y1,&z1,&vx1,&vy1,&vz1,
		    &nsph,&ndark,&nstar,&np1) ;
  }
  else if (strcmp(argv[9],"td")==0)  { 
    getposvel_tipsy(fp,(double)time1,'d',&x1,&y1,&z1,&vx1,&vy1,&vz1,
		    &nsph,&ndark,&nstar,&np1) ;
  }
  else if (strcmp(argv[9],"tg")==0)  { 
    getposvel_tipsy(fp,(double)time1,'g',&x1,&y1,&z1,&vx1,&vy1,&vz1,
		    &nsph,&ndark,&nstar,&np1) ;
  }
  else if (strcmp(argv[9],"ts")==0)  { 
    getposvel_tipsy(fp,(double)time1,'s',&x1,&y1,&z1,&vx1,&vy1,&vz1,
		    &nsph,&ndark,&nstar,&np1) ;
  }
  else {
    die("illegal type for file1") ;
  }
  fclose(fp) ;

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
   * Check to see if we want to do redshift distortions
   */
  if(argc>12)
    {
      redshift_flag=atoi(argv[12]);
      if(redshift_flag>0)iredshift=1;
    }
  if(argc>13)
    velocity_axis=atoi(argv[13]);
  else
    velocity_axis=3;

  /* Change it so that the lin-of-sight direction (the velocity axis)
   * is always the z array.
   */
  if(velocity_axis==1)
    {
      temp_arr=vx1;
      vx1=vz1;
      vz1=temp_arr;

      temp_arr=x1;
      x1=z1;
      z1=temp_arr;
    }
  if(velocity_axis==2)
    {
      temp_arr=vy1;
      vy1=vz1;
      vz1=temp_arr;

      temp_arr=y1;
      y1=z1;
      z1=temp_arr;
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
  rscale=rcube/(xmax-xmin) ;
  rhalf=rcube/2.0 ;
  j=0 ;
  for (i=0;i<np1;i++)  {
    if (weight1[i]>0.0)  {
      x1[j]=rscale*(x1[i]-xmin) ;
      y1[j]=rscale*(y1[i]-xmin) ;
      z1[j]=rscale*(z1[i]-xmin) ;
      vx1[j]=vscale*vx1[i] ;
      vy1[j]=vscale*vy1[i] ;
      vz1[j]=vscale*vz1[i] ;
      weight1[j]=weight1[i] ;
      avgweight1+=weight1[j] ;
      j++ ;
    }
  }
  np1=j ;
  avgweight1/=(float)np1 ;


  if (iauto==0)  {
    j=0 ;
    avgweight2=0.0 ;
    for (i=0;i<np2;i++)  {
      if (weight2[i]>0.0)  {
	x2[j]=rscale*(x2[i]-xmin) ;
	y2[j]=rscale*(y2[i]-xmin) ;
	z2[j]=rscale*(z2[i]-xmin) ;
	vx2[j]=vscale*vx2[i] ;
	vy2[j]=vscale*vy2[i] ;
	vz2[j]=vscale*vz2[i] ;
	weight2[j]=weight2[i] ;
	avgweight2+=weight2[j] ;
	j++ ;
      }
    }
    np2=j ;
    avgweight2/=(float)np2 ;
  }
  else  {
    np2=np1 ;
    avgweight2=avgweight1 ;
  }
  
  /* make the redhsift space distortions in the z-direction */
  if(iredshift)
    {
      redshift_distortions(z1,vz1,rcube,np1,redshift_flag);
      if(!iauto)
	redshift_distortions(z2,vz2,rcube,np2,redshift_flag);
    }

  /* check that no particles lie outside or exactly on boundary */
  cmin=1.e-5 ;
  cmax=0.99998*rcube ;
  nreset=0 ;
  for (i=0;i<np1;i++)  {
    if (x1[i]<=0.0 || x1[i]>=rcube || y1[i]<=0.0 || y1[i]>=rcube || 
	z1[i]<=0.0 || z1[i]>=rcube)  {
      fprintf(stderr,"covar3> warning, resetting coordinates of") ;
      fprintf(stderr," particle %d: %f %f %f\n",i,x1[i],y1[i],z1[i]) ;
      if (x1[i]<=0.0) x1[i]=cmin ;
      if (x1[i]>=rcube) x1[i]=cmax ;
      if (y1[i]<=0.0) y1[i]=cmin ;
      if (y1[i]>=rcube) y1[i]=cmax ;
      if (z1[i]<=0.0) z1[i]=cmin ;
      if (z1[i]>=rcube) z1[i]=cmax ;
      nreset++ ;
      if (nreset>15)  {
	die("too many particles reset, something's probably wrong") ;
      }
    }
  }
  

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
      }
    for (j=0;j<8;j++)  {
      weightsum[i][j]=0.0 ;
      vstream[i][j]=0.0 ;
      vpar2[i][j]=0.0 ;
      vtan2[i][j]=0.0 ;
      vxi2[i][j]=0.0 ;
    }
  }

  /* make linklist */
  if (HANDHOLD) 
    fprintf(stderr,"covar3>  calling linklist\n") ;
  linklist(np2,x2,y2,z2,0.0,rcube,rmax,&listbods,&lattice,&nlattice) ;

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

    /* Find particles within rmax.  For autocorrelation, find only particles
     * with index j>i (using rfind2); for cross-correlation, find all particles.
     *
     */
    nrmax=np2 ;
    if (iauto==1)  {	
      rfind2(x2,y2,z2,0.0,rcube,rmax,listbods,lattice,nlattice,
	     xp,yp,zp,indx,rsqr,&nrmax,i)  ;
    }  else  {
      rfind(x2,y2,z2,0.0,rcube,rmax,listbods,lattice,nlattice,
	    xp,yp,zp,indx,rsqr,&nrmax)  ;
    }


    for (k=0;k<nrmax;k++)  {

      /* find radial separation bin for this particle pair */
      j=indx[k] ;

      /* eliminate pairs with identical positions so that the cross-correlation
       * of a distribution with itself is the same as the autocorrelation,
       * and to prevent numerical errors in velocity computations 
       */
      if (rsqr[k]==0.0)  {	
	continue ;
      }

      r=sqrt(rsqr[k]) ;
      kbin=binlookup[(int)(binfac*r)] ;
      dz=mabs(z1[i]-z2[j]);
      if(dz>rhalf)dz=rcube-dz;
      if(dz>zmax)continue;

      /* for an autocorrelation, choose the octant based on particle i half the
       * time (when j is odd) and particle j half the time, to avoid biases
       * in number of pairs per octant caused by only looking at pairs with
       * j > i 
       */
      if (iauto==1 && j%2==0)  {
	ioct=0 ;
	if (x2[j]>rhalf)  ioct++ ;
	if (y2[j]>rhalf)  ioct+=2 ;
	if (z2[j]>rhalf)  ioct+=4 ;
      }
      else  {
	ioct=ioct0 ;
      }
      
      /* add weight to total, which will be used for xi(r) and weighting vel stats */
      weight=weight1[i]*weight2[j] ;
      weightsum[kbin][ioct]+=weight ;
      rbar[kbin]+=weight*r ;
      npair[kbin]++ ;
      
    }
  }



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
    vol=pi*(rupp[kbin]*rupp[kbin]-rlow*rlow)*2*zmax ;
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
     */
    printf("%7.4f %7.4f %7.4f %10.3e %10.3e %10d %10.0f\n",
	   rlow,rupp[kbin],rbar[kbin],xi,errxi,npair[kbin],weightrandom*8.0/7.0) ;
    rlow=rupp[kbin] ;
  }					/* next radial bin */



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
