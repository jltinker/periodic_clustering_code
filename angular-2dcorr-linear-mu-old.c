#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"

#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define pi 3.1415926535898

#define muh(x) fprintf(stderr,"muh %d\n",x)
#define muf(x) fprintf(stderr,"muh %f\n",x)
#define ind(a,b,c) (a)*njack*njack+(b)*njack+(c)

#define HANDHOLD 0
#define NBINLOOKUP 500000     		/* entries in radial bin lookup table */
#define NBINMAX1 201			/* max. number of xi bins, +1 */
#define RBINMAX1 201

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
float func1(float);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
float qtrap(float (*func)(float), float a, float b);

/* Globals
 */
float *xx,*yy,*zz;
int nglobal;

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

  /*char *calloc() ;*/
  FILE *fp ;

  /** TINKER's workspace **/
  float rpi,rsigma,             /* coordinates in (pi,sigma) */
    density2d,                  /* surface density */
    density1d,                  /* linear density */
    los_low,                    /* lower bound of l-o-s bin */
    phi,                        /* angle between sigma,pi */
    *xi_mono,*xi_quad,**xi_mono_jack,**xi_quad_jack,*xi_bar,**xi_bar_jack,
    *xi_qp,*xi_20,**xi_qp_jack,**xi_20_jack,*xi_hexa,**xi_hexa_jack,*xi_wedge1,*xi_wedge2,
    *temp_arr,err_mono[100],err_quad[100],err_qp[100],err_20[100],err_hexa[100];
 
  double  phibar[NBINMAX1][RBINMAX1], /* average value of angle in bins [i][j] */
    pibar[NBINMAX1][RBINMAX1];  /* average value of radial seperation */
  int iredshift=0;              /* flag for distorting z-positions */
  float DELTA_PHI,phi_low,xi_array[NBINMAX1][RBINMAX1];
  int redshift_flag,            /* Output file number */
    velocity_axis;              /* which direction is line-f-sight */

  FILE *jackfile[200];
  char jackname[100];

  FILE *covarfile;
  char covarname[100],aa[100];
  double **covar_mono,**covar_qp,**covar_20;
  float ***xi_jack;
  int njack=5,njack3,ix,iy,iz,ijack,ijack2,iout;

  if (argc < 13 )  {
    fprintf(stderr,"usage: 2dcorr rmin rmax nbin rcube xmin xmax vscale ") ;
    fprintf(stderr,
	    "file1 format1 time1 weight1 los_axis->0,1,2,3 njack [covarfilename]\n");
    exit(0);
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

  sscanf(argv[13],"%d",&njack) ;
  njack3 = njack*njack*njack;
  if(argc>14)
    sscanf(argv[14],"%s",covarname);
  else
    sprintf(covarname,"xi2d.covar");
  
  /* Set parameters for the angular bins
   */
  nrbin=60;
  DELTA_PHI=1.0/nrbin;

  xi_jack = f3tensor(0,njack3-1,0,nbin-1,0,nrbin-1);

  for(i=0;i<njack3;++i)
    for(j=0;j<nbin;++j)
      for(k=0;k<nrbin;++k)
	xi_jack[i][j][k]=0;
  
  covar_mono = dmatrix(0,nbin-1,0,nbin-1);
  covar_qp = dmatrix(0,nbin-1,0,nbin-1);
  covar_20 = dmatrix(0,nbin-1,0,nbin-1);

  for(i=0;i<nbin;++i)
    for(j=0;j<nbin;++j)
      covar_mono[i][j]=
	covar_qp[i][j] = 
	covar_20[i][j] = 0;



  /* Allocate arrays for multipoles.
   */
  xi_mono=vector(0,nbin-1);
  xi_quad=vector(0,nbin-1);
  xi_hexa=vector(0,nbin-1);
  xi_bar=vector(0,nbin-1);
  xi_qp=vector(0,nbin-1);
  xi_20=vector(0,nbin-1);
  xi_wedge1=vector(0,nbin-1);
  xi_wedge2=vector(0,nbin-1);
  xi_mono_jack=matrix(0,njack3-1,0,nbin-1);
  xi_quad_jack=matrix(0,njack3-1,0,nbin-1);
  xi_hexa_jack=matrix(0,njack3-1,0,nbin-1);
  xi_bar_jack=matrix(0,njack3-1,0,nbin-1);
  xi_qp_jack=matrix(0,njack3-1,0,nbin-1);
  xi_20_jack=matrix(0,njack3-1,0,nbin-1);


  /* set the line-of-sight range equal to the input range */
  /*
  zmin=rmin;
  zmax=rmax;
  nrbin=nbin;
  */

  /* read position and velocity data for first object list */
  fp=fopen(argv[8],"r") ;
  sscanf(argv[10],"%f",&time1) ;
  if (fp==NULL)  
    die("couldn't open first input file") ;
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
  redshift_flag = iredshift = 1;
  if(argc>13)
    velocity_axis=atoi(argv[12]);
  else
    velocity_axis=3;
  if(!velocity_axis)
    {
      redshift_flag=0;
      velocity_axis=3;
    }

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
      if (nreset>50)  {
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

  /* Idit's bins are 0.2 in log10
   */
  //  lrstep = rmin;
  //  rmin = pow(10.0,-0.49+lrstep*2);
  //rmax = pow(10.0,-0.49+lrstep*(nbin+1));

  fprintf(stderr,"rmin = %f rmax = %f\n",rmin,rmax);

  for (i=0;i<=NBINLOOKUP;i++)  {
    r=rmax*i/NBINLOOKUP ;
    if (r>0)  {
      //      kbin=(int)floor(log10(r/rmin)/lrstep+1.0) ;
      kbin=(int)floor(log(r/rmin)/lrstep+1.0) ;

      /* if bins are linear use this 
       */
      kbin=(int)((r-rmin)/(rmax-rmin)*nbin);
      
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
	phibar[i][j]=0;
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
  linklist3d(np2,x2,y2,z2,0.0,rcube,rmax,&listbods,&lattice,&nlattice) ;

  fprintf(stderr,"Done with linklist\n");

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

  iout = np1/100;

  /* loop over particles in 1st list */
  for (i=0;i<np1;i++)  {		    
    xp=x1[i] ;
    yp=y1[i] ;
    zp=z1[i] ;

    // print out the progress
    if(i%iout==0)
      fprintf(stderr,"%f\n",i*1./np1);

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
      rfind2_3d(x2,y2,z2,0.0,rcube,rmax,listbods,lattice,nlattice,
	     xp,yp,zp,indx,rsqr,&nrmax,i)  ;
    }  else  {
      rfind_3d(x2,y2,z2,0.0,rcube,rmax,listbods,lattice,nlattice,
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
      if(r<rmin)continue;
      kbin=binlookup[(int)(binfac*r)] ;
      if(kbin>=nbin)continue;
      dz=mabs(z1[i]-z2[j]);      
      if(dz>rhalf)dz=rcube-dz;

      dx=x2[j]-x1[i];
      dy=y2[j]-y1[i];
      if (dx>rhalf)  dx -= rcube ;
      if (dx<-rhalf) dx += rcube ;
      if (dy>rhalf)  dy -= rcube ;
      if (dy<-rhalf) dy += rcube ;
      rsigma=sqrt(dx*dx+dy*dy);

      phi=fabs(dz/r);
      rbin=phi/DELTA_PHI;
      if(rbin==nrbin)rbin--;

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

      if(rbin<nrbin)
	{
	  weightsum2d[kbin][rbin][ioct]+=weight;
	  npair2d[kbin][rbin]++;
	  phibar[kbin][rbin]+=asin(phi);

	  ix = x1[i]*(njack/rcube);
	  iy = y1[i]*(njack/rcube);
	  iz = z1[i]*(njack/rcube);
	  if(ix==njack)ix--;
	  if(iy==njack)iy--;
	  if(iz==njack)iz--;      
	  ijack = ind(ix,iy,iz);
	  
	  ix = x2[j]*(njack/rcube);
	  iy = y2[j]*(njack/rcube);
	  iz = z2[j]*(njack/rcube);
	  if(ix==njack)ix--;
	  if(iy==njack)iy--;
	  if(iz==njack)iz--;
	  ijack2 = ind(ix,iy,iz);
	  
	  if(ijack>=njack3)
	    {
	      fprintf(stderr,"here %d %d %d %d\n",ix,iy,iz,ijack);
	      exit(0);
	    }
	  if(j%2==0)
	    xi_jack[ijack][kbin][rbin]+=weight;
	  else
	    xi_jack[ijack2][kbin][rbin]+=weight;
	}

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
  
  rlow=0;

  for (kbin=0;kbin<nbin;kbin++)  {	/* loop over proj-sep bins */

    /* First get the average value of the radius 
     * (averaged over all the angular bins.
     */
    los_low=0;
    weight0=0.0;
    for (ioct=0;ioct<8;ioct++)  {  
      weight0+=weightsum[kbin][ioct] ;
    }
    if(weight0>0.0) {
      fac=1./weight0;
      rbar[kbin] *= fac;
    }
    else {
      fac=0.0;
      rbar[kbin]=(rupp[kbin]+rlow)*0.5;
    }

    /* loop over the angular bins
     */
    xi_mono[kbin]=xi_quad[kbin]=0;
    xi_wedge1[kbin] = xi_wedge2[kbin] = 0;
    for(i=0;i<njack3;++i)
      xi_mono_jack[i][kbin] = xi_quad_jack[i][kbin] = 0;

    for (rbin=0;rbin<nrbin;rbin++)
      {
	weight0=0.0 ;

	for (ioct=0;ioct<8;ioct++)  {		/* add up values from octants */
	  weight0+=weightsum2d[kbin][rbin][ioct] ;
	}
	if (weight0>0.0)  {
	  fac=1./weight0 ;
	  phibar[kbin][rbin] *= fac ;
	}
	else  {					/* avoid errors in empty bins */
	  fac=0.0 ;
	  phibar[kbin][rbin]=(rbin+0.5)*DELTA_PHI;
	}
    
	/* Volume of wedge of spherical shell 
	 */
	phi_low=rbin*DELTA_PHI;
	vol=2./3.*pi*(rupp[kbin]*rupp[kbin]*rupp[kbin]-rlow*rlow*rlow);
	//	printf("%d %d %f %f\n",rbin,kbin,rupp[kbin],rlow);
	vol*=DELTA_PHI*2.0;
	/*vol*=(cos(phi_low)-cos((phi_low+DELTA_PHI)))*2.0;*/
	weightrandom=np1*avgweight1*density*vol ;	
	xi=weight0/weightrandom-1 ;

	xi_array[kbin][rbin]=xi;

	/* Calculate the monopole and quadrupole for this bin.
	 */
	phi=pi/2 - phibar[kbin][rbin];
	//phi=pi/2 - DELTA_PHI*(rbin+0.5);
   
	xi_mono[kbin]+=xi*DELTA_PHI;
	xi_quad[kbin]+=5*(3*cos(phi)*cos(phi)-1)/2.0*DELTA_PHI*xi; 
	xi_hexa[kbin]+=9*(35*pow(cos(phi),4.0)-30*cos(phi)*cos(phi)+3)/8.*DELTA_PHI*xi; 
	if(rbin<nrbin/2)xi_wedge1[kbin]+=xi*DELTA_PHI/0.5;
	else xi_wedge2[kbin]+=xi*DELTA_PHI/0.5;


	for(i=0;i<njack3;++i)
	  {
	    jackweight = weight0 - xi_jack[i][kbin][rbin];
	    xi = jackweight/(weightrandom*(njack3-1.)/njack3) - 1;
	    xi_mono_jack[i][kbin] += xi*DELTA_PHI;
	    xi_quad_jack[i][kbin] += 5*(3*cos(phi)*cos(phi)-1)/2.0*
	      DELTA_PHI*xi; 
	    xi_hexa_jack[i][kbin] += 9*(35*pow(cos(phi),4.0)-30*cos(phi)*cos(phi)+3)/8.
	      *DELTA_PHI*xi; 
	  }	    

	/*
	printf("%d %d %f %f %e\n",kbin,rbin,rbar[kbin],phi,xi);
	fflush(stdout);
	*/
      }	 /* next angular bin */

    

    /* Calculate diagonal errors on the multipoles
     */
    err_hexa[kbin] = err_mono[kbin] = err_quad[kbin] = 0;
    for(i=0;i<njack3;++i)
      {
	diff = xi_mono[kbin] - xi_mono_jack[i][kbin];
	err_mono[kbin] += diff*diff;
	diff = xi_quad[kbin] - xi_quad_jack[i][kbin];
	err_quad[kbin] += diff*diff;
	diff = xi_hexa[kbin] - xi_hexa_jack[i][kbin];
	err_hexa[kbin] += diff*diff;
      }
    err_mono[kbin] = sqrt(err_mono[kbin]*(njack3-1.)/njack3);
    err_quad[kbin] = sqrt(err_quad[kbin]*(njack3-1.)/njack3);
    err_hexa[kbin] = sqrt(err_hexa[kbin]*(njack3-1.)/njack3);


    /* Output the values of the multipoles.
     */
    
    printf("%e %e %e %e %e %e\n",rbar[kbin],xi_mono[kbin],xi_quad[kbin],xi_hexa[kbin],xi_wedge1[kbin],xi_wedge2[kbin]);
    fflush(stdout);
    
    rlow=rupp[kbin] ;

  } /* next radial bin */
  exit(0);

  xx = vector(1,nbin);
  yy = vector(1,nbin);
  zz = vector(1,nbin);
  nglobal = nbin;

  /* Calculate xi_bar for the overall sample.
   */
  for(i=0;i<nbin;++i)
    {
      xx[i+1] = log(rbar[i]);
      yy[i+1] = log(xi_mono[i]);
    }
  spline(xx,yy,nbin,1.0E+30,1.0E+30,zz);
  for(i=0;i<nbin;++i)
    xi_bar[i] = qtrap(func1,log(0.1),log(rbar[i]))*3/(rbar[i]*rbar[i]*rbar[i] - 0.1*0.1*0.1);


  for(i=0;i<-nbin;++i)
    printf("%e %e %e %e %e\n",rbar[i],xi_quad[i]/(xi_mono[i]-xi_bar[i]),xi_mono[i],xi_quad[i],xi_bar[i]);

  /* Calculate xi_bar for each jackknife sample
   */
  for(j=0;j<njack3;++j)
    {
      for(i=0;i<nbin;++i)
	{
	  xx[i+1] = log(rbar[i]);
	  yy[i+1] = log(xi_mono_jack[j][i]);
	}
      spline(xx,yy,nbin,1.0E+30,1.0E+30,zz);
      for(i=0;i<nbin;++i)
	xi_bar_jack[j][i] = qtrap(func1,log(0.1),log(rbar[i]))*
	  3/(rbar[i]*rbar[i]*rbar[i]-0.1*0.1*0.1);
    }

  /* Calculate covariance matrix for monopole
   */
  for(i=0;i<nbin;++i)
    for(j=0;j<nbin;++j)
      for(k=0;k<njack3;++k)
	covar_mono[i][j] += (xi_mono[i] - xi_mono_jack[k][i])*
	  (xi_mono[j] - xi_mono_jack[k][j])*(njack3-1.)/njack3;

  /* Output covariance matrix
   */
  sprintf(aa,"%s.mono",covarname);
  fp = fopen(aa,"w");
  for(i=0;i<nbin;++i)
    for(j=0;j<nbin;++j)
      fprintf(fp,"%e\n",covar_mono[i][j]);
  fclose(fp);

  /* Calculate the xi_2/xi_0 ratio.
   */
  for(i=0;i<nbin;++i)
    {
      xi_20[i] = xi_quad[i]/xi_mono[i];
      for(j=0;j<njack3;++j)
	xi_20_jack[j][i] = xi_quad_jack[j][i]/xi_mono_jack[j][i];
    }

  /* Calculate covariance matrix for 2/0 ratio
   */
  for(i=0;i<nbin;++i)
    for(j=0;j<nbin;++j)
      for(k=0;k<njack3;++k)
	covar_20[i][j] += (xi_20[i] - xi_20_jack[k][i])*
	  (xi_20[j] - xi_20_jack[k][j])*(njack3-1.)/njack3;

  /* Output covariance matrix
   */
  sprintf(aa,"%s.20",covarname);
  fp = fopen(aa,"w");
  for(i=0;i<nbin;++i)
    for(j=0;j<nbin;++j)
      fprintf(fp,"%e\n",covar_20[i][j]);
  fclose(fp);

  /* Calculate the quadrupole moment
   */
  for(i=0;i<nbin;++i)
    {
      xi_qp[i] = xi_quad[i]/(xi_mono[i] - xi_bar[i]);
      for(j=0;j<njack3;++j)
	xi_qp_jack[j][i] = xi_quad_jack[j][i]/(xi_mono_jack[j][i]-xi_bar_jack[j][i]);
    }
  
  /* Calculate covariance matrix for the quadrupole moment
   */
  for(i=0;i<nbin;++i)
    for(j=0;j<nbin;++j)
      for(k=0;k<njack3;++k)
	covar_qp[i][j] += (xi_qp[i] - xi_qp_jack[k][i])*
	  (xi_qp[j] - xi_qp_jack[k][j])*(njack3-1.)/njack3;

  /* Output covariance matrix
   */
  sprintf(aa,"%s.qp",covarname);
  fp = fopen(aa,"w");
  for(i=0;i<nbin;++i)
    for(j=0;j<nbin;++j)
      fprintf(fp,"%e\n",covar_qp[i][j]);
  fclose(fp);

  /* Output avg values to stdout.
   */
  for(i=0;i<nbin;++i)
    fprintf(stdout,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
	    rbar[i],xi_mono[i],err_mono[i],xi_quad[i],err_quad[i],xi_bar[i],
	    xi_20[i],sqrt(covar_20[i][i]),xi_qp[i],sqrt(covar_qp[i][i]),xi_hexa[i],err_hexa[i]);

}

float func1(float r)
{
  float a;
  splint(xx,yy,zz,nglobal,r,&a);
  r=exp(r);
  return(exp(a)*r*r*r);
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
