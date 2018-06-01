/* covar3 -- compute covariance function or cross-correlation and velocity stats
   covar3  rmin rmax nbin rcube xmin xmax vscale file1 format1 time1 weight1 
           file2 format2 time2 weight2
     * rmin = outer radius of smallest bin (in Mpc/h)
     * rmax = outer radius of largest bin
     * nbin = number of bins (logarithmically spaced in r)
     * rcube = size of cube in Mpc/h
     * xmin, xmax = min. and max. coordinate values of particles
     * vscale = scaling from velocity units to km/s
     * file1 = name of first object file
     * format1 = format of first object file: 
		 f=fastfood
                 a=ascii lines with x y z vx vy vz, no header lines
		 t=tipsy, all particles
		 td=tipsy, dark particles only
		 tg=tipsy, gas particles only
		 ts=tipsy, star particles only
     * time1 = desired output step (only relevant for t or f)
     * weight1 = file containing weights for objects in list 1, a straight
                 ascii list, 1 line per object; if weight is given as '1'
		 then no file is read and all objects are assigned equal weight
     * file2 etc. : if file2 is given as "auto", then an autocorrelation 
                    analysis is performed for objects in list 1; otherwise, 
		    file2, format2, time2, and weight2 give the analogous 
		    information for objects in a second list, and the 
		    cross-correlation is computed (though particles with
                    identical positions are not included in a cross-correlation)
     > rlow, rupp, rbar, xi, err, vstream, err, vpar2, err, vtan2, err, 
         vxi, err, npair
       - rlow = lower limit of radius in bin
       - rupp = upper limit of radius in bin 
       - rbar = mean pair separation in bin
       - xi = value of two-point correlation function
       - vstream = mean pairwise streaming motion 
       - vpar2 = dispersion about mean streaming motion along line of sight
       - vtan2 = pairwise dispersion tangential to line of sight (reduced
		 to 1-d, i.e. divided by root(2))
       - vxi = dot product velocity correlation 
       - npair = number of pairs in bin
       :: radii are in Mpc/h (or whatever the units of rcube are)
          velocities are in km/s (or whatever units are given by vscale)
	  vpar2, vtan2, and vxi are square-rooted -- units are km/s not (km/s)^2
	  error bars are calculated via jackknife using the eight octants of
	    the cube
*/
#include <stdio.h>
#include <math.h>
#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define pi 3.1415926535898

#define HANDHOLD 0
#define NBINLOOKUP 50000     		/* entries in radial bin lookup table */
#define NBINMAX1 101			/* max. number of xi bins, +1 */

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

main(argc,argv)
int argc ;
char *argv[] ;

{
    int np1,np2,		/* particle numbers in 2 samples */
	nlattice,		/* chaining mesh lattice size */
	i, j, k,		/* assorted integers */
	nbin,			/* number of radial bins */
	ibin,kbin,		/* bin counters */
	nsph,ndark,nstar,	/* for use in tipsy binary reads */
	iauto,			/* 1->autocorrelation, 0->cross-correlation */
	ioct,ioct0,		/* octant counter */
        nreset,			/* number of particles with reset coordinates */
	nrmax ;			/* number of particles within rmax */

    float rmin,rmax,rmax2,	/* radial range, and square of rmax */
	  rcube,		/* size of cube in Mpc/h */
	  rhalf,		/* half of rcube */
	  xmin,xmax,		/* range of particle coordinates */
	  vscale,		/* velocity scaling factor */
	  time1,time2,		/* times to read from input files */
	  rscale,		/* scaling factor for positions */
	  avgweight1,		/* mean weight of particles in set 1 */
	  avgweight2,		/* mean weight of particles in set 2 */
	  cmin,cmax,		/* used to reset particles on boundary */
	  lrstep,		/* logarithmic radial step */
	  binfac,		/* used in radial bin calculation */
	  xp,yp,zp,		/* position of reference particle */
	  vxp,vyp,vzp,		/* velocity of reference particle */
	  dx,dy,dz,		/* difference in position */
	  dvx,dvy,dvz,		/* difference in velocity */
	  r,   			/* radial separation */
	  v12,			/* line-of-sight velocity */
	  rlow,			/* lower radius of bin */
	  density,		/* mean weighted particle density in cube */
	  weight0,		/* sum of weight over octants */
	  vstream0,		/* sum of vstream over octants */
	  vpar0,		/* sum of vpar2 over octants */
	  vtan0,		/* sum of vtan2 over octants */
	  vxi0,			/* sum of vxi2 over octants */
	  fac,			/* 1/(summed weight in bin) */
	  vol,			/* volume of radial bin */
	  weightrandom,		/* weight expected for random distribution */
	  xi,			/* correlation function */
	  ovstream,		/* output value of streaming velocity */
	  ovpar,		/* output value of parallel dispersion */
	  ovtan,		/* output value of tangential dispersion */
	  ovxi,			/* output value of velocity correlation */
	  errxi,		/* error in correlation function */
	  errvstream,		/* error in streaming velocity */
	  errvpar,		/* error in parallel dispersion */
	  errvtan,		/* error in tangential dispersion */
	  errvxi,		/* error in velocity correlation */
	  jackweight,		/* weight in jackknife subsample */
	  diff,			/* difference of jackknife value from mean */
	  vstream1,		/* streaming velocity in jackknife sample */
	  vxi1 ;		/* velocity correlation in jackknife sample */

    double weight ;		/* product of weights of two particles */

    int *listbods, ***lattice, 	/* linklist arrays */
        *indx,			/* index array for particles within rmax */
        *binlookup ;		/* radial bin given separation */

    float *rupp,		/* upper radius of bin */
          *rsqr,		/* squared distance for particles within rmax */
          *x1, *y1, *z1, 	/* positions for first particle set */
	  *vx1, *vy1, *vz1,     /* velocities for first particle set */
	  *weight1,		/* weights for first particle set */
          *x2, *y2, *z2, 	/* positions for second particle set */
	  *vx2, *vy2, *vz2,	/* velocities for second particle set */
	  *weight2 ;		/* weights for second particle set */

/* note that second subscript in subsequent arrays is for values in 
 * individual octants 
 */
    int npair[NBINMAX1] ;		/* number of pairs in radial bin */
    double rbar[NBINMAX1],		/* mean separation in radial bin */
           weightsum[NBINMAX1][8],      /* summed weight product */
	   vstream[NBINMAX1][8],	/* mean streaming velocity */
	   vpar2[NBINMAX1][8],		/* parallel pairwise dispersion */
	   vtan2[NBINMAX1][8],		/* tangential pairwise dispersion */
	   vxi2[NBINMAX1][8] ;		/* velocity correlation function */

    char *calloc() ;
    FILE *fp ;

    if (argc != 13 && argc != 16)  {
	fprintf(stderr,"usage: covar3 rmin rmax nbin rcube xmin xmax vscale ") ;
	fprintf(stderr,
                   "file1 format1 time1 weight1 file2 format2 time2 weight2\n");
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

    if (strcmp(argv[12],"auto")==0)  {		/* auto-correlation */
	iauto=1 ;
	x2=x1 ;
	y2=y1 ;
	z2=z1 ;
	vx2=vx1 ;
	vy2=vy1 ;
	vz2=vz1 ;
	weight2=weight1 ;
	np2=np1 ;
    }
    else  {		/* read data and weights for second object list */
	iauto=0 ;
	fp=fopen(argv[12],"r") ;
	sscanf(argv[14],"%f",&time2) ;
	if (fp==NULL)  
	    die("couldn't open second input file") ;
	if (strcmp(argv[13],"a")==0)  { 
	    getposvel_ascii(fp,&x2,&y2,&z2,&vx2,&vy2,&vz2,&np2) ;
	}
	else if (strcmp(argv[13],"f")==0)  { 
	    getposvel_fastfood(fp,time2,&x2,&y2,&z2,&vx2,&vy2,&vz2,&np2) ;
	}
	else if (strcmp(argv[13],"t")==0)  { 
	    getposvel_tipsy(fp,(double)time2,'a',&x2,&y2,&z2,&vx2,&vy2,&vz2,
			    &nsph,&ndark,&nstar,&np2) ;
	}
	else if (strcmp(argv[13],"td")==0)  { 
	    getposvel_tipsy(fp,(double)time2,'d',&x2,&y2,&z2,&vx2,&vy2,&vz2,
			    &nsph,&ndark,&nstar,&np2) ;
	}
	else if (strcmp(argv[13],"tg")==0)  { 
	    getposvel_tipsy(fp,(double)time2,'g',&x2,&y2,&z2,&vx2,&vy2,&vz2,
			    &nsph,&ndark,&nstar,&np2) ;
	}
	else if (strcmp(argv[13],"ts")==0)  { 
	    getposvel_tipsy(fp,(double)time2,'s',&x2,&y2,&z2,&vx2,&vy2,&vz2,
			    &nsph,&ndark,&nstar,&np2) ;
	}
	else {
	    die("illegal type for file2") ;
	}
	fclose(fp) ;

	weight2 = (float *) calloc(np2,sizeof(float)) ;	
	if (strcmp(argv[15],"1")!=0)  {
	    fp=fopen(argv[15],"r") ;
	    if (fp==NULL)
		die("couldn't open second weight file") ;
	    for (i=0;i<np2;i++)  {
		if (fscanf(fp,"%f",&weight2[i])==EOF)  {
		    die("premature end of second weight file") ;
		}
	    }
	    fclose(fp) ;
	}
	else  {
	    for (i=0;i<np2;i++)  {
		weight2[i]=1.0 ;
	    }
	}
    }

/* rescale positions and velocities, culling weight=0 objects along the 
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
    if (iauto==0)  {
	for (i=0;i<np2;i++)  {
	    if (x2[i]<=0.0 || x2[i]>=rcube || y2[i]<=0.0 || y2[i]>=rcube || 
	      z2[i]<=0.0 || z2[i]>=rcube)  {
		fprintf(stderr,"covar3> warning, resetting coordinates of") ;
		fprintf(stderr," particle %d: %f %f %f\n",i,x2[i],y2[i],z2[i]) ;
	        if (x2[i]<=0.0) x2[i]=cmin ;
		if (x2[i]>=rcube) x2[i]=cmax ;
		if (y2[i]<=0.0) y2[i]=cmin ;
		if (y2[i]>=rcube) y2[i]=cmax ;
		if (z2[i]<=0.0) z2[i]=cmin ;
		if (z2[i]>=rcube) z2[i]=cmax ;
	    }
	}
    }


/* Create lookup table so that we can go from r to the appropriate bin
 * without taking a log.  The bins have a nearly but not exactly logarithmic
 * spacing.  rupp[ibin] is the largest radius that goes into bin ibin (more
 * exactly, it's the smallest radius that goes into bin ibin+1).
 * In non-exhaustive tests, this is about 14% faster than taking the log
 * to get exactly logarithmic bins, as done by covar3b.
 */
    rupp = (float *) calloc(nbin+1,sizeof(float)) ;
    lrstep=log(rmax/rmin)/(float)(nbin-1) ;
    binlookup=(int *)calloc(NBINLOOKUP+2,sizeof(int)) ;
    ibin=0 ;
    for (i=0;i<=NBINLOOKUP;i++)  {
	r=rmax*i/NBINLOOKUP ;
        if (r>0)  {
	    kbin=(int)floor(log(r/rmin)/lrstep+1.0) ;
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

/* do the work */
    for (i=0;i<np1;i++)  {		/* loop over particles in 1st list */

        xp=x1[i] ;
	yp=y1[i] ;
	zp=z1[i] ;
	vxp=vx1[i] ;
	vyp=vy1[i] ;
	vzp=vz1[i] ;

/* compute octant index based on particle position */
	ioct0=0 ;
	if (xp>rhalf)  ioct0++ ;
	if (yp>rhalf)  ioct0+=2 ;
	if (zp>rhalf)  ioct0+=4 ;

/* Find particles within rmax.  For autocorrelation, find only particles
 * with index j>i (using rfind2); for cross-correlation, find all particles.
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

/* compute contribution to velocity statistics */
	    dx=x2[j]-xp ;
	    dy=y2[j]-yp ;
	    dz=z2[j]-zp ;
	    if (dx>rhalf)  dx -= rcube ;
	    if (dx<-rhalf) dx += rcube ;
	    if (dy>rhalf)  dy -= rcube ;
	    if (dy<-rhalf) dy += rcube ;
	    if (dz>rhalf)  dz -= rcube ;
	    if (dz<-rhalf) dz += rcube ;		
	    dvx=vx2[j]-vxp ;
	    dvy=vy2[j]-vyp ;
	    dvz=vz2[j]-vzp ;
	    v12=(dvx*dx+dvy*dy+dvz*dz)/r ;
	    vstream[kbin][ioct]+=weight*v12 ;
	    vpar2[kbin][ioct]+=weight*v12*v12 ;
	    vtan2[kbin][ioct]+=weight*(dvx*dvx+dvy*dvy+dvz*dvz-v12*v12) ; 
	    vxi2[kbin][ioct]+=
	             weight*(vx2[j]*vxp+vy2[j]*vyp+vz2[j]*vzp) ;

        }				/* next inner loop particle */

    }					/* next outer loop particle */

/* output general info */
    if (iauto==1)  {
         printf("# auto: np = %d ; rcube = %f ; vscale = %f\n",
		 np1,rcube,vscale) ;
    }
    else  {
         printf("# cross: np1 = %d ; np2 = %d ; rcube = %f ; vscale = %f\n",
		 np1,np2,rcube,vscale) ;
    }

/* compute the output quantities and error bars */
    rlow=0.0 ;
    density=np2*avgweight2/(rcube*rcube*rcube) ;
    if (iauto==1)  {		/* pairs not double counted */
        density=density/2. ;
    }

    for (kbin=0;kbin<nbin;kbin++)  {	/* loop over radial bins */
	weight0=0.0 ;
	vstream0=0.0 ;
	vpar0=0.0 ;
	vtan0=0.0 ;
	vxi0=0.0 ;

	for (ioct=0;ioct<8;ioct++)  {		/* add up values from octants */
	    weight0+=weightsum[kbin][ioct] ;
	    vstream0+=vstream[kbin][ioct] ;
	    vpar0+=vpar2[kbin][ioct] ;
	    vtan0+=vtan2[kbin][ioct] ;
	    vxi0+=vxi2[kbin][ioct] ;
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
	vol=4.*pi/3.*(rupp[kbin]*rupp[kbin]*rupp[kbin]-rlow*rlow*rlow) ;
	weightrandom=np1*avgweight1*density*vol ;	
	xi=weight0/weightrandom-1 ;

/* compute output values of velocity statistics */
	ovstream = vstream0*fac ;
	ovpar=sqrt(vpar0*fac-ovstream*ovstream) ;
	ovtan=sqrt(vtan0*fac/2.) ;
	if (vxi0>0.0)  {
	    ovxi=sqrt(vxi0*fac) ;
	}
	else  {
	    ovxi= -sqrt(-vxi0*fac) ;
	}

/* compute the jackknife errors by subtracting each octant from the total
 * in turn, computing the difference between this jackknife subsample (seven
 * of the eight octants) and the result for the whole cube, and summing
 * the differences in quadrature.
 */
	weightrandom *= 7.0/8.0 ;	/* each jack sample is 7/8 of volume */
	errxi=0.0 ;
	errvstream=0.0 ;
	errvpar=0.0 ;
	errvtan=0.0 ;
	errvxi=0.0 ;
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
	    vstream1=(vstream0-vstream[kbin][ioct])*fac ;
	    diff=vstream1-ovstream ;
	    errvstream+=diff*diff ;
	    diff=sqrt((vpar0-vpar2[kbin][ioct])*fac-vstream1*vstream1*0.9999)
                 -ovpar ;    /* factor 0.9999 prevents sqrt(negative) error */
	    errvpar+=diff*diff ;
	    diff=sqrt((vtan0-vtan2[kbin][ioct])*fac/2.)-ovtan ;
	    errvtan+=diff*diff ;
	    vxi1=vxi0-vxi2[kbin][ioct] ;
	    if (vxi1>=0.0)  {
	        diff=sqrt(vxi1*fac)-ovxi ;
	    }
	    else  {
	        diff=-sqrt(-vxi1*fac)-ovxi ;
	    }
	    errvxi+=diff*diff ;
	}

	errxi=sqrt(errxi*8./7.) ;
	errvstream=sqrt(errvstream*8./7.) ;
	errvpar=sqrt(errvpar*8./7.) ;
	errvtan=sqrt(errvtan*8./7.) ;
	errvxi=sqrt(errvxi*8./7.) ;

	printf("%7.4f %7.4f %7.4f %10.3e %10.3e %7.4f %6.4f",
	       rlow,rupp[kbin],rbar[kbin],xi,errxi,ovstream,errvstream) ;
	printf(" %7.4f %6.4f %7.4f %6.4f %7.4f %6.4f %10d\n",
	       ovpar,errvpar,ovtan,errvtan,ovxi,errvxi,npair[kbin]) ;

	rlow=rupp[kbin] ;
    }					/* next radial bin */

}


void die(char errmsg[])  		/* report fatal error and quit */
{
    fprintf(stderr,"covar3> %s\n",errmsg) ;
    fprintf(stderr,"exiting\n") ;
    exit(-1) ;
}
