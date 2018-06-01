/* rfind2 -- find all particles within rmax of a specified position
   rfind2(x,y,z,smin,smax,rmax,linklist,lattice,nlattice,xpos,ypos,zpos,
	 indx,rsqr,nrmax,pmin) ;
     * x, y, z = arrays of particle x, y, z coordinates
     * smin, smax: particles are assumed to occupy a box running from
			(smin,smin,smin) to (smax,smax,smax)
     * rmax = maximum radius within which to identify neighboring particles
     * linklist = particle link list, produced by "linklist"
     * lattice = lattice pointer array, produced by "linklist"
     * nlattice = dimension of lattice array in each dimension
     * xpos, ypos, zpos = coordinates defining center of search (within box)
     * indx = on return, contains list of particles within rmax of search center
	      allocate as (int *)calloc(nrmax,sizeof(int))
     * rsqr = on return, rsqr[i] is the squared distance of particle indx[i]
	      allocate as (float *)calloc(nrmax,sizeof(float))
     * nrmax = on input: dimension of indx[] and rsqr[], for error checking
	       on return: number of particles within rmax
     * pmin = minimum particle index; only particles with p>pmin will be
	      returned
     
     - rfind2 should be called after a call to "linklist" with the same values
       of the first 9 parameters
     - currently, rfind2 assumes a periodic boundary condition

     Version 2 includes a minimum particle index, to be used in calculations
     like the correlation function where there is a symmetry over pairs
*/
#include <stdio.h>
#include <math.h>
#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define pi 3.1415926535898

void rfind2(x,y,z,smin,smax,rmax,linklist,lattice,nlattice,xpos,ypos,zpos,indx,
	   rsqr,nrmax,pmin) 
float *x, *y, *z ;
float smin, smax, rmax ;
int *linklist, **lattice ;
int nlattice ;
float xpos, ypos, zpos ;
int *indx ;
float *rsqr ;
int *nrmax ;
int pmin ;

{
    int ix,iy,iz,
	iix,iiy,iiz,
	iiix,iiiy,iiiz,
	nr,p ;
    float rmax2,r2,
	  side,side2,sinv,
	  dx,dy,dz ;

    side=(smax-smin) ;
    side2=(smax-smin)/2. ;
    sinv=1./side ;

    ix=(int)(nlattice*(xpos-smin)*sinv) ;
    iy=(int)(nlattice*(ypos-smin)*sinv) ;
    if (ix>=nlattice || iy>=nlattice || ix<0 || iy<0 )  {
	fprintf(stderr,"rfind>  error in position or lattice parameters\n");
	fprintf(stderr,"rfind>  nlattice, ix, iy  = %d %d %d %f %f\n",
		nlattice,ix,iy,xpos,ypos) ;
	exit(-1) ;
    }

    rmax2=rmax*rmax ;
    nr=0 ;

    for (iix=-1;iix<=1;iix++)
      for (iiy=-1;iiy<=1;iiy++) {
	iiix=(ix+iix+nlattice)%nlattice ;
	iiiy=(iy+iiy+nlattice)%nlattice ;
	p=lattice[iiix][iiiy] ;
	while (p>=0)  {
	    if (p>pmin)  {
		dx=mabs(xpos-x[p]) ;
		dy=mabs(ypos-y[p]) ;
		if (dx>side2)  dx=side-dx ;
		if (dy>side2)  dy=side-dy ;
		r2=dx*dx+dy*dy ;
		if (r2<=rmax2)  {
		    indx[nr]=p ;
		    rsqr[nr]=r2 ;
		    nr++ ;
		}
		if (nr>*nrmax)  {
		    fprintf(stderr,"rfind>  too many particles in indx list\n");
		    fprintf(stderr,"rfind>  reset nrmax and try again\n");
		    exit(-1) ;
		}
	    }
	    p=linklist[p] ;
	}
    }
    *nrmax=nr ;
}
