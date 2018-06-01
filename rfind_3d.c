/* rfind -- find all particles within rmax of a specified position
   rfind(x,y,z,smin,smax,rmax,linklist,lattice,nlattice,xpos,ypos,zpos,
	 indx,rsqr,nrmax) ;
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
     
     - rfind should be called after a call to "linklist" with the same values
       of the first 9 parameters
     - currently, rfind assumes a periodic boundary condition
*/
#include <stdio.h>
#include <math.h>
#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define pi 3.1415926535898

void rfind_3d(x,y,z,smin,smax,rmax,linklist,lattice,nlattice,xpos,ypos,zpos,indx,
	   rsqr,nrmax) 
float *x, *y, *z ;
float smin, smax, rmax ;
int *linklist, ***lattice ;
int nlattice ;
float xpos, ypos, zpos ;
int *indx ;
float *rsqr ;
int *nrmax ;

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
    iz=(int)(nlattice*(zpos-smin)*sinv) ;
    if (ix>=nlattice||iy>=nlattice||iz>=nlattice ||ix<0||iy<0||iz<0)  {
	fprintf(stderr,"rfind>      error in position or lattice parameters\n");
	fprintf(stderr,"rfind>      nlattice, ix, iy, iz = %d %d %d %d\n",
		nlattice,ix,iy,iz) ;
	fprintf(stderr,"rfind>     position %f %f %f\n",xpos,ypos,zpos) ;
	if (ix>=nlattice)  {
	    ix=nlattice-1 ;
	    fprintf(stderr,"setting ix=%d\n",ix) ;
	}
	if (iy>=nlattice)  {
	    iy=nlattice-1 ;
	    fprintf(stderr,"setting iy=%d\n",iy) ;
	}
	if (iz>=nlattice)  {
	    iz=nlattice-1 ;
	    fprintf(stderr,"setting iz=%d\n",iz) ;
	}
    }

    rmax2=rmax*rmax ;
    nr=0 ;

    for (iix=-1;iix<=1;iix++)
    for (iiy=-1;iiy<=1;iiy++)
    for (iiz=-1;iiz<=1;iiz++)  {
	iiix=(ix+iix+nlattice)%nlattice ;
	iiiy=(iy+iiy+nlattice)%nlattice ;
	iiiz=(iz+iiz+nlattice)%nlattice ;
	p=lattice[iiix][iiiy][iiiz] ;
	while (p>=0)  {
	    dx=mabs(xpos-x[p]) ;
	    dy=mabs(ypos-y[p]) ;
	    dz=mabs(zpos-z[p]) ;
	    if (dx>side2)  dx=side-dx ;
	    if (dy>side2)  dy=side-dy ;
	    if (dz>side2)  dz=side-dz ;
	    r2=dx*dx+dy*dy+dz*dz ;
	    if (r2<=rmax2)  {
		indx[nr]=p ;
		rsqr[nr]=r2 ;
		nr++ ;
	    }
	    if (nr>*nrmax)  {
		fprintf(stderr,"rfind>      too many particles in indx list\n");
		fprintf(stderr,"rfind>      reset nrmax and try again\n");
		exit(-1) ;
	    }
	    p=linklist[p] ;
	}
    }
    *nrmax=nr ;
}
