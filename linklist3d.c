/* linklist -- create a linked list for efficient neighbor searching
   linklist(np,x,y,z,smin,smax,rmax,listbods,lattice,nlattice)
     * np = number of particles
     * x, y, z = arrays of particle x, y, z coordinates
     * smin, smax: particles are assumed to occupy a box running from
		   (smin,smin,smin) to (smax,smax,smax)
     * rmax = maximum radius that will be searched for neighbors (used to 
	      determine size of lattice mesh cells
     * listbods = on return, pointer to the particle link list
     * lattice = on return, pointer to array of pointers to the first 
			    particle in each cell
     * nlattice = on return, dimension of lattice array in each dimension
*/
#include <stdio.h>
#include <math.h>
#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define pi 3.1415926535898

#define	ALLOC3D(type,dim0,dim1,dim2) \
	(type***)a3alloc((unsigned)dim0,(unsigned)dim1,(unsigned)dim2,(unsigned)sizeof(type))
#define NLATMAX 1200			/* maximum lattice dimension */
int ***i3tensor_2(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

void linklist3d(np,x,y,z,smin,smax,rmax,listbods,lattice,nlattice)
int np ;
float *x, *y, *z ;
float smin, smax, rmax ;
int **listbods, ****lattice ;
int *nlattice ;

{
    int i,ix,iy,iz,
	j,k,
	nmesh ;
    float sinv ;
    char *calloc() ;

    nmesh=(int)((smax-smin)/rmax) ;
    if (nmesh<4)  {
	fprintf(stderr,
		"linklist>   nlattice = %d, why use linklist?\n",nmesh) ;
	/*exit(-1) ;*/
    }
    if (nmesh>NLATMAX)  nmesh=NLATMAX ;
    /*
    *lattice=ALLOC3D(int,nmesh,nmesh,nmesh) ;
    */
    *lattice=(int ***)i3tensor_2(0,nmesh-1,0,nmesh-1,0,nmesh-1);
    *listbods=(int *)calloc(np,sizeof(int)) ;
    for (i=0;i<np;i++)  (*listbods)[i]=-1 ;
    for (i=0;i<nmesh;i++)
    for (j=0;j<nmesh;j++)
    for (k=0;k<nmesh;k++)
	(*lattice)[i][j][k]=-1 ;

    sinv=1./(smax-smin) ;
    for (i=0;i<np;i++)  {
	ix=(int)(nmesh*(x[i]-smin)*sinv) ;
	iy=(int)(nmesh*(y[i]-smin)*sinv) ;
	iz=(int)(nmesh*(z[i]-smin)*sinv) ;
	if (ix>nmesh-1)  ix=nmesh-1 ;	 /* this shouldn't happen, but . . . */
	if (iy>nmesh-1)  iy=nmesh-1 ;
	if (iz>nmesh-1)  iz=nmesh-1 ;
	(*listbods)[i]=(*lattice)[ix][iy][iz] ;
	(*lattice)[ix][iy][iz]=i ;
    }
    *nlattice=nmesh ;
}
