/* getposvel_fastfood -- get positions and velocities from fastfood binary file
   void getposvel_fastfood(fp,zstep,*x,*y,*z,*vx,*vy,*vz,np) 
     * fp = file pointer to binary file
     * zstep = redshift of desired output 
     * x,y,z = pointers to position arrays
     * vx,vy,vz = pointers to velocity arrays
     * np = number of objects with returned positions
*/

#include <stdio.h>
#include <math.h>
#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define pi 3.1415926535898

void getposvel_fastfood(FILE *fp, float zstep, 
			float **x, float **y, float **z,
			float **vx, float **vy, float **vz, int *np) 
{
    int i,
	np1,
	ibyte,
	narray=6,
	idat[5] ;
    float znow,
	  xdat[9] ;

    ftread(idat,sizeof(int),5,fp) ;		/* read header data */
    ftread(xdat,sizeof(float),9,fp) ;

    *np = idat[1] ;				

    np1= *np ;

    *x=(float *)malloc(np1*sizeof(float)) ;
    *y=(float *)malloc(np1*sizeof(float)) ;
    *z=(float *)malloc(np1*sizeof(float)) ;
    if (*x==NULL || *y==NULL || *z==NULL)  {
        fprintf(stderr,"getposvel_fastfood> no memory for position arrays\n") ;
	exit(-1) ;
    }

    *vx=(float *)malloc(np1*sizeof(float)) ;
    *vy=(float *)malloc(np1*sizeof(float)) ;
    *vz=(float *)malloc(np1*sizeof(float)) ;
    if (*vx==NULL || *vy==NULL || *vz==NULL)  {
        fprintf(stderr,"getposvel_fastfood> no memory for velocity arrays\n") ;
	exit(-1) ;
    }

    while (fread(&ibyte,sizeof(int),1,fp) == 1)  {	/* read next timestep */
	fread(&znow,sizeof(float),1,fp) ;
	fread(&ibyte,sizeof(int),1,fp) ;
	if (mabs(zstep-znow)>0.05)  {
	    for (i=0;i<narray;i++)
	        ftread(*x,sizeof(float),np1,fp) ;
	    continue ;
	}
	if (mabs(zstep-znow)>0.05)  {
	    fprintf(stderr,"getposvel_fastfood> using redshift %f\n",znow) ;
	}
	ftread(*x,sizeof(float),np1,fp) ;
	ftread(*y,sizeof(float),np1,fp) ;
	ftread(*z,sizeof(float),np1,fp) ;
	ftread(*vx,sizeof(float),np1,fp) ;
	ftread(*vy,sizeof(float),np1,fp) ;
	ftread(*vz,sizeof(float),np1,fp) ;
    }
}
