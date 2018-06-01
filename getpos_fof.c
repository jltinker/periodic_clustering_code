/* getposvel_ascii -- get positions and velocities from an ascii file
   void getposvel_ascii(fp,*x,*y,*z,*vx,*vy,*vz,np) 
     * fp = file pointer to ascii file
     * x,y,z = pointers to position arrays
     * vx,vy,vz = pointers to velocity arrays
     * np = number of objects with returned positions
     - assumed format for ascii file is 
	x y z vx vy vz
       one object per line, no header or trailing lines
*/

#include <stdio.h>
#include <math.h>

void getposvel_ascii(FILE *fp, 
		     float **x, float **y, float **z,
		     float **vx, float **vy, float **vz, int *np) 
{
    int i,
	np1 ;
    float xdum,ydum,zdum,vxdum,vydum,vzdum ;

    np1=0 ;
    while(fscanf(fp,"%f %f %f %f %f %f",&xdum,&ydum,&vxdum,&vydum,&vzdum)
          !=EOF) {
	np1++ ;
    }
    rewind(fp) ;

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

    *np=np1 ;
    for (i=0;i<np1;i++)  {
	fscanf(fp,"%f %f %f %f %f %f",
	       &(*x)[i],&(*y)[i],&(*z)[i],&(*vx)[i],&(*vy)[i],&(*vz)[i]) ;
    }
}
