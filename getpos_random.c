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
#include <stdlib.h>

float ran1(long *idum);

void getposvel_random(float rcube, 
		     float **x, float **y, float **z,
		     float **vx, float **vy, float **vz, int *np, int iseed) 
{
  char a[1000];
  int i,
    np1=1e6 ;
  float x1,x2,x3,x4,x5,x6 ;
  FILE *fp;

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

    if(fp = fopen("iseed","r"))
      {
	fscanf(fp,"%d",&iseed);
	fclose(fp);
      }

    srand48(iseed);

    np1=*np; ;
    for (i=0;i<np1;i++)  {
      (*x)[i]=drand48()*rcube;
      (*y)[i]=drand48()*rcube;
      (*z)[i]=drand48()*rcube;
      (*vx)[i]=0;
      (*vy)[i]=0;
      (*vz)[i]=0;
    }

}
