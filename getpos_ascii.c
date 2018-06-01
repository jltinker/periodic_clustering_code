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

void getposvel_ascii(FILE *fp, 
		     float **x, float **y, float **z,
		     float **vx, float **vy, float **vz, int *np) 
{
  char a[1000];
  int i,
    np1=0 ;
  float x1,x2,x3,x4,x5,x6 ;
    
    while(!feof(fp))
      {
	fgets(a,1000,fp);
	np1++;
      }
    np1--;

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
      fscanf(fp,"%f %f %f %f %f %f",&x1,&x2,&x3,&x4,&x5,&x6);
      (*x)[i]=x1;
      (*y)[i]=x2;
      (*z)[i]=x3;
      (*vx)[i]=x4;
      (*vy)[i]=x5;
      (*vz)[i]=x6;
      fgets(a,1000,fp);
    }

}
