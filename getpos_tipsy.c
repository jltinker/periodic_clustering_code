/* getpos_tipsy -- get positions from a tipsy binary file
   void getpos_tipsy(fp,time,type,*x,*y,*z,nsph,ndark,nstar,np) 
     * fp = file pointer to binary file
     * time = desired output time
     * type = d, g, s, b, a for dark, gas, star, baryon, or all positions
   variables set on return:
     * x,y,z = pointers to position arrays
     * nsph = number of sph particles in tipsy file
     * ndark = number of dark particles in tipsy file
     * nstar = number of sph particles in tipsy file
     * np = number of objects with returned positions
*/

#include <stdio.h>
#include <math.h>
#include "tipsydefs.h"

void getpos_tipsy(FILE *fp, double time, char type, 
                 float **x, float **y, float **z,
                 int *nsph, int *ndark, int *nstar, int *np) 
{
    int i ;
    struct gas_particle *gas_particles, *gp, *lastgp;
    struct dark_particle *dark_particles, *dp, *lastdp;
    struct star_particle *star_particles, *sp, *lastsp;
    struct dump header;
    double currtime = 0.0;
    long currpos = 0L ;
    long lastpos = 0L ;
    int np1 ;
    char *calloc() ;

    if ((float)currtime > (float)time){
	fseek(fp,0L,0);
	currtime=0.0;
	currpos=0;
    }

/* find file position for this output time */
    while(1)  {			
       if (fread((char *)&header,sizeof(header),1,fp) != 1) {
	   fprintf(stderr,"getpos_tipsy> time too large, using %f\n",
		          (float)currtime) ;
	    break ;
        }
	currtime = header.time ;
	currpos = ftell(fp) - sizeof(header);
	if ( (float)header.time >= (float)time ) 
	    break ;
	fseek(fp,
	      sizeof(gas_particles[0])*header.nsph +
	      sizeof(dark_particles[0])*header.ndark +
	      sizeof(star_particles[0])*header.nstar,
	      1) ;
    }

/* read the particle data */
    fseek(fp,currpos,0) ;
    fread((char *)&header,sizeof(header),1,fp) ;
    if(header.nsph != 0) {
	gas_particles = (struct gas_particle *)
			    malloc(header.nsph*sizeof(*gas_particles));
	if(gas_particles == NULL) {
	    fprintf(stderr,"getpos_tipsy> no memory for gas particles\n") ;
	    exit(-1) ;
	}
    }
    if(header.ndark != 0) {
	dark_particles = (struct dark_particle *)
			    malloc(header.ndark*sizeof(*dark_particles));
	if(dark_particles == NULL) {
	    fprintf(stderr,"getpos_tipsy> no memory for dark particles\n") ;
	    exit(-1) ;
	}
    }
    if(header.nstar != 0) {
	star_particles = (struct star_particle *)
			    malloc(header.nstar*sizeof(*star_particles));
	if(star_particles == NULL) {
	    fprintf(stderr,"getpos_tipsy> no memory for star particles\n") ;
	    exit(-1) ;
	}
    }

    fread((char *)gas_particles,sizeof(struct gas_particle),
		     header.nsph,fp) ;
    fread((char *)dark_particles,sizeof(struct dark_particle),
		     header.ndark,fp) ;
    fread((char *)star_particles,sizeof(struct star_particle),
		     header.nstar,fp) ;
    currpos = lastpos ;
    fseek(fp,currpos,0) ;
    currtime = header.time ;
    if ((float)time != (float)currtime){
	fprintf(stderr,"getpos_tipsy> used time %f\n",(float)currtime);
    }

    lastgp = gas_particles + header.nsph ;
    lastdp = dark_particles + header.ndark ;
    lastsp = star_particles + header.nstar ;

    *ndark = header.ndark ;
    *nstar = header.nstar ;
    *nsph = header.nsph ;

    if (type=='d') *np=header.ndark ;
    if (type=='g') *np=header.nsph ;
    if (type=='s') *np=header.nstar ;
    if (type=='b') *np=header.nsph+header.nstar ;
    if (type=='a') *np=header.nsph+header.nstar+header.ndark ;

    np1= *np ;

    *x=(float *)malloc(np1*sizeof(float)) ;
    *y=(float *)malloc(np1*sizeof(float)) ;
    *z=(float *)malloc(np1*sizeof(float)) ;

    if (*x==NULL || *y==NULL || *z==NULL)  {
        fprintf(stderr,"getpos_tipsy> no memory for position arrays\n") ;
	exit(-1) ;
    }

    i=0 ;
    if (type=='g' || type=='b' || type=='a')  {
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    (*x)[i]=gp->pos[0] ;
	    (*y)[i]=gp->pos[1] ;
	    (*z)[i]=gp->pos[2] ;
	    i++ ;
	}
    }

    if (type=='d' || type=='a')  {
	for(dp=dark_particles ; dp < lastdp ; dp++) {
	    (*x)[i]=dp->pos[0] ;
	    (*y)[i]=dp->pos[1] ;
	    (*z)[i]=dp->pos[2] ;
	    i++ ;
	}
    }

    if (type=='s' || type=='b' || type=='a')  {
	for(sp=star_particles; sp < lastsp ; sp++) {
	    (*x)[i]=sp->pos[0] ;
	    (*y)[i]=sp->pos[1] ;
	    (*z)[i]=sp->pos[2] ;
	    i++ ;
	}
    }

    if (i!= *np)  {
	fprintf(stderr,"getpos_tipsy> np=%d, i=%d\n",np,i) ;
	exit(-1) ;
    }

    free(gas_particles) ;
    free(dark_particles) ;
    free(star_particles) ;
}

/* getposvel_tipsy -- get positions and velocities from a tipsy binary file
   void getposvel_tipsy(fp,time,type,*x,*y,*z,*vx,*vy,*vz,nsph,ndark,nstar,np) 
     * fp = file pointer to binary file
     * time = desired output time
     * type = d, g, s, b, a for dark, gas, star, baryon, or all positions
   variables set on return:
     * x,y,z = pointers to position arrays
     * vx,vy,vz = pointers to velocity arrays
     * nsph = number of sph particles in tipsy file
     * ndark = number of dark particles in tipsy file
     * nstar = number of sph particles in tipsy file
     * np = number of objects with returned positions
*/

void getposvel_tipsy(FILE *fp, double time, char type, 
                 float **x, float **y, float **z,
		 float **vx, float **vy, float **vz,
                 int *nsph, int *ndark, int *nstar, int *np) 
{
    int i ;
    struct gas_particle *gas_particles, *gp, *lastgp;
    struct dark_particle *dark_particles, *dp, *lastdp;
    struct star_particle *star_particles, *sp, *lastsp;
    struct dump header;
    double currtime = 0.0;
    long currpos = 0L ;
    long lastpos = 0L ;
    int np1 ;
    char *calloc() ;

    if ((float)currtime > (float)time){
	fseek(fp,0L,0);
	currtime=0.0;
	currpos=0;
    }

/* find file position for this output time */
    while(1)  {			
       if (fread((char *)&header,sizeof(header),1,fp) != 1) {
	   fprintf(stderr,"getpos_tipsy> time too large, using %f\n",
		          (float)currtime) ;
	    break ;
        }
	currtime = header.time ;
	currpos = ftell(fp) - sizeof(header);
	if ( (float)header.time >= (float)time ) 
	    break ;
	fseek(fp,
	      sizeof(gas_particles[0])*header.nsph +
	      sizeof(dark_particles[0])*header.ndark +
	      sizeof(star_particles[0])*header.nstar,
	      1) ;
    }

/* read the particle data */
    fseek(fp,currpos,0) ;
    fread((char *)&header,sizeof(header),1,fp) ;
    if(header.nsph != 0) {
	gas_particles = (struct gas_particle *)
			    malloc(header.nsph*sizeof(*gas_particles));
	if(gas_particles == NULL) {
	    fprintf(stderr,"getpos_tipsy> no memory for gas particles\n") ;
	    exit(-1) ;
	}
    }
    if(header.ndark != 0) {
	dark_particles = (struct dark_particle *)
			    malloc(header.ndark*sizeof(*dark_particles));
	if(dark_particles == NULL) {
	    fprintf(stderr,"getpos_tipsy> no memory for dark particles\n") ;
	    exit(-1) ;
	}
    }
    if(header.nstar != 0) {
	star_particles = (struct star_particle *)
			    malloc(header.nstar*sizeof(*star_particles));
	if(star_particles == NULL) {
	    fprintf(stderr,"getpos_tipsy> no memory for star particles\n") ;
	    exit(-1) ;
	}
    }

    fread((char *)gas_particles,sizeof(struct gas_particle),
		     header.nsph,fp) ;
    fread((char *)dark_particles,sizeof(struct dark_particle),
		     header.ndark,fp) ;
    fread((char *)star_particles,sizeof(struct star_particle),
		     header.nstar,fp) ;
    currpos = lastpos ;
    fseek(fp,currpos,0) ;
    currtime = header.time ;
    if ((float)time != (float)currtime){
	fprintf(stderr,"getpos_tipsy> used time %f\n",(float)currtime);
    }

    lastgp = gas_particles + header.nsph ;
    lastdp = dark_particles + header.ndark ;
    lastsp = star_particles + header.nstar ;

    *ndark = header.ndark ;
    *nstar = header.nstar ;
    *nsph = header.nsph ;

    if (type=='d') *np=header.ndark ;
    if (type=='g') *np=header.nsph ;
    if (type=='s') *np=header.nstar ;
    if (type=='b') *np=header.nsph+header.nstar ;
    if (type=='a') *np=header.nsph+header.nstar+header.ndark ;

    np1= *np ;

    *x=(float *)malloc(np1*sizeof(float)) ;
    *y=(float *)malloc(np1*sizeof(float)) ;
    *z=(float *)malloc(np1*sizeof(float)) ;
    if (*x==NULL || *y==NULL || *z==NULL)  {
        fprintf(stderr,"getpos_tipsy> no memory for position arrays\n") ;
	exit(-1) ;
    }

    *vx=(float *)malloc(np1*sizeof(float)) ;
    *vy=(float *)malloc(np1*sizeof(float)) ;
    *vz=(float *)malloc(np1*sizeof(float)) ;
    if (*vx==NULL || *vy==NULL || *vz==NULL)  {
        fprintf(stderr,"getpos_tipsy> no memory for velocity arrays\n") ;
	exit(-1) ;
    }

    i=0 ;
    if (type=='g' || type=='b' || type=='a')  {
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    (*x)[i]=gp->pos[0] ;
	    (*y)[i]=gp->pos[1] ;
	    (*z)[i]=gp->pos[2] ;
	    (*vx)[i]=gp->vel[0] ;
	    (*vy)[i]=gp->vel[1] ;
	    (*vz)[i]=gp->vel[2] ;
	    i++ ;
	}
    }

    if (type=='d' || type=='a')  {
	for(dp=dark_particles ; dp < lastdp ; dp++) {
	    (*x)[i]=dp->pos[0] ;
	    (*y)[i]=dp->pos[1] ;
	    (*z)[i]=dp->pos[2] ;
	    (*vx)[i]=dp->vel[0] ;
	    (*vy)[i]=dp->vel[1] ;
	    (*vz)[i]=dp->vel[2] ;
	    i++ ;
	}
    }

    if (type=='s' || type=='b' || type=='a')  {
	for(sp=star_particles; sp < lastsp ; sp++) {
	    (*x)[i]=sp->pos[0] ;
	    (*y)[i]=sp->pos[1] ;
	    (*z)[i]=sp->pos[2] ;
	    (*vx)[i]=sp->vel[0] ;
	    (*vy)[i]=sp->vel[1] ;
	    (*vz)[i]=sp->vel[2] ;
	    i++ ;
	}
    }

    if (i!= *np)  {
	fprintf(stderr,"getposvel_tipsy> np=%d, i=%d\n",np,i) ;
	exit(-1) ;
    }

    free(gas_particles) ;
    free(dark_particles) ;
    free(star_particles) ;
}
