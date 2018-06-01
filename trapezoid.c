#include <stdlib.h>
#include <stdio.h>
#include <math.h>

float ANSWER;
int MODE;

float f1(float x);
void ans(void);


main(int argc, char **argv)
{
  float dx,x,sum,x_low,x_high,prev_sum;
  float error=10,tol;
  int Nsteps,i,flag=1;

  if(argc<4)exit(0);
  x_low=atof(argv[1]);
  x_high=atof(argv[2]);
  Nsteps=atoi(argv[3]);
  if(argc<5)tol=0.001;
  else tol=atof(argv[4]);

  if(argc<6)MODE=1;
  else MODE=atoi(argv[5]);
  ans();

  dx=2*(x_high-x_low)/Nsteps;
  sum=0;

  while(error>tol)
    {
      prev_sum=sum;
      dx=dx*0.5;
      Nsteps=(x_high-x_low)/dx;
      
      sum=0.5*dx*f1(x_low);
      for(i=2,x=x_low+dx;i<=Nsteps;++i,x+=dx)
	sum+=dx*f1(x);
      sum+=0.5*(dx)*f1(x_high);

      if(flag)
	{ flag=0; continue; }

      error=(prev_sum-sum)/prev_sum;
      error=sqrt(error*error);
      printf("%d\t%f  %f\n",Nsteps,sum,error);
    }
}

float f1(float x)
{
  float y;

  switch(MODE) {
  case 1:
    y=pow(x,-1.5);
    return(y);
  case 2:
    y=1/(pow(x,1.5)*(1+pow(x,1.5)));
    return(y);
  case 3:
    y=sin(x)/x;
    return(y);
  case 4:
    y=sin(x)*sin(x)/sqrt(x);
    return(y);
  }
}

void ans()
{
  switch(MODE){
  case 1:
    ANSWER=1.105572809;
    return;
  case 2:
    ANSWER=0.2240291;
    return;
  case 3:
    ANSWER=0.616142;
    return;
  case 4:
    ANSWER=30.8399;
    return;
  }
}

      
      



