#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define FUNC(x) ((*func)(x))
#define alpha 1.0E-5
#define dx_lim 1

float integrate_los(float (*func)(float), float a, float b, float center, float dv)
{
  int i,j,n=0;
  double sum1=0;
  double x,dx,fprev,f,xprev,slope;

  /* First, take numeric derivative near the center. (problems at mode of PDF)
   */
  slope=fabs((FUNC(center+dv)-FUNC(center-dv))/(2*dv)); 

  /* Integrate from center to [b]
   */
  dx=alpha/slope;
  if(dx>dx_lim)dx=dx_lim;
  fprev=FUNC(center);
  xprev=center;
  x=center;
  sum1=0;

  printf("%f %e %f %f %e\n",dx,slope,dv,center,FUNC(100.0));

  while(x+dx<b)
    {
      n++;
      x+=dx;
      f=FUNC(x);
      sum1+=0.5*(f+fprev)*dx;

      slope=(f-fprev)/(x-xprev);
      dx=alpha/fabs(slope);
      if(dx>dx_lim && fabs(x-center)<100)dx=dx_lim;

      fprev=f;
      xprev=x;
    }
  f=FUNC(b);
  sum1+=0.5*(f+fprev)*(b-xprev);

  /* Integrate from center to [a]
   */
  dx=alpha/(FUNC(center+dv)-FUNC(center-dv))/(2*dv);  
  dx=dx_lim;
  fprev=FUNC(center);
  xprev=center;
  x=center;

  while(x-dx>a)
    {
      n++;
      x-=dx;
      f=FUNC(x);
      sum1+=0.5*(f+fprev)*dx;

      slope=(f-fprev)/(x-xprev);
      dx=alpha/fabs(slope);
      if(dx>dx_lim && fabs(x-center)<100)dx=dx_lim;

      fprev=f;
      xprev=x;
    }
  f=FUNC(a);
  sum1+=0.5*(f+fprev)*fabs(xprev-a);
  
  printf("%d\n",n);
  return(sum1);
}
