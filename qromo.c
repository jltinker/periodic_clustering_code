#include <math.h>
#include <stdio.h>
#include "nrutil.h"
#define EPS 1.0e-6
#define JMAX 14
#define JMAXP (JMAX+1)
#define K 5

float qromo(float (*func)(float), float a, float b,
	float (*choose)(float(*)(float), float, float, int))
{
	void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
	void nrerror(char error_text[]);
	int j;
	float ss,dss,h[JMAXP+1],s[JMAXP+1];

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=(*choose)(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) < EPS*fabs(ss)) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=h[j]/9.0;
	}
	/*nrerror("Too many steps in routing qromo");*/
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K

#define NRANSI

void polint(float xa[], float ya[], int n, float x, float *y, float *dy)
{
	int i,m,ns=1;
	float den,dif,dift,ho,hp,w;
	float *c,*d;
	
	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);

}
#undef NRANSI

#define FUNC(x) ((*func)(x))

float midpnt(float (*func)(float), float a, float b, int n)
{
	float x,tnm,sum,del,ddel;
	static float s;
	int it,j;

	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}
#undef FUNC

#define FUNC(x) ((*funk)(1.0/(x))/((x)*(x)))

float midinf(float (*funk)(float), float aa, float bb, int n)
{
	float x,tnm,sum,del,ddel,b,a;
	static float s;
	int it,j;

	b=1.0/aa;
	a=1.0/bb;
	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		return (s=(s+(b-a)*sum/tnm)/3.0);
	}
}
#undef FUNC
