#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N 101

void linearize(float xi[N][N], sigma[N][N], pi[N][N], int nx, int ny)
{
  int i,j,k;
  float linxi[N][N],linx[N],liny[N];

  for(i=0;i<nx;++i)
    liny[i]=linx[i]=1.0+10.0*i/(float)nx;

  for(i=0;i<nx;++i)
    for(j=0;j<nx;++j)
      {
	

}
