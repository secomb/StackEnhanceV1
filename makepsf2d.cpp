/************************************************************************
makepsf2d.cpp
Generate 2D Gaussian point spread function with variable rotation theta
dpsf is the diameter of the psf in pixels, should be an odd integer
TWS - March 2017
***********************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void makepsf2d(float **psf2d, float **psf2dx, float sigmax, float sigmay, int dpsf, float theta)
{
	int i,j,i1,j1,rpsf;
	float x,y,sum = 0.;
	rpsf = (1 + dpsf)/2;	
	for(i=1; i<=dpsf; i++) for(j=1; j<=dpsf; j++){
		psf2d[i][j] = 0.;
        i1 = i - rpsf;
        j1 = j - rpsf;
		if(SQR(i1) + SQR(j1) <= SQR(rpsf)){
			x = i1*cos(theta) - j1*sin(theta);
			y = i1*sin(theta) + j1*cos(theta);
			psf2d[i][j] = exp(-SQR(x/sigmax)/2.-SQR(y/sigmay)/2.);
			psf2dx[i][j] = x*psf2d[i][j];
			sum += psf2d[i][j];
		}
	}
	if(sum != 0.) for(i=1; i<=dpsf; i++) for(j=1; j<=dpsf; j++) psf2d[i][j] = psf2d[i][j]/sum;
}
