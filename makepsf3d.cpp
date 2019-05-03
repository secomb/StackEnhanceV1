/************************************************************************
makepsf3d.cpp
Generate 3D Gaussian point spread function with variable rotation
defined by a vector from the file SpherePoints28.txt
dpsf is the diameter of the psf in pixels, should be an odd integer
TWS - March 2017
***********************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void makepsf3d(float ***psf3d, float ***psf3dx, float sigmax, float sigmay, int dpsf3d, float *dirvector)
{
	int i,j,k,i1,j1,k1,rpsf3d;
	float x,y2,sum = 0.;
	rpsf3d = (1 + dpsf3d)/2;	
	for(i=1; i<=dpsf3d; i++) for(j=1; j<=dpsf3d; j++) for(k=1; k<=dpsf3d; k++){
		psf3d[i][j][k] = 0.;
        i1 = i - rpsf3d;
        j1 = j - rpsf3d;
        k1 = k - rpsf3d;
		if(SQR(i1) + SQR(j1) + SQR(k1) <= SQR(rpsf3d)){
			//parallel component - dot product
			x = i1*dirvector[1] + j1*dirvector[2] + k1*dirvector[3];
			//perpendicular component - cross product
			y2 = SQR(j1*dirvector[3] - k1*dirvector[2])
				+ SQR(k1*dirvector[1] - i1*dirvector[3])
				+ SQR(i1*dirvector[2] - j1*dirvector[1]);
			psf3d[i][j][k] = exp(-SQR(x/sigmax)/2.-y2/SQR(sigmay)/2.);
			psf3dx[i][j][k] = x*psf3d[i][j][k];
			sum += psf3d[i][j][k];
		}
	}
	if(sum != 0.) for(i=1; i<=dpsf3d; i++) for(j=1; j<=dpsf3d; j++) for(k=1; k<=dpsf3d; k++){
		psf3d[i][j][k] = psf3d[i][j][k]/sum;
		psf3dx[i][j][k] = psf3dx[i][j][k]/sum;
	}
}
