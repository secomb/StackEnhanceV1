/************************************************************************
convolute2d.cpp
Basic 2D convolution routine for real arguments
Trims result to dimensions of input array
Assumes PSF has odd dimensions, with origin at midpoint
TWS - March 2017
***********************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void convolute2d(float ***input, float **psf2d, float ***output, int w, int h, int dpsf)
{
	int i,j,ipsf,jpsf,ipsfmin,jpsfmin,ipsfmax,jpsfmax;
	int dpsf2 = (dpsf + 1)/2;
	float psf2dsum = 0.;

	for(j=1; j<=h; j++) for(i=1; i<=w; i++) output[i][j][1] = 0.;

	for(j=1; j<=h; j++){
		jpsfmin = (j <= h - dpsf2 + 1 ? -dpsf2 + 1 : j - h); 
		jpsfmax = (j >= dpsf2 ? dpsf2 - 1 : j - 1);
		for(jpsf=jpsfmin; jpsf<=jpsfmax; jpsf++){
			for(i=1; i<=w; i++){
				ipsfmin = (i <= w - dpsf2 + 1 ? -dpsf2 + 1 : i - w); 
				ipsfmax = (i >= dpsf2 ? dpsf2 - 1 : i - 1);
				for(ipsf=ipsfmin; ipsf<=ipsfmax; ipsf++){
					output[i][j][1] += input[i-ipsf][j-jpsf][1]*psf2d[dpsf2+ipsf][dpsf2+jpsf];
					psf2dsum += psf2d[dpsf2+ipsf][dpsf2+jpsf];
				}
			}
		}
	}
}