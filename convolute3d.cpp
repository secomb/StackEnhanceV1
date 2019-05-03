/************************************************************************
convolute3d.cpp
Basic 3D convolution routine for real arguments
Trims result to dimensions of input array
Assumes PSF has odd dimensions, with origin at midpoint
TWS - December 2016
***********************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void convolute3d(float ***input, float ***psf, float ***output, int w, int h, int d, int dpsf)
{
	int i,j,k,ipsf,jpsf,kpsf,ipsfmin,jpsfmin,kpsfmin,ipsfmax,jpsfmax,kpsfmax;
	int dpsf2 = (dpsf + 1)/2;

	for(k=1; k<=d; k++) for(j=1; j<=h; j++) for(i=1; i<=w; i++) output[i][j][k] = 0.;

	for(k=1; k<=d; k++){
		kpsfmin = (k <= d - dpsf2 + 1  ? -dpsf2 + 1 : k - d); 
		kpsfmax = (k >= dpsf2 ? dpsf2 - 1 : k - 1);
		for(kpsf=kpsfmin; kpsf<=kpsfmax; kpsf++){
			for(j=1; j<=h; j++){
				jpsfmin = (j <= h - dpsf2 + 1 ? -dpsf2 + 1 : j - h); 
				jpsfmax = (j >= dpsf2 ? dpsf2 - 1 : j - 1);
				for(jpsf=jpsfmin; jpsf<=jpsfmax; jpsf++){
					for(i=1; i<=w; i++){
						ipsfmin = (i <= w - dpsf2 + 1 ? -dpsf2 + 1 : i - w); 
						ipsfmax = (i >= dpsf2 ? dpsf2 - 1 : i - 1);
						for(ipsf=ipsfmin; ipsf<=ipsfmax; ipsf++){
							output[i][j][k] += input[i-ipsf][j-jpsf][k-kpsf]*psf[dpsf2+ipsf][dpsf2+jpsf][dpsf2+kpsf];
						}
					}
				}
			}
		}
	}
}
