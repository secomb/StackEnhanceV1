/************************************************************************
percentilefilter2d.cpp
Perform percentile filter on stack using specified 2D window
Also includes option to rescale intensities from 0 to 1
TWS - December 2016
***********************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void percentilefilter2d(float lowpercent, float highpercent, float minfactor, float maxfactor, int iwindow, int jwindow)
{
	extern int ***image;
	extern int w,h,d;
	int npixsum,intmin,intmax,intdiff,pixlim,npix;
	int ii,jj,iimin,iimax,jjmin,jjmax;
	int i,j,k;
	int iwindow2 = (iwindow-1)/2,jwindow2 = (jwindow-1)/2;
	float intfactor;
	int *npixint,***image1,***image2;
	npixint = ivector(0,255);
	image1 = i3tensor(1,w,1,h,1,d);	
	image2 = i3tensor(1,w,1,h,1,d);	

	for(k=1; k<=d; k++){
		printf("%i ",k);
		for(j=1; j<=h; j++) for(i=1; i<=w; i++){	
			for(ii=0; ii<=255; ii++) npixint[ii] = 0;	//create distribution of pixel intensities
			npix = 0;
			iimin = IMAX(i-iwindow2,1);
			iimax = IMIN(i+iwindow2,w);
			jjmin = IMAX(j-jwindow2,1);
			jjmax = IMIN(j+jwindow2,h);
			for(jj=jjmin; jj<=jjmax; jj+=2) for(ii=iimin; ii<=iimax; ii+=2){
				npixint[image[ii][jj][k]]++;
				npix++;
			}
			pixlim = npix*lowpercent/100.;	//lower cutoff
			npixsum = 0;
			intmin = 0;
			do{			//find intensity with at least pixlim pixels at or below it
				npixsum += npixint[intmin];
				intmin++;
			}
			while(npixsum < pixlim && intmin <= 255);
			intmin--;
			image1[i][j][k] = intmin;	//upper cutoff
			pixlim = npix*(1. - highpercent/100.);
			npixsum = 0;
			intmax = 255;
			do{			//find intensity with at least pixlim pixels at or above it
				npixsum += npixint[intmax];
				intmax--;
			}
			while(npixsum < pixlim && intmax >= 0);
			intmax++;
			image2[i][j][k] = intmax;
		}
		for(j=1; j<=h; j++) for(i=1; i<=w; i++){
			intdiff = image2[i][j][k] - image1[i][j][k];
			if(intdiff == 0) intfactor = maxfactor;	//limit to a maximum factor of intensification
			else{
				intfactor = FMIN(maxfactor,255./intdiff);
				intfactor = FMAX(minfactor,255./intdiff);
			}
			image[i][j][k] = (image[i][j][k] - image1[i][j][k])*intfactor;
			image[i][j][k] =  IMIN(image[i][j][k],255);
			image[i][j][k] =  IMAX(image[i][j][k],0);
		}
	}
	printf("\n");
	free_ivector(npixint,0,255);
	free_i3tensor(image1,1,w,1,h,1,d);
	free_i3tensor(image2,1,w,1,h,1,d);
}
