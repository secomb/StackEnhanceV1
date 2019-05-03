/************************************************************************
percentilefilter2dGPU.cpp
Perform 2D percentile filter (one image at a time)
TWS - December 2016
***********************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <cutil_inline.h>
#include "nrutil.h"

extern "C" void percentilefilter2dGPU(int *d_image, int *d_image1, int *d_image2, int w, int h, float lowpercent, float highpercent, int iwindow, int jwindow);
extern "C" void rescale2dGPU(int *d_image, int *d_image1, int *d_image2, int w, int h, float minfactor, float maxfactor);

void percentilefilter2dGPUc(float lowpercent, float highpercent, float minfactor, float maxfactor, int iwindow, int jwindow)
{
	extern int ***image;
	extern int w,h,d;
	extern int *h_image,*d_image,*d_image1,*d_image2;

	int i,j,k;

	h_image = ivector(0,w*h-1);
	cudaMalloc((void **) &d_image, w*h*sizeof(int));
	cudaMalloc((void **) &d_image1, w*h*sizeof(int));
	cudaMalloc((void **) &d_image2, w*h*sizeof(int));

	for(k=1; k<=d; k++){
		printf("%i ",k);
		for(j=1; j<=h; j++) for(i=1; i<=w; i++) h_image[(i-1) + (j-1)*w] = image[i][j][k];
		cudaMemcpy(d_image, h_image, w*h*sizeof(int), cudaMemcpyHostToDevice);

		percentilefilter2dGPU(d_image,d_image1,d_image2,w,h,lowpercent,highpercent,iwindow,jwindow);
		rescale2dGPU(d_image,d_image1,d_image2,w,h,minfactor,maxfactor);
	
		cudaMemcpy(h_image, d_image, w*h*sizeof(int), cudaMemcpyDeviceToHost);
		for(j=1; j<=h; j++) for(i=1; i<=w; i++)	image[i][j][k] = h_image[(i-1) + (j-1)*w];
	}
	printf("\n");
	free_ivector(h_image,0,w*h-1);
		
	cudaFree(d_image);
	cudaFree(d_image1);
	cudaFree(d_image2);
}
