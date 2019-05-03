/**************************************************************************
percentile3dSamplingGPUc.cpp
GPU version of 3D percentilefiltering
Uses predefined 3D set of sampling points
Bohan Li, 2019
**************************************************************************/
#include <shrUtils.h>
#include <cutil_inline.h>
#include <cublas.h>
#include <time.h>
#include "nrutil.h"

extern "C" void percentilefilter3dSamplingGPU(int *d_image, int *d_image1, int *d_image2, int w, int h, int d, float lowpercent, float highpercent, int *samples, int nsamples);
extern "C" void rescale3dGPU(int *d_image, int *d_image1, int *d_image2, int w, int h, int d, float minfactor, float maxfactor);

void percentilefilter3dSamplingGPUc(float lowpercent, float highpercent, float minfactor, float maxfactor, int **samples, int nsamples)
{
	extern int ***image;
	extern int w,h,d,nvoxels;
	extern int *h_image,*d_image,*d_image1,*d_image2;

	int i,j,k,isample;
	int *h_samples,*d_samples;

	clock_t start, end;
	start = clock();
	end = clock();

	h_image = ivector(0,nvoxels-1);
	h_samples = ivector(0,nsamples*3-1);

	cudaMalloc((void **) &d_image, nvoxels*sizeof(int));
	cudaMalloc((void **) &d_image1, nvoxels*sizeof(int));
	cudaMalloc((void **) &d_image2, nvoxels*sizeof(int));
	cudaMalloc((void **) &d_samples, nsamples*3*sizeof(int));

	//initialize input and output
	for(k=1; k<=d; k++) for(j=1; j<=h; j++) for(i=1; i<=w; i++)	h_image[(i-1) + (j-1)*w + (k-1)*w*h] = image[i][j][k];
	for (isample=0;isample<nsamples;isample++) for (i=1;i<=3;i++) h_samples[isample*3+i-1] = samples[isample+1][i];

	cudaMemcpy(d_image, h_image, nvoxels*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_samples, h_samples, nsamples*3*sizeof(int), cudaMemcpyHostToDevice);

	printf("Running 3d percentilefilteringsampling\n");
	percentilefilter3dSamplingGPU(d_image,d_image1,d_image2,w,h,d,lowpercent,highpercent,d_samples,nsamples);
	rescale3dGPU(d_image,d_image1,d_image2,w,h,d,minfactor,maxfactor);

	cudaMemcpy(h_image, d_image, nvoxels*sizeof(int), cudaMemcpyDeviceToHost);
	for(k=1; k<=d; k++) for(j=1; j<=h; j++) for(i=1; i<=w; i++)	image[i][j][k] = h_image[(i-1) + (j-1)*w + (k-1)*w*h];

	free_ivector(h_image,0,nvoxels-1);
	free_ivector(h_samples,0,nsamples*3-1);
		
	cudaFree(d_image);
	cudaFree(d_image1);
	cudaFree(d_image2);
	cudaFree(d_samples);

	end = clock();
	printf("Percentile filtering took %0.4lf seconds\n", (double) (end-start)/CLOCKS_PER_SEC);
}