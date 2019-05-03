/**************************************************************************
percentile3dGPUc.cpp
GPU version of 3D percentilefiltering
**************************************************************************/
#include <shrUtils.h>
#include <cutil_inline.h>
#include <cublas.h>
#include <time.h>
#include "nrutil.h"

extern "C" void percentilefilter3dGPU(int *d_image, int *d_image1, int *d_image2, int w, int h, int d, float lowpercent, float highpercent, int iwindow, int jwindow, int kwindow);
extern "C" void rescale3dGPU(int *d_image, int *d_image1, int *d_image2, int w, int h, int d, float minfactor, float maxfactor);

void percentilefilter3dGPUc(float lowpercent, float highpercent, float minfactor, float maxfactor, int iwindow, int jwindow, int kwindow)
{
	extern int ***image;
	extern int w,h,d,nvoxels;
	extern int *h_image,*d_image,*d_image1,*d_image2;

	int i,j,k;

	clock_t start, end;
	start = clock();
	end = clock();

	h_image = ivector(0,nvoxels-1);

	cudaMalloc((void **) &d_image, nvoxels*sizeof(int));
	cudaMalloc((void **) &d_image1, nvoxels*sizeof(int));
	cudaMalloc((void **) &d_image2, nvoxels*sizeof(int));

	//initialize input and output
	for(k=1; k<=d; k++) for(j=1; j<=h; j++) for(i=1; i<=w; i++)	h_image[(i-1) + (j-1)*w + (k-1)*w*h] = image[i][j][k];

	cudaMemcpy(d_image, h_image, nvoxels*sizeof(int), cudaMemcpyHostToDevice);

	printf("Running 3d percentilefiltering\n");
	percentilefilter3dGPU(d_image,d_image1,d_image2,w,h,d,lowpercent,highpercent,iwindow,jwindow,kwindow);
	rescale3dGPU(d_image,d_image1,d_image2,w,h,d,minfactor,maxfactor);

	cudaMemcpy(h_image, d_image, nvoxels*sizeof(int), cudaMemcpyDeviceToHost);

	for(k=1; k<=d; k++) for(j=1; j<=h; j++) for(i=1; i<=w; i++) image[i][j][k] = h_image[(i-1) + (j-1)*w + (k-1)*w*h];

	free_ivector(h_image,0,nvoxels-1);
		
	cudaFree(d_image);
	cudaFree(d_image1);
	cudaFree(d_image2);

	end = clock();
	printf("Percentile filtering took %0.4lf seconds\n", (double) (end-start)/CLOCKS_PER_SEC);
}