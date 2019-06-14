/***********************************************************
percentilefilter3dSamplingGPU.cu
GPU version of 3D percentilefiltering
Uses predefined 3D set of sampling points
Uses 2D grid of blocks
Bohan Li, 2019
************************************************************/
#include <stdio.h>
#include <cutil_inline.h>
#include <math.h>

__global__ void percentilefilter3dSamplingGPUKernel(int *d_image, int *d_image1, int *d_image2, int w, int h, int d, float lowpercent, float highpercent, int *d_samples, int nsamples)
{
	int ipp0 = threadIdx.x + threadIdx.y*blockDim.x;
	int ipp2 = blockIdx.x + blockIdx.y*gridDim.x;
	int ipp = ipp0 + ipp2*blockDim.x*blockDim.y;
	int npixsum,intmin,intmax,pixlim,npix;
	int ii,jj,kk,i,j,k,ipp1,isample;
	int d_npixint[256];

	if(ipp < w*h*d){
		i = ipp%w + 1;
		j = (ipp/w)%h + 1;
		k = ipp/(w*h) + 1;
		for(ii=0; ii<=255; ii++) d_npixint[ii] = 0;	//create distribution of pixel intensities
		npix = 0;
		for(isample=0; isample<nsamples; isample++){
			ii = d_samples[isample*3] + i;
			jj = d_samples[isample*3+1] + j;
			kk = d_samples[isample*3+2] + k;
			if (ii<1 || ii>w || jj<1 || jj>h ||kk<1 || kk>d) continue;
			ipp1 = (ii-1) + (jj-1)*w + (kk-1)*w*h;
			d_npixint[d_image[ipp1]]++;
			npix++;
		}
		pixlim = npix*lowpercent/100.;	//lower cutoff
		npixsum = 0;
		intmin = 0;
		do{			//find intensity with at least pixlim pixels at or below it
			npixsum += d_npixint[intmin];
			intmin++;
		}
		while(npixsum < pixlim && intmin <= 255);
		intmin--;
		d_image1[ipp] = intmin;
		pixlim = npix*(1. - highpercent/100.);		//upper cutoff
		npixsum = 0;
		intmax = 255;
		do{			//find intensity with at least pixlim pixels at or above it
			npixsum += d_npixint[intmax];
			intmax--;
		}
		while(npixsum < pixlim && intmax >= 0);
		intmax++;
		d_image2[ipp] = intmax;
	}
}

extern "C" void percentilefilter3dSamplingGPU(int *d_image, int *d_image1, int *d_image2, int w, int h, int d, float lowpercent, float highpercent, int *d_samples, int nsamples)
{
	dim3 threadsPerBlock(16, 16);
	dim3 blocksPerGrid((w+threadsPerBlock.x-1)/threadsPerBlock.x, (h*d+threadsPerBlock.x-1)/threadsPerBlock.y);
	percentilefilter3dSamplingGPUKernel<<<blocksPerGrid, threadsPerBlock>>>(d_image,d_image1,d_image2,w,h,d,lowpercent,highpercent,d_samples,nsamples);
}