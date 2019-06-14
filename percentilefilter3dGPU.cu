/***********************************************************
percentilefilter3dGPU.cu
GPU kernel for 3D percentile filtering
Uses 2D grid of blocks
TWS December 2016 - July 2017.
************************************************************/
#include <stdio.h>
#include <cutil_inline.h>

__global__ void percentilefilter3dGPUKernel(int *d_image, int *d_image1, int *d_image2, int w, int h, int d, float lowpercent, float highpercent, int iwindow, int jwindow, int kwindow)
{
	//int ipp = (blockIdx.x*blockDim.x + threadIdx.x) + (blockIdx.y*blockDim.y + threadIdx.y)*blockDim.x*(w+blockDim.x-1)/blockDim.x; //pixel number
	
	int ipp0 = threadIdx.x + threadIdx.y*blockDim.x;
	int ipp2 = blockIdx.x + blockIdx.y*gridDim.x;
	int ipp = ipp0+ipp2*blockDim.x*blockDim.y;	
	int npixsum,intmin,intmax,pixlim,npix;
	int ii,jj,kk,iimin,iimax,jjmin,jjmax,kkmin,kkmax;
	int i,j,k,ipp1;
	int iwindow2 = (iwindow-1)/2, jwindow2 = (jwindow-1)/2, kwindow2 = (kwindow-1)/2;
	int d_npixint[256];	//added July 2017

    if(ipp < w*h*d){
		i = ipp%w + 1;
		j = (ipp/w)%h + 1;
		k = ipp/(w*h) + 1;
		for(ii=0; ii<=255; ii++) d_npixint[ii] = 0;	//create distribution of pixel intensities
		npix = 0;
		iimin = (i-iwindow2 > 1 ? i-iwindow2 : 1);
		iimax = (i+iwindow2 < w ? i+iwindow2 : w);
		jjmin = (j-jwindow2 > 1 ? j-jwindow2 : 1);
		jjmax = (j+jwindow2 < h ? j+jwindow2 : h);
		kkmin = (k-kwindow2 > 1 ? k-kwindow2 : 1);
		kkmax = (k+kwindow2 < d ? k+kwindow2 : d);
		
		for(jj=jjmin; jj<=jjmax; jj+=1) for(ii=iimin; ii<=iimax; ii+=1) for (kk=kkmin; kk<=kkmax; kk+=1){
			ipp1 = (ii-1) + (jj-1)*w + (kk-1)*w*h;
			d_npixint[d_image[ipp1]]++;
			npix++;
		}
		pixlim = npix*lowpercent/100.;		//lower cutoff
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

extern "C" void percentilefilter3dGPU(int *d_image, int *d_image1, int *d_image2, int w, int h, int d, float lowpercent, float highpercent, int iwindow, int jwindow, int kwindow)
{
	dim3 threadsPerBlock(16, 16);
	dim3 blocksPerGrid((w+threadsPerBlock.x-1)/threadsPerBlock.x, (h*d+threadsPerBlock.x-1)/threadsPerBlock.y);
	percentilefilter3dGPUKernel<<<blocksPerGrid, threadsPerBlock>>>(d_image,d_image1,d_image2,w,h,d,lowpercent,highpercent,iwindow,jwindow,kwindow);
}