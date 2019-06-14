/***********************************************************
percentilefilter2dGPU.cu
GPU kernel for 2D percentile filtering (one image at a time).
TWS December 2016 - July 2017.
************************************************************/
#include <stdio.h>
#include <cutil_inline.h>

__global__ void percentilefilter2dGPUKernel(int *d_image, int *d_image1, int *d_image2, int w, int h, float lowpercent, float highpercent, int iwindow, int jwindow)
{
	int ipp = blockDim.x * blockIdx.x + threadIdx.x;	//pixel number
	int npixsum,intmin,intmax,pixlim,npix;
	int ii,jj,iimin,iimax,jjmin,jjmax;
	int i,j,ipp1;	
	int iwindow2 = (iwindow-1)/2, jwindow2 = (jwindow-1)/2;
	int d_npixint[256];

    if(ipp < w*h){
		i = ipp%w + 1;
		j = (ipp/w)%h + 1;
		for(ii=0; ii<=255; ii++) d_npixint[ii] = 0;	//create distribution of pixel intensities
		npix = 0;
		iimin = (i-iwindow2 > 1 ? i-iwindow2 : 1);
		iimax = (i+iwindow2 < w ? i+iwindow2 : w);
		jjmin = (j-jwindow2 > 1 ? j-jwindow2 : 1);
		jjmax = (j+jwindow2 < h ? j+jwindow2 : h);
		for(jj=jjmin; jj<=jjmax; jj+=1) for(ii=iimin; ii<=iimax; ii+=1){
			ipp1 = (ii-1) + (jj-1)*w;
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
		pixlim = npix*(1. - highpercent/100.);	//upper cutoff
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

extern "C" void percentilefilter2dGPU(int *d_image, int *d_image1, int *d_image2, int w, int h, float lowpercent, float highpercent, int iwindow, int jwindow)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (w*h + threadsPerBlock - 1) / threadsPerBlock;
	percentilefilter2dGPUKernel<<<blocksPerGrid, threadsPerBlock>>>(d_image,d_image1,d_image2,w,h,lowpercent,highpercent,iwindow,jwindow);
}