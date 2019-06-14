/***********************************************************
rescale2dGPU.cu
GPU kernel for 2D elementwise rescaling 
TWS December 2016
************************************************************/
#include <stdio.h>
#include <cutil_inline.h>

__global__ void rescale2dGPUKernel(int *d_image, int *d_image1, int *d_image2, int w, int h, float minfactor, float maxfactor)
{
	int ipp = blockDim.x * blockIdx.x + threadIdx.x;	//pixel number
	int intdiff;
	float intfactor;
    if(ipp < w*h){
		intdiff = d_image2[ipp] - d_image1[ipp];
		if(intdiff == 0) intfactor = maxfactor;	//limit to a maximum factor of intensification
		else{
			intfactor = 255./intdiff;
			intfactor = (intfactor < maxfactor ? intfactor : maxfactor);
			intfactor = (intfactor > minfactor ? intfactor : minfactor);
		}
		d_image[ipp] = (d_image[ipp] - d_image1[ipp])*intfactor;
		d_image[ipp] = (d_image[ipp] < 255 ? d_image[ipp] : 255);
		d_image[ipp] = (d_image[ipp] > 0 ? d_image[ipp] : 0);
	}
}

extern "C" void rescale2dGPU(int *d_image, int *d_image1, int *d_image2, int w, int h, float minfactor, float maxfactor)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (w*h + threadsPerBlock - 1) / threadsPerBlock;
	rescale2dGPUKernel<<<blocksPerGrid, threadsPerBlock>>>(d_image,d_image1,d_image2,w,h,minfactor,maxfactor);
}