/***********************************************************
rescale3dGPU.cu
GPU kernel for 3D elementwise rescaling
Uses 2D grid of blocks 
TWS December 2016
************************************************************/
#include <stdio.h>
#include <cutil_inline.h>

__global__ void rescale3dGPUKernel(int *d_image, int *d_image1, int *d_image2, int w, int h, int d, float minfactor, float maxfactor)
{
	int ipp0 = threadIdx.x + threadIdx.y*blockDim.x;
	int ipp2 = blockIdx.x + blockIdx.y*gridDim.x;
	int ipp = ipp0 + ipp2*blockDim.x*blockDim.y;
	int intdiff;
	float intfactor;

	if(ipp < w*h*d){
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

extern "C" void rescale3dGPU(int *d_image, int *d_image1, int *d_image2, int w, int h, int d, float minfactor, float maxfactor)
{
	dim3 threadsPerBlock(16, 16);
	dim3 blocksPerGrid((w+threadsPerBlock.x-1)/threadsPerBlock.x, (h*d+threadsPerBlock.x-1)/threadsPerBlock.y);
	rescale3dGPUKernel<<<blocksPerGrid, threadsPerBlock>>>(d_image,d_image1,d_image2,w,h,d,minfactor,maxfactor);
}