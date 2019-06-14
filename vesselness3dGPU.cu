/***********************************************************
vesselness3dGPU.cu
GPU kernel for 3D vesselness filter
TWS May 2017
************************************************************/
#include <stdio.h>
#include <cutil_inline.h>

__global__ void vesselness3dGPUKernel(float *d_input, float *d_output, float *d_psf3d, float *d_psf3dx, int w, int h, int d, int dpsf3d, float fac)
{
	int ipp0 = threadIdx.x + threadIdx.y*blockDim.x;
	int ipp2 = blockIdx.x + blockIdx.y*gridDim.x;
	int ipp = ipp0 + ipp2*blockDim.x*blockDim.y;	//voxel number

	int i,j,k,ipsf,jpsf,kpsf,ipsfmin,jpsfmin,kpsfmin,ipsfmax,jpsfmax,kpsfmax;
	int ippinput,ipppsf;
	int dpsf3d2 = (dpsf3d+1)/2;
	float temp, temp1;

    if(ipp < w*h*d){
		temp = 0.;
		temp1 = 0.;
		i = ipp%w + 1;
		j = (ipp/w)%h + 1;
		k = ipp/(w*h) + 1;
		ipsfmin = (i <= w - dpsf3d2 + 1 ? -dpsf3d2 + 1 : i - w); 
		ipsfmax = (i >= dpsf3d2 ? dpsf3d2 - 1 : i - 1);
		jpsfmin = (j <= h - dpsf3d2 + 1 ? -dpsf3d2 + 1 : j - h); 
		jpsfmax = (j >= dpsf3d2 ? dpsf3d2 - 1 : j - 1);
		kpsfmin = (k <= d - dpsf3d2 + 1 ? -dpsf3d2 + 1 : k - d); 
		kpsfmax = (k >= dpsf3d2 ? dpsf3d2 - 1 : k - 1);
		for(kpsf=kpsfmin; kpsf<=kpsfmax; kpsf++){
			for(jpsf=jpsfmin; jpsf<=jpsfmax; jpsf++){
				for(ipsf=ipsfmin; ipsf<=ipsfmax; ipsf++){
					ippinput = (i-ipsf-1) + (j-jpsf-1)*w + (k-kpsf-1)*w*h;
					ipppsf = (dpsf3d2+ipsf-1) + (dpsf3d2+jpsf-1)*dpsf3d + (dpsf3d2+kpsf-1)*dpsf3d*dpsf3d;
					temp += d_input[ippinput]*d_psf3d[ipppsf];	//convolute input with kernel
					temp1 += d_input[ippinput]*d_psf3dx[ipppsf];	//convolute input with kernel derivative
				}
			}
		}
		temp -= fac*temp1*temp1;	//avoid directions with steep gradient - gives sharper image
		d_output[ipp] = (temp >= d_output[ipp] ? temp : d_output[ipp]);	//if the value is higher for this direction, save this value
	}
}

extern "C" void vesselness3dGPU(float *d_input, float *d_output, float *d_psf3d, float *d_psf3dx, int w, int h, int d, int dpsf3d, float fac)
{
	dim3 threadsPerBlock(16, 16);
	dim3 blocksPerGrid((w+threadsPerBlock.x-1)/threadsPerBlock.x, (h*d+threadsPerBlock.x-1)/threadsPerBlock.y);
	vesselness3dGPUKernel<<<blocksPerGrid, threadsPerBlock>>>(d_input,d_output,d_psf3d,d_psf3dx,w,h,d,dpsf3d,fac);
}