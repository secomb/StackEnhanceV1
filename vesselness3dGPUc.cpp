/**************************************************************************
vesselness3dGPUc.cpp
GPU version of 3D vesselness filter
TWS, May 2017
**************************************************************************/
#include <shrUtils.h>
#include <cutil_inline.h>
#include "nrutil.h"

extern "C" void vesselness3dGPU(float *d_input, float *d_output, float *d_psf3d, float *d_psf3dx, int w, int h, int d, int dpsf3d, float fac);
void makepsf3d(float ***psf3d, float ***psf3dx, float sigmax, float sigmay, int dpsf3d, float *dirvector);

void vesselness3dGPUc(float sigmax3d, float sigmay3d, int dpsf3d, float fac)
{
	extern int w,h,d,nvoxels;
	extern float ***input;
	extern float ***psf3d,***psf3dx,***input,***output,**spherepoints28;
	extern float *h_psf3d,*h_psf3dx,*h_input,*h_output;
	extern float *d_psf3d,*d_psf3dx,*d_input,*d_output;

	int nvoxelspsf = dpsf3d*dpsf3d*dpsf3d;
	int i,j,k,iangle,ii;
	float *dirvector;
	dirvector = vector(1,3);

	cudaMalloc((void **) &d_input, nvoxels*sizeof(float));
	cudaMalloc((void **) &d_output, nvoxels*sizeof(float));
	cudaMalloc((void **) &d_psf3d, nvoxelspsf*sizeof(float));
	cudaMalloc((void **) &d_psf3dx, nvoxelspsf*sizeof(float));

	h_input = vector(0,nvoxels-1);
	h_output = vector(0,nvoxels-1);
	h_psf3d = vector(0,nvoxelspsf-1);
	h_psf3dx = vector(0,nvoxelspsf-1);

	//initialize input and output
	for(k=1; k<=d; k++) for(j=1; j<=h; j++) for(i=1; i<=w; i++){
		h_input[(i-1) + (j-1)*w + (k-1)*w*h] = input[i][j][k];
		h_output[(i-1) + (j-1)*w + (k-1)*w*h] = 0.;
	}
	cudaMemcpy(d_input, h_input, nvoxels*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_output, h_output, nvoxels*sizeof(float), cudaMemcpyHostToDevice);

	for(iangle=1; iangle<=28; iangle++){
		printf("%i ",iangle);
		for(ii=1; ii<=3; ii++) dirvector[ii] = spherepoints28[iangle][ii];
		makepsf3d(psf3d, psf3dx, sigmax3d, sigmay3d, dpsf3d, dirvector);	//make kernels for filter with given orientation
		for(k=1; k<=dpsf3d; k++) for(j=1; j<=dpsf3d; j++) for(i=1; i<=dpsf3d; i++){
			h_psf3d[(i-1) + (j-1)*dpsf3d + (k-1)*dpsf3d*dpsf3d] = psf3d[i][j][k];
			h_psf3dx[(i-1) + (j-1)*dpsf3d + (k-1)*dpsf3d*dpsf3d] = psf3dx[i][j][k];
		}
		cudaMemcpy(d_psf3d, h_psf3d, nvoxelspsf*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_psf3dx, h_psf3dx, nvoxelspsf*sizeof(float), cudaMemcpyHostToDevice);
		//******* call to GPU routine *******
		vesselness3dGPU(d_input,d_output,d_psf3d,d_psf3dx,w,h,d,dpsf3d,fac);
		//***********************************
	}
	printf("\n");
	cudaMemcpy(h_output, d_output, nvoxels*sizeof(float), cudaMemcpyDeviceToHost);
	for(k=1; k<=d; k++) for(j=1; j<=h; j++) for(i=1; i<=w; i++) output[i][j][k] = h_output[(i-1) + (j-1)*w + (k-1)*w*h];

	free_vector(h_input,0,nvoxels-1);
	free_vector(h_output,0,nvoxels-1);
	free_vector(h_psf3d,0,nvoxelspsf-1);
	free_vector(h_psf3dx,0,nvoxelspsf-1);
	
	cudaFree(d_input);
	cudaFree(d_output);
	cudaFree(d_psf3d);
	cudaFree(d_psf3dx);
}