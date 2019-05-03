/************************************************************************
"StackEnhanceV1"
Program to enhance vascular network image stacks
Assumes all images in stack have same dimensions
Uses percentilefiltering and vesselness filtering
w = width of images
h = height of images
d = depth of stack
GPU version
TWS, December 2016 - April 2019
**************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"
#include "tiffio.h"
#include <cutil_inline.h>
#include <time.h>
#include <Windows.h>	//needed for System

void writefile(int*** image, int w, int h, int kimage, const char name[]);
void writeslice(int*** image, int w, int jimage, int d, const char name[]);
void percentilefilter2d(float lowpercent, float highpercent, float minfactor, float maxfactor, int iwindow, int jwindow);
void makepsf2d(float **psf2d, float **psf2dx, float sigmax, float sigmay, int dpsf2d, float theta);
void makepsf3d(float ***psf3d, float ***psf3dx, float sigmax, float sigmay, int dpsf3d, float *dirvector);
void convolute3d(float ***input, float ***psf, float ***output, int w, int h, int d, int dpsf3d);
void convolute2d(float ***input, float **psf2d, float ***output, int w, int h, int dpsf2d);

void percentilefilter2dGPUc(float lowpercent, float highpercent, float minfactor, float maxfactor, int iwindow, int jwindow);
void percentilefilter3dGPUc(float lowpercent, float highpercent, float minfactor, float maxfactor, int iwindow, int jwindow, int kwindow);
void percentilefilter3dSamplingGPUc(float lowpercent, float highpercent, float minfactor, float maxfactor, int **samples, int nsamples);
void vesselness3dGPUc(float sigmax3d, float sigmay3d, int dpsf, float fac);

int w,h,d,nvoxels;
float **spherepoints28;
int ***image;
float ***temp,***temp1,***input,***output;
float **psf2d,**psf2dx,***psf3d,***psf3dx;
float sigmax,sigmay,sigmaz;
float pi1 = atan(1.)*4.;

//Needed for GPU version
int *h_image,*d_image,*d_image1,*d_image2,*d_npixint;
float *h_input,*h_output;
float *d_input,*d_output;
float *d_maxerr,*d_errorx,*d_errory,*d_errorz;
float *h_psf2d,*d_psf2d;
float *h_psf3d,*d_psf3d,*h_psf3dx,*d_psf3dx;

int main()
{
	//**************************
	int useGPU = 3;	// use GPU
	if(useGPU) cudaSetDevice( useGPU-1 );
	//**************************
	int run_percentilefilter2d = 0;	//2D percentiling
	float lowpercent2d = 50.;
	float highpercent2d = 99.;
	float minfactor2d = 1.;
	float maxfactor2d = 8.;
	int iwindow2d = 51, jwindow2d = 51;
	if(iwindow2d%2 == 0) iwindow2d++;	//window must have odd dimensions
	if(jwindow2d%2 == 0) jwindow2d++;
	//**************************
	int run_percentilefilter3d = 0;	//3D percentiling 
	float lowpercent3d = 50.;
	float highpercent3d = 99.;
	float minfactor3d = 1.;
	float maxfactor3d = 8.;
	int iwindow3d = 51, jwindow3d = 51, kwindow3d = 1;	// (2d for testing purposes)
	if(iwindow3d%2 == 0) iwindow3d++;	//window must have odd dimensions
	if(jwindow3d%2 == 0) jwindow3d++;
	if(kwindow3d%2 == 0) kwindow3d++;
	//**************************
	int run_percentilefilter3dsampling = 1;	//3D percentiling with predefined sample scheme
	float lowpercent3dsamp = 50.;
	float highpercent3dsamp = 99.;
	float minfactor3dsamp = 1.;
	float maxfactor3dsamp = 8.;
	int discradius = 50;	//sampling region is a disk plus a sphere
	int sphereradius = 10;
	//**************************
	int run_fillin3d = 1;	//2D "fillin" filter
	int fillinradius = 25;
	//**************************
	int run_vesselness2d = 0;	//2D "vesselness" filter
	float sigmax2d = 3.;
	float sigmay2d = 0.5;	//spread in pixels
	int dpsf2d = 13;
	if(dpsf2d%2 == 0) dpsf2d++;	//domain must have odd dimensions
	//**************************
	int run_vesselness3d = 1;	//3D "vesselness" filter
	float sigmax3d = 3.0;
	float sigmay3d = 0.5; //0.5;	//spread in voxels
	int dpsf3d = 13;
	if(dpsf3d%2 == 0) dpsf3d++;	//domain must have odd dimensions
	//**************************

	int i,j,k,ii,jj,kk;
	char fname[120];
	TIFF* tif;
	FILE *ifp;
	uint32 ww,hh,*raster;
	size_t npixels;
	clock_t start, end;

	//Create a NewStack folder if it does not already exist
	DWORD ftyp = GetFileAttributesA("NewStack\\");
	if (ftyp != FILE_ATTRIBUTE_DIRECTORY) system("mkdir NewStack");

	psf2d = matrix(1,dpsf2d,1,dpsf2d);
	psf2dx = matrix(1,dpsf2d,1,dpsf2d);
	psf3d = f3tensor(1,dpsf3d,1,dpsf3d,1,dpsf3d);
	psf3dx = f3tensor(1,dpsf3d,1,dpsf3d,1,dpsf3d);
	//*********************************
	//read an image stack
	d = 331;
	//note: can't easily move this read procedure to a separate function
	//because values of w and h are needed to define image array 
	for(k=1; k<=d; k++){		
		sprintf(fname,"E:\\BrainBoas2015\\stackImages\\im%03i.tif",k); // d = 331
		//sprintf(fname,"E:\\BrainBoas2015\\StackEnhance2019\\vesselness3d%03i.tif",k); // d = 331
		tif = TIFFOpen(fname, "r");
		if(k == 1){
			TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &ww);
			TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &hh);
			raster = (uint32*) _TIFFmalloc(ww*hh*sizeof(uint32));
			w = ww;
			h = hh;
			npixels = w*h;
			nvoxels = w*h*d;
			image = i3tensor(1,w,1,h,1,d);	//image in integer values			
		}
		TIFFReadRGBAImage(tif, ww, hh, raster, 0);
		TIFFClose(tif);
		//use middle 8 bits to get intensity - works for grayscale and for green images
		//need to reverse sequence of j values (vertical flip) at this point
		//for(j=1; j<=h; j++) for(i=1; i<=w; i++) image[i][j][k] = raster[(h-j)*ww + (i-1)]%256;//works for grayscale
		for(j=1; j<=h; j++) for(i=1; i<=w; i++) image[i][j][k] = (raster[(hh-j)*ww + (i-1)]/256)%256;
	}
	_TIFFfree(raster);
	//********************************
	temp = f3tensor(1,w,1,h,1,d);		//images using float values
	temp1 = f3tensor(1,w,1,h,1,d);
	input = f3tensor(1,w,1,h,1,d);
	output = f3tensor(1,w,1,h,1,d);

	//**************************
	//use percentile filtering and rescaling to remove background variations
	if(run_percentilefilter2d){
		printf("Percentile filter low %g%% high %g%%, minfactor %g, maxfactor %g, window %ix%i\n", lowpercent2d, highpercent2d, minfactor2d, maxfactor2d, iwindow2d, jwindow2d);
		if(useGPU) percentilefilter2dGPUc(lowpercent2d,highpercent2d,minfactor2d,maxfactor2d,iwindow2d,jwindow2d);
		else percentilefilter2d(lowpercent2d,highpercent2d,minfactor2d,maxfactor2d,iwindow2d,jwindow2d);		
		for(k=1; k<=d; k++) writefile(image,w,h,k,"NewStack\\percentilefiltered2d");
	}
	
	if (run_percentilefilter3d){
		printf("Percentile filter 3d low %g%% high %g%%, minfactor %g, maxfactor %g, window %ix%ix%i\n", lowpercent3d, highpercent3d, minfactor3d, maxfactor3d, iwindow3d, jwindow3d, kwindow3d);
		percentilefilter3dGPUc(lowpercent3d,highpercent3d,minfactor3d,maxfactor3d,iwindow3d,jwindow3d,kwindow3d);
		for(k=1; k<=d; k++) writefile(image,w,h,k,"NewStack\\percentilefiltered3d");
		printf("\n");
	}

	if (run_percentilefilter3dsampling){	//sampling region is a disk plus a sphere
		int nsamplesmax=pi1*(SQR(discradius) + 4./3.*SQR(sphereradius)*sphereradius);
		int **samples=imatrix(1,nsamplesmax,1,3);

		int nsamples=0;
		for (ii=-sphereradius;ii<=sphereradius;ii++) for (jj=-sphereradius;jj<=sphereradius;jj++) for (kk=-sphereradius;kk<=sphereradius;kk++) {
			if(ii*ii + jj*jj + kk*kk <= SQR(sphereradius)){
				nsamples++;
				samples[nsamples][1]=ii;
				samples[nsamples][2]=jj;
				samples[nsamples][3]=kk;
			}
		}
		for (ii=-discradius;ii<=discradius;ii++) for (jj=-discradius;jj<=discradius;jj++) {	//sample only every second point
			if((ii + jj)%2 == 0 && ii*ii + jj*jj <= SQR(discradius) && ii*ii + jj*jj > SQR(sphereradius)){
				nsamples++;
				samples[nsamples][1]=ii;
				samples[nsamples][2]=jj;
				samples[nsamples][3]=0;
			}
		}
		printf("Percentile filter 3d low %g%% high %g%%, minfactor %g, maxfactor %g\n", lowpercent3dsamp, highpercent3dsamp, minfactor3dsamp, maxfactor3dsamp);
		printf("Sampling region is a disk with radius %i and a sphere with radius %i\n", discradius,sphereradius);
		percentilefilter3dSamplingGPUc(lowpercent3dsamp,highpercent3dsamp,minfactor3dsamp,maxfactor3dsamp,samples,nsamples);
		for(k=1; k<=d; k++) writefile(image,w,h,k,"NewStack\\percentilefiltered3dsampling");
		printf("\n");
	}
	//**************************
	//2D "vesselness" filter - no GPU version yet
	if(run_vesselness2d){
		float theta = 0.;
		int itheta;
		float fac = 0.02;
		for(j=1; j<=h; j++) for(i=1; i<=w; i++){
			input[i][j][1] = 1./255.*image[i][j][1];
			output[i][j][1] = 0.;
		}
		for(itheta=0; itheta<10; itheta++){		//create rotated 2D point spread function
			theta = itheta*pi1/10.;
			if(itheta == 0) printf("2D Gaussian point spread function, s.d. = (%g, %g), diameter %i\n",sigmax2d,sigmay2d,dpsf2d,theta);
			printf("%i ",itheta); 
			makepsf2d(psf2d, psf2dx, sigmax2d, sigmay2d, dpsf2d, theta);
			convolute2d(input, psf2d, temp, w, h, dpsf2d);
			convolute2d(input, psf2dx, temp1, w, h, dpsf2d);	//x-derivative of psf
			for(j=1; j<=h; j++) for(i=1; i<=w; i++){
				temp[i][j][1] -= fac*SQR(temp1[i][j][1]);	//avoid directions with steep overall gradient - gives sharper image
				output[i][j][1] = FMAX(temp[i][j][1],output[i][j][1]);
			}
		}
		printf("\n");
		for(j=1; j<=h; j++) for(i=1; i<=w; i++) image[i][j][1] = output[i][j][1]*255.;
		writefile(image,w,h,1,"NewStack\\vesselness2d");
	}	
	//**************************
	//3D "vesselness" filter, including GPU version
	if(run_vesselness3d){
		start = clock();
		printf("3D Gaussian point spread function, s.d. = (%g, %g), diameter %i\n",sigmax3d,sigmay3d,dpsf3d);
		int iangle;
		float fac = 1;//0.02;
		float *dirvector;
		spherepoints28 = matrix(1,28,1,3);
		dirvector = vector(1,3);
		//create rotated 3D point spread
		ifp = fopen("SpherePoints28.txt","r");	//read components of direction vectors to scan a sphere
		if(ifp) {
			for(iangle=1; iangle<=28; iangle++)	for(ii=1; ii<=3; ii++) fscanf(ifp,"%f", &spherepoints28[iangle][ii]);
			fclose(ifp);
		}
		else printf("*** Error: missing SpherePoints file ***\n");

		for(i=1; i<=w; i++) for(j=1; j<=h; j++) for(k=1; k<=d; k++){
			input[i][j][k] = 1./255.*image[i][j][k];
			output[i][j][k] = 0.;
		}

		if(useGPU) vesselness3dGPUc(sigmax3d, sigmay3d, dpsf3d, fac);
		else{
			for(iangle=1; iangle<=28; iangle++){
				printf("%i ",iangle);
				for(ii=1; ii<=3; ii++) dirvector[ii] = spherepoints28[iangle][ii];
				makepsf3d(psf3d, psf3dx, sigmax3d, sigmay3d, dpsf3d, dirvector);//make kernels for filter with given orientation
				convolute3d(input, psf3d, temp, w, h, d, dpsf3d);
				convolute3d(input, psf3dx, temp1, w, h, d, dpsf3d);	//derivative of kernel along its length
				for(i=1; i<=w; i++) for(j=1; j<=h; j++) for(k=1; k<=d; k++){
					temp[i][j][k] -= fac*SQR(temp1[i][j][k]);	//avoid directions with steep gradient along kernel - gives sharper image
					output[i][j][k] = FMAX(temp[i][j][k],output[i][j][k]);
				}
			}
		}
		printf("\n");
		for(j=1; j<=h; j++) for(i=1; i<=w; i++) for(k=1; k<=d; k++) image[i][j][k] = output[i][j][k]*255.;
		for(k=1; k<=d; k++) writefile(image,w,h,k,"NewStack\\vesselness3d");
		end = clock();
		printf("3D vesselness filtering took %0.4lf seconds\n", (double) (end-start)/CLOCKS_PER_SEC);
	}
	//**************************
	//3D "fillin" filter
	if(run_fillin3d){
		start = clock();
		printf("3D fill-in filter, radius = %i\n",fillinradius);
		float max1, max2, max3, max4, max5, max6,summax2,mean,coeffv;
		for(i=1; i<=w; i++) for(j=1; j<=h; j++) for(k=1; k<=d; k++){
			input[i][j][k] = image[i][j][k];
			output[i][j][k] = input[i][j][k];
		}
		for(i=1; i<=w; i++) for(j=1; j<=h; j++) for(k=1; k<=d; k++) if(input[i][j][k] > 4.) {	//don't do this for black regions		
			max1 = input[i][j][k];
			max2 = input[i][j][k];
			max3 = input[i][j][k];
			max4 = input[i][j][k];
			max5 = input[i][j][k];
			max6 = input[i][j][k];
			for(ii=1; ii<=fillinradius; ii++){
				if(i + ii <= w) max1 = FMAX(max1,input[i+ii][j][k]);
				if(i - ii >= 1) max2 = FMAX(max2,input[i-ii][j][k]);
				if(j + ii <= h) max3 = FMAX(max3,input[i][j+ii][k]);
				if(j - ii >= 1) max4 = FMAX(max4,input[i][j-ii][k]);
				if(k + ii <= d) max5 = FMAX(max3,input[i][j][k+ii]);
				if(k - ii >= 1) max6 = FMAX(max4,input[i][j][k-ii]);
			}
			mean = (max1 + max2 + max3 + max4 + max5 + max6)/6.;
			if(mean > input[i][j][k]) {	
				summax2 = SQR(max1) + SQR(max2) + SQR(max3) + SQR(max4) + SQR(max5) + SQR(max6);
				coeffv = sqrt(summax2/SQR(mean)/6 - 1.);
				//the more consistent the high values are on the six test lines, the greater the intensity boost
				output[i][j][k] = input[i][j][k] + (mean - input[i][j][k])*exp(-coeffv);
			}
		}
		for(j=1; j<=h; j++) for(i=1; i<=w; i++) for(k=1; k<=d; k++)	image[i][j][k] = output[i][j][k];
		for(k=1; k<=d; k++)	writefile(image,w,h,k,"NewStack\\fillin3d");

		end = clock();
		printf("Fill-in filtering took %0.4lf seconds\n", (double) (end-start)/CLOCKS_PER_SEC);
	}
	//**************************
	//3D "vesselness" filter - 2nd run
	if(run_vesselness3d){
		start = clock();
		printf("3D Gaussian point spread function, s.d. = (%g, %g), diameter %i\n",sigmax3d,sigmay3d,dpsf3d);
		float fac = 1.;
		for(i=1; i<=w; i++) for(j=1; j<=h; j++) for(k=1; k<=d; k++){
			input[i][j][k] = 1./255.*image[i][j][k];
			output[i][j][k] = 0.;
		}
		if(useGPU) vesselness3dGPUc(sigmax3d, sigmay3d, dpsf3d, fac);

		printf("\n");
		for(j=1; j<=h; j++) for(i=1; i<=w; i++) for(k=1; k<=d; k++) image[i][j][k] = output[i][j][k]*255.;
		for(k=1; k<=d; k++) writefile(image,w,h,k,"NewStack\\vesselness3d");
		end = clock();
		printf("3D vesselness filtering took %0.4lf seconds\n", (double) (end-start)/CLOCKS_PER_SEC);
	}
	//**************************	
	free_matrix(psf2d,1,dpsf2d,1,dpsf2d);
	free_matrix(psf2dx,1,dpsf2d,1,dpsf2d);
	free_f3tensor(psf3d,1,dpsf3d,1,dpsf3d,1,dpsf3d);
	free_f3tensor(psf3dx,1,dpsf3d,1,dpsf3d,1,dpsf3d);
	free_i3tensor(image,1,w,1,h,1,d);
	free_f3tensor(temp,1,w,1,h,1,d);
	free_f3tensor(temp1,1,w,1,h,1,d);
	free_f3tensor(input,1,w,1,h,1,d);
	free_f3tensor(output,1,w,1,h,1,d);
	return 0;
}
