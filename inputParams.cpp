/************************************************************
input - reads input parameters for StackEnhance. TWS June2019
*************************************************************/
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cutil_inline.h>
#include "nrutil.h"

void inputParams(void)
{
	extern int d, useGPU, run_percentilefilter2d, iwindow2d, jwindow2d, run_percentilefilter3d, iwindow3d, jwindow3d, kwindow3d;
	extern int run_percentilefilter3dsampling, discradius, sphereradius;
	extern int run_fillin3d, fillinradius, run_vesselness2d, run_vesselness3d, dpsf2d, dpsf3d, run_vesselness3d_again;

	extern float lowpercent2d, highpercent2d, minfactor2d, maxfactor2d;
	extern float lowpercent3d, highpercent3d, minfactor3d, maxfactor3d;
	extern float lowpercent3dsamp, highpercent3dsamp, minfactor3dsamp, maxfactor3dsamp;
	extern float sigmax2d, sigmay2d, sigmax3d, sigmay3d;

	extern char fname0[200];

	int max = 200;
	char bb[200], *pos;
	FILE *ifp;

	ifp = fopen("StackEnhanceParams.dat", "r");
	fgets(bb, max, ifp);
	printf("%s\n", bb);
	fscanf(ifp,"%i %*[^\n]", &d);					//number of images in stack
	fgets(bb, max, ifp);
	fgets(bb, max, ifp);
	fgets(fname0, max, ifp);
	if ((pos = strchr(fname0, '\n')) != NULL) *pos = '\0';	//remove newline character from string
	fscanf(ifp,"%i%*[^\n]", &useGPU);
	if(useGPU) cudaSetDevice( useGPU-1 );
	//**************************
	fscanf(ifp,"%i%*[^\n]", &run_percentilefilter2d);	//2D percentiling
	fscanf(ifp,"%f%*[^\n]", &lowpercent2d);
	fscanf(ifp,"%f%*[^\n]", &highpercent2d);
	fscanf(ifp,"%f%*[^\n]", &minfactor2d);
	fscanf(ifp,"%f%*[^\n]", &maxfactor2d);
	fscanf(ifp,"%i %i%*[^\n]", &iwindow2d, &jwindow2d);
	if(iwindow2d%2 == 0) iwindow2d++;			//window must have odd dimensions
	if(jwindow2d%2 == 0) jwindow2d++;
	//**************************
	fscanf(ifp,"%i%*[^\n]", &run_percentilefilter3d);	//3D percentiling
	fscanf(ifp,"%f%*[^\n]", &lowpercent3d);
	fscanf(ifp,"%f%*[^\n]", &highpercent3d);
	fscanf(ifp,"%f%*[^\n]", &minfactor3d);
	fscanf(ifp,"%f%*[^\n]", &maxfactor3d);
	fscanf(ifp,"%i %i %i%*[^\n]", &iwindow3d, &jwindow3d, &kwindow3d);
	if(iwindow3d%2 == 0) iwindow3d++;			//window must have odd dimensions
	if(jwindow3d%2 == 0) jwindow3d++;
	if(kwindow3d%2 == 0) kwindow3d++;
	//**************************
	fscanf(ifp,"%i%*[^\n]", &run_percentilefilter3dsampling);	//3D percentiling with predefined sample scheme
	fscanf(ifp,"%f%*[^\n]", &lowpercent3dsamp);
	fscanf(ifp,"%f%*[^\n]", &highpercent3dsamp);
	fscanf(ifp,"%f%*[^\n]", &minfactor3dsamp);
	fscanf(ifp,"%f%*[^\n]", &maxfactor3dsamp);
	fscanf(ifp,"%i%*[^\n]", &discradius);		//sampling region is a disk plus a sphere
	fscanf(ifp,"%i%*[^\n]", &sphereradius);
	//**************************
	fscanf(ifp,"%i%*[^\n]", &run_vesselness2d);	//2D "vesselness" filter
	fscanf(ifp,"%f%*[^\n]", &sigmax2d);
	fscanf(ifp,"%f%*[^\n]", &sigmay2d);			//spread in pixels
	fscanf(ifp,"%i%*[^\n]", &dpsf2d);
	if(dpsf2d%2 == 0) dpsf2d++;					//domain must have odd dimensions
	//**************************
	fscanf(ifp,"%i%*[^\n]", &run_vesselness3d);	//3D "vesselness" filter
	fscanf(ifp,"%f%*[^\n]", &sigmax3d);
	fscanf(ifp,"%f%*[^\n]", &sigmay3d);			//spread in voxels
	fscanf(ifp,"%i%*[^\n]", &dpsf3d);
	if(dpsf3d%2 == 0) dpsf3d++;					//domain must have odd dimensions
	//**************************
	fscanf(ifp,"%i%*[^\n]", &run_fillin3d);		//3D "fillin" filter
	fscanf(ifp,"%i%*[^\n]", &fillinradius);
	//**************************
	fscanf(ifp,"%i%*[^\n]", &run_vesselness3d_again);	//3D "vesselness" filter, optional second run
	fclose(ifp);	
}

