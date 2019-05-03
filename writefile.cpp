/************************************************************************
writefile.cpp
write tif file from numerical data in one slice of a stack
TWS - November 2016
*******************************************************************
Note: for strictly grayscale output, can make file smaller:
create raster1 in uint8 and use
TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISWHITE);
TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
TIFFWriteEncodedStrip(tif, 0, &raster1[0], w*h);
***********************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"
#include "tiffio.h"

void writefile(int*** image, int w, int h, int kimage, const char name[]) 
{
	TIFF* tif;
	size_t npixels = w*h;
	//uint32 *raster;
	//raster = (uint32*) _TIFFmalloc(npixels*sizeof(uint32));
	uint8 *raster;
	raster = (uint8*) _TIFFmalloc(npixels*sizeof(uint8));
	int i,j;
	char fname[120];
	
	//recreate image - grayscale
	//for(j=1; j<=h; j++)	for(i=1; i<=w; i++) raster[(j-1)*w + (i-1)] = 0xff000000 + (1+256+256*256)*image[i][j][kimage];	
	for(j=1; j<=h; j++)	for(i=1; i<=w; i++) raster[(j-1)*w + (i-1)] = image[i][j][kimage];	
	//modify image - keep only red channel in this example
	//for(j=1; j<=h; j++)	for(i=1; i<=w; i++) raster[(j-1)*w + (i-1)] = 0xff000000 + image[i][j][k];
	sprintf(fname, "%s%03i.tif",name,kimage);
	//sprintf(fname, "%s.tif",name);
	tif = TIFFOpen(fname, "w");
	TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, (uint32) w);
	TIFFSetField(tif, TIFFTAG_IMAGELENGTH, (uint32) h);
	TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
	TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
	//TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
	//TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 4);
	TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, h); 
	TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	//TIFFWriteEncodedStrip(tif, 0, &raster[0], 4*w*h);
	TIFFWriteEncodedStrip(tif, 0, &raster[0], w*h);
	TIFFClose(tif);
	_TIFFfree(raster);	
}
