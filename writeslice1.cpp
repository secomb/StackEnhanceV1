/************************************************************************
writeslice.cpp
write tif file from numerical data in one vertical slice of a stack
TWS - December 2016
***********************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"
#include "tiffio.h"

void writeslice(int*** image, int w, int jimage, int d, const char name[]) 
{
	TIFF* tif;
	size_t npixels = w*d;
	uint8 *raster;
	raster = (uint8*) _TIFFmalloc(npixels*sizeof(uint8));
	int i,k;
	char fname[80];	

	for(k=1; k<=d; k++)	for(i=1; i<=w; i++) raster[(k-1)*w + (i-1)] = image[i][jimage][k];	
	sprintf(fname, "NewStack\\%s%03i.tif",name,jimage);
	tif = TIFFOpen(fname, "w");
	TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, (uint32) w);
	TIFFSetField(tif, TIFFTAG_IMAGELENGTH, (uint32) d);
	TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
	TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, d); 
	TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFWriteEncodedStrip(tif, 0, &raster[0], w*d);
	TIFFClose(tif);
	_TIFFfree(raster);	
}
