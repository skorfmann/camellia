/***************************************
 *
 *  Camellia Image Processing Library
 *

    The Camellia Image Processing Library is an open source low-level image processing library.
    As it uses the IplImage structure to describe images, it is a good replacement to the IPL (Intel) library
    and a good complement to the OpenCV library. It includes a lot of functions for image processing
    (filtering, morphological mathematics, labeling, warping, loading/saving images, etc.),
    some of them being highly optimized; It is also cross-platform and robust. It is doxygen-documented
    and examples of use are provided.

    This software library is an outcome of the Camellia european project (IST-2001-34410).
    It was developped by the Ecole des Mines de Paris (ENSMP), in coordination with
    the other partners of the project.

  ==========================================================================

    Copyright (c) 2002-2008, Ecole des Mines de Paris - Centre de Robotique
    All rights reserved.

    Redistribution and use in integral and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

        * Redistributions of integral code must retain the above copyright
          notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
          notice, this list of conditions and the following disclaimer
          in the documentation and/or other materials provided with the distribution.
        * Neither the name of the Ecole des Mines de Paris nor the names of
          its contributors may be used to endorse or promote products
          derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
    THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
    PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  ==========================================================================
*/

/* CamKeypoints descriptor / signature implementation
 * C code */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795 
#endif
#include "camellia.h"
#include "camellia_internals.h"
#ifdef __SSE2__
#include <emmintrin.h>
#endif

extern int camPatchSizeParam;
extern int camSigmaParam;
int camAbsCoeff = 8;
int camColorCoeff = 8;
int camHaarFilterSizeParam = 5;

int camBuildGaussianFilter(CamImage *image, double sigma);

static int camKPNbAttPoints[20 * 20] = {-1}, camKPAttPoint[20 * 20 * 4], camKPCoeff[20 * 20 * 4];

void camKeypointsInternalsPrepareDescriptor()
{
    int x, y, i;

    // Test whether this initialization has already been done before
    if (camKPNbAttPoints[0] == -1) {

	for (i = 0; i < 20 * 20 * 4; i++) {
	    camKPAttPoint[i] = 0;
	    camKPCoeff[i] = 0;
	}
	// For all the camKPAttPoints, find the attraction camKPAttPoints
	// and compute the bilinear interpolation camKPCoefficients
	for (y = 0, i = 0; y < 20; y++) {
	    for (x = 0; x < 20; x++, i++) {
		if (x < 2) {
		    if (y <= 2) {
			camKPNbAttPoints[i] = 1;
			camKPAttPoint[i * 4] = 0;
			camKPCoeff[i * 4] =  25;
		    } else if (y >= 17) {
			camKPNbAttPoints[i] = 1;
			camKPAttPoint[i * 4] = 12;
			camKPCoeff[i * 4] = 25;
		    } else if ((y - 2) % 5 == 0) {
			camKPNbAttPoints[i] = 1;
			camKPAttPoint[i * 4] = ((y - 2) / 5) * 4;
			camKPCoeff[i * 4] = 25;
		    } else {
			camKPNbAttPoints[i] = 2;
			camKPAttPoint[i * 4] = ((y - 2) / 5) * 4;
			camKPAttPoint[i * 4 + 1] = ((y - 2) / 5 + 1) * 4;
			camKPCoeff[i * 4] = 25 - 5 * ((y - 2) % 5);
			camKPCoeff[i * 4 + 1] = 5 * ((y - 2) % 5);
		    }
		} else if (x > 17) {
		    if (y <= 2) {
			camKPNbAttPoints[i] = 1;
			camKPAttPoint[i * 4] = 3;
			camKPCoeff[i * 4] =  25;
		    } else if (y >= 17) {
			camKPNbAttPoints[i] = 1;
			camKPAttPoint[i * 4] = 15;
			camKPCoeff[i * 4] = 25;
		    } else if ((y - 2) % 5 == 0) {
			camKPNbAttPoints[i] = 1;
			camKPAttPoint[i * 4] = ((y - 2) / 5) * 4 + 3;
			camKPCoeff[i * 4] = 25;
		    } else {
			camKPNbAttPoints[i] = 2;
			camKPAttPoint[i * 4] = ((y - 2) / 5) * 4 + 3;
			camKPAttPoint[i * 4 + 1] = ((y - 2) / 5 + 1) * 4 + 3;
			camKPCoeff[i * 4] = 25 - 5 * ((y - 2) % 5);
			camKPCoeff[i * 4 + 1] = 5 * ((y - 2) % 5);
		    }
		} else if ((x - 2) % 5 == 0) {
		    if (y <= 2) {
			camKPNbAttPoints[i] = 1;
			camKPAttPoint[i * 4] = (x - 2) / 5;
			camKPCoeff[i * 4] =  25;
		    } else if (y >= 17) {
			camKPNbAttPoints[i] = 1;
			camKPAttPoint[i * 4] = 12 + (x - 2) / 5;
			camKPCoeff[i * 4] = 25;
		    } else if ((y - 2) % 5 == 0) {
			camKPNbAttPoints[i] = 1;
			camKPAttPoint[i * 4] = ((y - 2) / 5) * 4 + (x - 2) / 5;
			camKPCoeff[i * 4] = 25;
		    } else {
			camKPNbAttPoints[i] = 2;
			camKPAttPoint[i * 4] = ((y - 2) / 5) * 4 + (x - 2) / 5;
			camKPAttPoint[i * 4 + 1] = ((y - 2) / 5 + 1) * 4 + (x - 2) / 5;
			camKPCoeff[i * 4] = 25 - 5 * ((y - 2) % 5);
			camKPCoeff[i * 4 + 1] = 5 * ((y - 2) % 5);
		    }
		} else {
		    if (y <= 2) {
			camKPNbAttPoints[i] = 2;
			camKPAttPoint[i * 4] = (x - 2) / 5;
			camKPAttPoint[i * 4 + 1] = (x - 2) / 5 + 1;
			camKPCoeff[i * 4] = 25 - 5 * ((x - 2) % 5);
			camKPCoeff[i * 4 + 1] = 5 * ((x - 2) % 5);
		    } else if (y >= 17) {
			camKPNbAttPoints[i] = 2;
			camKPAttPoint[i * 4] = 12 + (x - 2) / 5;
			camKPAttPoint[i * 4 + 1] = 13 + (x - 2) / 5;
			camKPCoeff[i * 4] = 25 - 5 * ((x - 2) % 5);
			camKPCoeff[i * 4 + 1] = 5 * ((x - 2) % 5);
		    } else if ((y - 2) % 5 == 0) {
			camKPNbAttPoints[i] = 2;
			camKPAttPoint[i * 4] = ((y - 2) / 5) * 4 + (x - 2) / 5;
			camKPAttPoint[i * 4 + 1] = ((y - 2) / 5) * 4 + (x - 2) / 5 + 1;
			camKPCoeff[i * 4] = 25 - 5 * ((x - 2) % 5);
			camKPCoeff[i * 4 + 1] = 5 * ((x - 2) % 5);
		    } else {
			camKPNbAttPoints[i] = 4;
			camKPAttPoint[i * 4] = ((y - 2) / 5) * 4 + (x - 2) / 5;
			camKPAttPoint[i * 4 + 1] = ((y - 2) / 5) * 4 + (x - 2) / 5 + 1;
			camKPAttPoint[i * 4 + 2] = ((y - 2) / 5 + 1) * 4 + (x - 2) / 5;
			camKPAttPoint[i * 4 + 3] = ((y - 2) / 5 + 1) * 4 + (x - 2) / 5 + 1;
			camKPCoeff[i * 4] = (5 - ((y - 2) % 5)) * (5 - ((x - 2) % 5));
			camKPCoeff[i * 4 + 1] = (5 - ((y - 2) % 5)) * ((x - 2) % 5);
			camKPCoeff[i * 4 + 2] = ((y - 2) % 5) * (5 - ((x - 2) % 5));
			camKPCoeff[i * 4 + 3] = (5 - ((y - 2) % 5)) * (5 - ((x - 2) % 5));
		    }
		}
	    }
	}
    }
}

int camKeypointDescriptor(CamKeypoint *point, CamImage *source, CamImage *filter, int options)
{
    CamWarpingParams params;
    CamImage rotated, filtered_h, filtered_v;
    CamArithmParams params2;

    int i, j, x, y, sum;
    int costheta = (int)(cos(-point->angle * 2 * M_PI / 360) * 65536.0);
    int sintheta = (int)(sin(-point->angle * 2 * M_PI / 360) * 65536.0);
    int scale, channel, start, end;
    CamROI roix;

    // For SURF-like descriptor :
    int sdx, sdy, sabs_dx, sabs_dy;
    CamROI roi;

    // For bi-linear interpolated descriptor :
    int coeff, idx, dxt[16], dyt[16], abs_dxt[16], abs_dyt[16];
    CAM_INT16 val_h, val_v, val2_h, val2_v, *imptr_h, *imptr_v, *tmpimptr_h, *tmpimptr_v; 

    const int xpp[4] = {-1, 1, 1, -1};
    const int ypp[4] = {-1, -1, 1, 1};

    int dx[20][20], dy[20][20];
    int w, sx, sy, xp, yp, xs, ys;
    CAM_INT32 *ptry[20], incx[20], *ptr;
#if 0
    static int nb = 0;
    char filename[256];
    FILE *handle;
#endif
   
    if (source->depth == 32) {
	// This is an integral image
	w = (point->scale * camPatchSizeParam) >> 4;
	sx = (w * camHaarFilterSizeParam) / 100; 
	for (i = 0; i < 20; i++) {
	    ptry[i] = (CAM_INT32*)(source->imageData + (point->y + w * ((i - 10) * 2 + 1) / 40) * source->widthStep);
	    incx[i] = point->x + w * ((i - 10) * 2 + 1) / 40;
	}
	sy = sx * (source->widthStep >> 2);

	// Compute the gradients
	for (y = 0; y < 20; y++) {
	    for (x = 0; x < 20; x++) {
		ptr = ptry[y] + incx[x];
		dx[x][y] = *(ptr + sx + sy) - *(ptr + sy) - *(ptr + sx - sy) + *(ptr - sy) -
		    (*(ptr + sy) - *(ptr - sx + sy) - *(ptr - sy) + *(ptr - sx - sy));
		dy[x][y] = *(ptr + sx + sy) - *(ptr + sx) - *(ptr - sx + sy) + *(ptr - sx) -
		    (*(ptr + sx) - *(ptr - sx) - *(ptr + sx - sy) + *(ptr - sx - sy));
		// Filtrage gaussien
	    }
	}	    

	// Build the descriptor
	if (options & CAM_DESC_SURF_LIKE) {
	    /*if (channel == 0)*/ {
		i = 0;
		for (y = 0, ys = 0; y < 4; y++, ys += 5) {
		    for (x = 0, xs = 0; x < 4; x++, xs += 5) {
			sdx = 0; sdy = 0; sabs_dx = 0; sabs_dy = 0;	
			for (yp = 0; yp < 5; yp++) {
			    for (xp = 0; xp < 5; xp++) {
				j = dx[xs + xp][ys + yp];
				sdx += j;
				sabs_dx += abs(j);
				j = dy[xs + xp][ys + yp];
				sdy += j;
				sabs_dy += abs(j);
			    }
			}
			point->descriptor[i] = sdx;
			point->descriptor[32 + i] = sabs_dx; 
			i++;
			point->descriptor[i] = sdy;
			point->descriptor[32 + i] = sabs_dy;
			i++;
		    }
		}
		point->size = 64;
	    } /* else {
		for (y = 0; y < 4; y++) {
		    for (x = 0; x < 4; x++) {
			roi.xOffset = x * 5;
			roi.yOffset = y * 5;
			dx = camSumOfPixels(&filtered_v);
			dy = camSumOfPixels(&filtered_h);
			point->descriptor[point->size++] = dx;
			point->descriptor[point->size++] = dy;
		    }
		}
	    }*/
	}
    } else {	
	if (source->nChannels == 3) {
	    camAllocateYUVImage(&rotated, 20, 20);
	} else {
	    camAllocateImage(&rotated, 20, 20, source->depth);
	}

	// Rotate the image
	params.perspective=1;
	params.interpolation=1;
	scale = point->scale * camPatchSizeParam;
	for (i = 0; i < 4; i++) {
	    params.p[i].x = costheta * xpp[i] - sintheta * ypp[i];
	    params.p[i].y = sintheta * xpp[i] + costheta * ypp[i];
	    params.p[i].x = (params.p[i].x * scale) >> 5;
	    params.p[i].y = (params.p[i].y * scale) >> 5;
	    params.p[i].x += (point->x << 16) + 32768;
	    params.p[i].y += (point->y << 16) + 32768;
	}
	camWarping(source, &rotated, &params);

#if 0
	sprintf(filename, "output/rotated%d.pgm", nb++);
	camSavePGM(&rotated, filename);
#endif

	// We now have the rotated image
	// Let's filter it...
	camAllocateImage(&filtered_h, 20, 20, CAM_DEPTH_16S);
	camAllocateImage(&filtered_v, 20, 20, CAM_DEPTH_16S);
	roix.xOffset = 0; roix.yOffset = 0; roix.width = 20; roix.height = 20;
	rotated.roi=&roix;

	for (channel = 0; channel < source->nChannels; channel++) {
	    roix.coi = channel + 1;
	    camFixedFilter(&rotated, &filtered_h, CAM_SCHARR_H);
	    camFixedFilter(&rotated, &filtered_v, CAM_SCHARR_V);
	    params2.operation = CAM_ARITHM_MUL;
	    params2.c1 = 16;
	    camDyadicArithm(&filtered_h, filter, &filtered_h, &params2);
	    camDyadicArithm(&filtered_v, filter, &filtered_v, &params2);

	    /*
	       camAbs(&filtered_h, &filtered_h);
	       camAbs(&filtered_v, &filtered_v);
	       camSavePGM(&filtered_h, "output/filtered_h.pgm");
	       camSavePGM(&filtered_v, "output/filtered_v.pgm");
	       */

	    if (options & CAM_DESC_SURF_LIKE) {
		filtered_h.roi = &roi;
		filtered_v.roi = &roi;
		roi.width = 5;
		roi.height = 5;

		if (channel == 0) {
		    i = 0;
		    for (y = 0; y < 4; y++) {
			for (x = 0; x < 4; x++) {
			    roi.xOffset = x * 5;
			    roi.yOffset = y * 5;
			    sdx = camSumOfPixels(&filtered_v);
			    sdy = camSumOfPixels(&filtered_h);
			    sabs_dx = camAbs(&filtered_v, &filtered_v);
			    sabs_dy = camAbs(&filtered_h, &filtered_h); 
			    point->descriptor[i] = sdx;
			    point->descriptor[32 + i] = sabs_dx; 
			    i++;
			    point->descriptor[i] = sdy;
			    point->descriptor[32 + i] = sabs_dy;
			    i++;
			}
		    }
		    point->size = 64;
		} else {
		    for (y = 0; y < 4; y++) {
			for (x = 0; x < 4; x++) {
			    roi.xOffset = x * 5;
			    roi.yOffset = y * 5;
			    sdx = camSumOfPixels(&filtered_v);
			    sdy = camSumOfPixels(&filtered_h);
			    point->descriptor[point->size++] = sdx;
			    point->descriptor[point->size++] = sdy;
			}
		    }
		}
	    } else {
		for (i = 0; i < 16; i++) {
		    dxt[i] = 0; dyt[i] = 0; abs_dxt[i] = 0; abs_dyt[i] = 0;
		}

		tmpimptr_h = (CAM_INT16*)filtered_h.imageData;
		tmpimptr_v = (CAM_INT16*)filtered_v.imageData;
		for (y = 0, i = 0; y < 20; y++) {
		    imptr_h = tmpimptr_h;
		    imptr_v = tmpimptr_v;
		    for (x = 0; x < 20; x++, i++, imptr_h++, imptr_v++) {
			val_h = *imptr_h;
			val_v = *imptr_v;
			for (j = 0; j < camKPNbAttPoints[i]; j++) {
			    idx = camKPAttPoint[i * 4 + j];
			    coeff = camKPCoeff[i * 4 + j];
			    val2_h = coeff * val_h;
			    val2_v = coeff * val_v;
			    dxt[idx] += val2_v;
			    dyt[idx] += val2_h;
			    if (channel == 0) {
				abs_dxt[idx] += abs(val2_v);
				abs_dyt[idx] += abs(val2_h);
			    }
			}
		    }
		    tmpimptr_h = (CAM_INT16*)(((char*)tmpimptr_h) + filtered_h.widthStep);
		    tmpimptr_v = (CAM_INT16*)(((char*)tmpimptr_v) + filtered_v.widthStep);
		}

		if (channel == 0) {
		    for (j = 0, i = 0; j < 16; j++) {
			point->descriptor[i] = dxt[j] << 3;
			point->descriptor[32 + i] = abs_dxt[j] * camAbsCoeff; 
			i++;
			point->descriptor[i] = dyt[j] << 3;
			point->descriptor[32 + i] = abs_dyt[j] * camAbsCoeff;
			i++;
		    }
		    point->size = 64;
		} else {
		    for (j = 0; j < 16; j++) {
			point->descriptor[point->size++] = dxt[j] * camColorCoeff;
			point->descriptor[point->size++] = dyt[j] * camColorCoeff;
		    }
		}
	    }
	}
#if 0
	// Save the signature
	sprintf(filename, "output/signature%d.txt", nb++);
	handle = fopen(filename, "wt");
	for (i = 0; i < point->size; i++) {
	    fprintf(handle, "%d\n", point->descriptor[i]);
	}
	fclose(handle);
#endif
	camDeallocateImage(&rotated);
	camDeallocateImage(&filtered_h);
	camDeallocateImage(&filtered_v);
    }
  
    // Normalization
    if (options & CAM_DESC_SEP_NORM) {
	for (j = 0; j < ((source->nChannels == 1)?1:2); j++) {
	    start = j * 64;
	    end = start + 64;
	    sum = 0;
	    for (i = start; i < end; i++) {
		sum += abs(point->descriptor[i]);
	    }
	    if (sum != 0) sum = (1 << 30) / sum; // Division
	    for (i = start; i < end; i++) {
		point->descriptor[i] = (point->descriptor[i] * sum) >> 6;
	    }	
	}
    } else {
	start = 0;
	end = point->size;
	sum = 0;
	for (i = start; i < end; i++) {
	    sum += abs(point->descriptor[i]);
	}
	if (sum != 0) sum = (1 << 30) / sum; // Division
	for (i = start; i < end; i++) {
	    point->descriptor[i] = (point->descriptor[i] * sum) >> 6;
	}
    }	
    return 1;
}

int camKeypointsDescriptor(CamKeypoints *points, CamImage *source, int options)
{
    int i;
    CamImage filter;

    camAllocateImage(&filter, 20, 20, CAM_DEPTH_16S);
    camBuildGaussianFilter(&filter, camSigmaParam);
    camKeypointsInternalsPrepareDescriptor();

    // Get the signature from the selected keypoints
    for (i = 0; i < points->nbPoints; i++) {
	camKeypointDescriptor(points->keypoint[i], source, &filter, options);
    }

    camDeallocateImage(&filter);
    return 1;
}

