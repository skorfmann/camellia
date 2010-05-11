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

    Copyright (c) 2002-2010, Ecole des Mines de Paris - Centre de Robotique
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
#include <assert.h>
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795 
#endif
#include "camellia.h"
#include "camellia_internals.h"
#ifdef __SSE2__
#include <emmintrin.h>
#endif

int camPatchSizeParam = 32; //38; ///8*3;
double camSigmaParam = 7;
int camAbsCoeff = 5;
int camColorCoeff = 20;
int camHaarFilterSizeParam = 5;
int camHaarFilterSpaceParam = 0;
int camGradientSaturationParam = 100;
double camSigma2Param = 0.18; //2.2;

int camBuildGaussianFilter(CamImage *image, double sigma);

static int camKPNbAttPoints[20 * 20] = {-1}, camKPAttPoint[20 * 20 * 4], camKPCoeff[20 * 20 * 4];

/*
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
			camKPCoeff[i * 4 + 3] = ((y - 2) % 5) * ((x - 2) % 5);
		    }
		}
	    }
	}
    }
}
*/

void camKeypointsInternalsPrepareDescriptor()
{
    int x, y, i, j;

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
		j = 0;
		if (x >= 2 && y >= 2) {
		    camKPAttPoint[i * 4 + j] = ((y + 3) / 5 - 1) * 4 + (x + 3) / 5 - 1;
		    camKPCoeff[i * 4 + j] = (5 - ((y + 3) % 5)) * (5 - ((x + 3) % 5));
		    if (camKPCoeff[i * 4 + j]) j++;
		}
		if (x < 17 && y >= 2) {
		    camKPAttPoint[i * 4 + j] = ((y + 3) / 5 - 1) * 4 + (x + 3) / 5;
		    camKPCoeff[i * 4 + j] = (5 - ((y + 3) % 5)) * ((x + 3) % 5);
		    if (camKPCoeff[i * 4 + j]) j++;
		}
		if (x >= 2 && y < 17) {
		    camKPAttPoint[i * 4 + j] = ((y + 3) / 5) * 4 + (x + 3) / 5 - 1;
		    camKPCoeff[i * 4 + j] = ((y + 3) % 5) * (5 - ((x + 3) % 5));
		    if (camKPCoeff[i * 4 + j]) j++;
		}
		if (x < 17 && y < 17) {
		    camKPAttPoint[i * 4 + j] = ((y + 3) / 5) * 4 + (x + 3) / 5;
		    camKPCoeff[i * 4 + j] = ((y + 3) % 5) * ((x + 3) % 5);
		    if (camKPCoeff[i * 4 + j]) j++;
		}
		camKPNbAttPoints[i] = j;
		/*
		 * for (k = 0; k < j; k++) {
		    printf("(%d -> %d)", camKPAttPoint[i * 4 + k], camKPCoeff[i * 4 + k]);
		}
		printf("\n"); */
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
    CAM_INT16 *imptr_h, *imptr_v, *tmpimptr_h, *tmpimptr_v; 
    int val_h, val_v, val2_h, val2_v;

    const int xpp[4] = {-1, 1, 1, -1};
    const int ypp[4] = {-1, -1, 1, 1};

    int dx[22][22], dy[22][22];
    int w, sx, sy, xp, yp, xs, ys, l1, fsx, fsy, shift;
    CAM_INT32 *ptry[20], incx[20], *ptr;
//#define CHECK_CODE
#ifdef CHECK_CODE
    double sumd, cx, cy;
#endif
    
    static int nb = 0;
    char filename[256];
    FILE *handle;
   
    if (source->depth == 32) {
	// This is an integral image
	// w is the width of the patch to analyze, multiplied by 16
        w = point->scale * camPatchSizeParam;
        // sx is the size of the Haar filter, in pixels
	sx = ((w * camHaarFilterSizeParam) / 100) >> 4;
        if (sx <= 1) sx = 2;
	// fsx is the space between positive and negative parts of the Haar filter. It is most possibly 0.
        fsx = ((w * camHaarFilterSpaceParam) / 100) >> 4;	
	// Check boundaries
	if (point->x - (((19 * w) / 40) >> 4) - sx - fsx <= 1 || point->y - (((19 * w) / 40) >> 4) - sx - fsx <= 1 ||
	    point->x + (((19 * w) / 40) >> 4) + sx + fsx >= source->width ||
	    point->y + (((19 * w) / 40) >> 4) + sx + fsx >= source->height) {
	    //printf("Out of boundaries : %d, %d, %d\n", point->x, point->y, point->scale);
	    point->size = 0;
	    return 0;
	}
	l1 = source->widthStep >> 2;
	sy = sx * l1;
	fsy = fsx * l1;
	// Estimate shift : shift = log2(sx)
	shift = 0; i = sx;
	do { shift++; i >>= 1; } while (i != 0);
	shift *= 2;
	shift -= 8; if (shift < 0) shift = 0;

	for (channel = 0; channel < source->nChannels; channel++) {

	    for (i = 0; i < 20; i++) {
		ptry[i] = (CAM_INT32*)(source->imageData + (channel * source->height + point->y + ((w * ((i - 10) * 2 + 1) / 40) >> 4)) * source->widthStep);
		incx[i] = point->x + ((w * ((i - 10) * 2 + 1) / 40) >> 4);
	    }

	    // Compute the gradients
	    if (camHaarFilterSpaceParam > 0) {
		for (y = 0; y < 20; y++) {
		    for (x = 0; x < 20; x++) {
			ptr = ptry[y] + incx[x];
			dx[y][x] = *(ptr + sx + sy + fsx) - *(ptr + sy + fsx) - *(ptr + sx - sy + fsx) + *(ptr - sy + fsx) -
			    (*(ptr + sy - 1 - fsx) - *(ptr - sx + sy - 1 - fsx) - *(ptr - sy - 1 - fsx) + *(ptr - sx - sy - 1 - fsx));
			dy[y][x] = *(ptr + sx + sy + fsy) - *(ptr + sx + fsy) - *(ptr - sx + sy + fsy) + *(ptr - sx + fsy) -
			    (*(ptr + sx - l1 - fsy) - *(ptr - sx - l1 - fsy) - *(ptr + sx - sy - l1 - fsy) + *(ptr - sx - sy - l1 - fsy));
		    }
		}	    
	    } else {
		for (y = 0; y < 20; y++) {
		    for (x = 0; x < 20; x++) {
			ptr = ptry[y] + incx[x];
			dx[y][x] = *(ptr + sx + sy) - *(ptr + sy) - *(ptr + sx - sy - l1) + *(ptr - sy - l1) -
			    (*(ptr + sy - 1) - *(ptr - sx + sy - 1) - *(ptr - sy - 1 - l1) + *(ptr - sx - sy - 1 - l1));
			dy[y][x] = *(ptr + sx + sy) - *(ptr + sx) - *(ptr - sx + sy - 1) + *(ptr - sx - 1) -
			    (*(ptr + sx - l1) - *(ptr - sx - l1 - 1) - *(ptr + sx - sy - l1) + *(ptr - sx - sy - l1 - 1));
		    }
		}	    
	    }		
	    for (y = 0; y < 20; y++) {
		for (x = 0; x < 20; x++) {
		    dx[y][x] >>= shift;
		    dy[y][x] >>= shift;
		}
	    }
	    // Gaussian filtering
	    for (y = 0; y < 20; y++) {
		imptr_h = (CAM_INT16*)(filter->imageData + y * filter->widthStep);
		for (x = 0; x < 20; x++) {
		    assert(abs(dx[y][x]) < (1 << 20));
		    dx[y][x] *= ((*imptr_h) >> 6);
		    assert(abs(dx[y][x]) < (1 << 30));
		    dx[y][x] >>= 10;	
		    assert(abs(dy[y][x]) < (1 << 20));
		    dy[y][x] *= ((*imptr_h) >> 6);
		    assert(abs(dy[y][x]) < (1 << 30));
		    dy[y][x] >>= 10;	
		    imptr_h++;
		}
	    }	    

	    // Reduce effect of overweighted gradients
	    sum = 0;
	    for (y = 0; y < 20; y++) {
		for (x = 0; x < 20; x++) {
		    sum += abs(dx[y][x]) + abs(dy[y][x]);     
		}
	    }
	    i = (sum >> 10) * camGradientSaturationParam;
	    for (y = 0; y < 20; y++) {
		for (x = 0; x < 20; x++) {
		    if (dx[y][x] > i) dx[y][x] = i;
		    else if (dx[y][x] < -i) dx[y][x] = -i;
		    if (dy[y][x] > i) dy[y][x] = i;
		    else if (dy[y][x] < -i) dy[y][x] = -i;
		}
	    } 
	    
	    if (options & CAM_SIMPLE_SUM) {
		if (channel == 0) {
		    point->size = 0;
		    for (y = 0, ys = 0; y < 4; y++, ys += 5) {
			for (x = 0, xs = 0; x < 4; x++, xs += 5) {
			    sdx = 0; sdy = 0; sabs_dx = 0; sabs_dy = 0;	
			    for (yp = 0; yp < 5; yp++) {
				for (xp = 0; xp < 5; xp++) {
				    j = dx[ys + yp][xs + xp];
				    sdx += j;
				    sabs_dx += abs(j);
				    j = dy[ys + yp][xs + xp];
				    sdy += j;
				    sabs_dy += abs(j);
				}
			    }
			    point->descriptor[point->size++] = sdx;
			    point->descriptor[point->size++] = sdy;
			    point->descriptor[point->size++] = sabs_dx; 
			    point->descriptor[point->size++] = sabs_dy;
			}
		    }
		} else {
		    for (y = 0, ys = 0; y < 4; y++, ys += 5) {
			for (x = 0, xs = 0; x < 4; x++, xs += 5) {
			    sdx = 0; sdy = 0;	
			    for (yp = 0; yp < 5; yp++) {
				for (xp = 0; xp < 5; xp++) {
				    j = dx[ys + yp][xs + xp];
				    sdx += j;
				    j = dy[ys + yp][xs + xp];
				    sdy += j;
				}
			    }
			    point->descriptor[point->size++] = sdx;
			    point->descriptor[point->size++] = sdy;
			}
		    }
		}
	    } else {
		// Descriptor with bilinear interpolation
		for (i = 0; i < 16; i++) {
		    dxt[i] = 0; dyt[i] = 0; abs_dxt[i] = 0; abs_dyt[i] = 0;
		}

#if 0
		for (i = 0; i < 16; i++) {
		    cx = 2 + (i % 4) * 5;
		    cy = 2 + (i / 4) * 5;
		    for (y = 0; y < 20; y++) {
			for (x = 0; x < 20; x++) {
			    j = 255 * (1.0 - camSigma2Param * sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy)));
			    if (j < 0) j = 0;
			    //j = (255 * exp(-((x-cx) * (x-cx) + (y-cy) * (y-cy)) / (2 * camSigma2Param * camSigma2Param)));
			    val_v = j * dx[y][x];
			    val_h = j * dy[y][x];
			    dxt[i] += val_v;
			    dyt[i] += val_h;
			    if (channel == 0) {
				abs_dxt[i] += abs(val_v);
				abs_dyt[i] += abs(val_h);
			    }
			}
		    }
		}
#else
		for (y = 0, i = 0; y < 20; y++) {
		    for (x = 0; x < 20; x++, i++) {
			val_v = dx[y][x];
			val_h = dy[y][x];
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
		}
#endif

		if (channel == 0) {
		    if (options & CAM_DESC_SEP_TEXTURE) {
			for (j = 0, i = 0; j < 16; j++, i += 4) {
			    point->descriptor[i] = dxt[j] << 3;
			    point->descriptor[i + 1] = dyt[j] << 3;
			    point->descriptor[i + 2] = (abs_dxt[j] - abs(dxt[j])) * camAbsCoeff; 
			    point->descriptor[i + 3] = (abs_dyt[j] - abs(dyt[j])) * camAbsCoeff;
			}
		    } else {
			for (j = 0, i = 0; j < 16; j++, i += 4) {
			    point->descriptor[i] = dxt[j] << 3;
			    point->descriptor[i + 1] = dyt[j] << 3;
			    point->descriptor[i + 2] = abs_dxt[j] * camAbsCoeff; 
			    point->descriptor[i + 3] = abs_dyt[j] * camAbsCoeff;
			}
		    }
		    point->size = 64;
		    
		    // Test average on regions
/*
		    if (source->nChannels == 1) {
			for (y = 0; y < 20; y++) {
			    for (x = 0; x < 20; x++) {
				ptr = ptry[y] + incx[x];
				dx[y][x] = *(ptr + sx + sy) - *(ptr + sx - sy - l1) - *(ptr - sx - 1 + sy) + *(ptr - sx - 1 - sy - l1);
			    }
			}
			for (i = 0; i < 16; i++) dxt[i] = 0;
			for (y = 0, i = 0; y < 20; y++) {
			    for (x = 0; x < 20; x++, i++) {
				for (j = 0; j < camKPNbAttPoints[i]; j++) {
				    idx = camKPAttPoint[i * 4 + j];
				    coeff = camKPCoeff[i * 4 + j];
				    dxt[idx] += coeff * dx[y][x];
				}
			    }
			}
			// Add these averages to the descriptor, with a given weight...

		    }	
*/

		} else {
		    for (j = 0; j < 16; j++) {
			point->descriptor[point->size++] = dxt[j] * camColorCoeff;
			point->descriptor[point->size++] = dyt[j] * camColorCoeff;
		    }
		}
	    }
	}
    } else {	
	if (source->nChannels == 3) {
	    camAllocateYUVImage(&rotated, 20, 20);
	} else {
	    camAllocateImage(&rotated, 20, 20, source->depth);
	}

	// Rotate the image
	params.perspective=0;
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

	if ((options & CAM_DESC_DEBUG) && nb < 20) {
	    sprintf(filename, "output/rotated%d.pgm", nb++);
	    camSavePGM(&rotated, filename);
	}

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

	    if (options & CAM_SIMPLE_SUM) {
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
			    point->descriptor[i + 1] = sdy;
			    point->descriptor[i + 2] = sabs_dy;
			    point->descriptor[i + 3] = sabs_dx; 
			    i += 4;
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
		    if (options & CAM_DESC_SEP_TEXTURE) {
			for (j = 0, i = 0; j < 16; j++, i += 4) {
			    point->descriptor[i] = dxt[j] << 3;
			    point->descriptor[i + 1] = dyt[j] << 3;
			    point->descriptor[i + 2] = (abs_dxt[j] - abs(dxt[j])) * camAbsCoeff; 
			    point->descriptor[i + 3] = (abs_dyt[j] - abs(dyt[j])) * camAbsCoeff;
			}
		    } else {
			for (j = 0, i = 0; j < 16; j++, i += 4) {
			    point->descriptor[i] = dxt[j] << 3;
			    point->descriptor[i + 1] = dyt[j] << 3;
			    point->descriptor[i + 2] = abs_dxt[j] * camAbsCoeff; 
			    point->descriptor[i + 3] = abs_dyt[j] * camAbsCoeff;
			}
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
    // "à la OpenCV"
#if 0
    start = 0;
    end = point->size;
    for (i = start; i < end; i += 4) {
	sumd = 0;
	for (j = 0; j < 4; j++) {
	    sumd += (double)point->descriptor[i + j] * point->descriptor[i + j];
	}
	if (sumd == 0) {
	    for (j = 0; j < 4; j++) {
		point->descriptor[i + j] = 0;
	    }
	} else {
	    sumd = 1./sqrt(sumd);
	    for (j = 0; j < 4; j++) {
		point->descriptor[i + j] = point->descriptor[i + j] * sumd * 100000.;
	    }
	}
    }
#endif

#if 1 
#ifdef CHECK_CODE
//#define CHECK_L2_NORM
#ifdef CHECK_L2_NORM 
     if (options & CAM_DESC_SEP_NORM) {
	for (j = 0; j < ((source->nChannels == 1)?1:2); j++) {
	    start = j * 64;
	    end = start + 64;
	    sumd = 0;
	    for (i = start; i < end; i++) {
		sumd += (double)point->descriptor[i] * point->descriptor[i];
	    }
	    if (sumd == 0) {
		for (i = start; i < end; i++) {
		    point->descriptor[i] = 0;
		}
	    } else {
		sumd = 1./sqrt(sumd);
		for (i = start; i < end; i++) {
		    point->descriptor[i] = point->descriptor[i] * sumd * 100000.;
		}	
	    }
	}
    } else {
	start = 0;
	end = point->size;
	sumd = 0;
	for (i = start; i < end; i++) {
	    sumd += (double)point->descriptor[i] * point->descriptor[i];
	}
	if (sumd == 0) {
	    for (i = start; i < end; i++) {
		point->descriptor[i] = 0;
	    }
	} else {
	    sumd = 1./sqrt(sumd);
	    for (i = start; i < end; i++) {
		point->descriptor[i] = point->descriptor[i] * sumd * 100000.;
	    }
	}
    }
#else
     if (options & CAM_DESC_SEP_NORM) {
	for (j = 0; j < ((source->nChannels == 1)?1:2); j++) {
	    start = j * 64;
	    end = start + 64;
	    sum = 0;
	    for (i = start; i < end; i++) {
		sum += abs(point->descriptor[i]);
	    }
	    if (sum == 0) {
		for (i = start; i < end; i++) {
		    point->descriptor[i] = 0;
		}
	    } else {
		for (i = start; i < end; i++) {
		    point->descriptor[i] = ((double)point->descriptor[i] * (double)100000.0) / (double)sum;
		}	
	    }
	}
    } else {
	start = 0;
	end = point->size;
	sum = 0;
	for (i = start; i < end; i++) {
	    sum += abs(point->descriptor[i]);
	}
	if (sum == 0) {
	    for (i = start; i < end; i++) {
		point->descriptor[i] = 0;
	    }
	} else {
	    for (i = start; i < end; i++) {
		point->descriptor[i] = ((double)point->descriptor[i] * (double)100000.0) / (double)sum;
	    }
	}
    }
#endif // CHECK_L2_NORM
#else // CHECK_CODE
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
	//sum *= camColorCoeff;
	//sum >>= 3;
	if (sum != 0) sum = (1 << 30) / sum; // Division
	for (i = start; i < end; i++) {
	    point->descriptor[i] = (point->descriptor[i] * sum) >> 6;
	}
    }
#endif // CHECK_CODE  
#endif
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

