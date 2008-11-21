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

/* CamKeypoints implementation
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

#define POST_SCALING

extern int camPatchSizeParam;
extern double camSigmaParam;

int camKeypointsSetParameters(int patchSize, double sigma)
{
    camPatchSizeParam = patchSize;
    camSigmaParam = sigma;
    return 1;
}

static const int CamScale[] = {3, 5, 7, 9, 13, 17, 25, 33, 49, 65};
static const int CamOffset[][2] = {
	{2, 1}, {3, 2}, {4, 2}, {6, 3}, {9, 4}, {11, 6}, {17, 8}, {22, 11}, {33, 16}, {43, 22}		
    };
static const int CamSampling[] = {0, 0, 0, 0, 1, 1, 2, 2, 3, 3};
static const int CamSamplingOffset[] = {0, 0, 0, 0, 0, 0, 1, 1, 3, 3};

#include "cam_keypoints_hessian_code.c"
#define CAM_FAST_APPROX_HESSIAN
#include "cam_keypoints_hessian_code.c"
#undef CAM_FAST_APPROX_HESSIAN
// Keypoints allocation
int camAllocateKeypoints(CamKeypoints *fpoints, int nbPoints)
{
    CAM_CHECK_ARGS(camAllocateKeypoints, fpoints != NULL);
    fpoints->nbPoints = 0;
    if (nbPoints) {
	fpoints->keypoint = (CamKeypoint**)malloc(sizeof(CamKeypoint*) * nbPoints);
	if (fpoints->keypoint == NULL) {
	    fpoints->allocated = 0;
	    camError("camAllocateKeypoints", "Memory allocation error");
	    return 0;
	}
    } else fpoints->keypoint = NULL;
    fpoints->bag = NULL;
    fpoints->allocated = nbPoints;
    return 1;
}

// Keypoints reallocation
int camReallocateKeypoints(CamKeypoints *fpoints, int nbPoints)
{
    CAM_CHECK_ARGS(camKeypointsReallocate, fpoints != NULL);
    if (fpoints->keypoint == NULL) {
	return camAllocateKeypoints(fpoints, nbPoints);
    }
    fpoints->keypoint = (CamKeypoint**)realloc(fpoints->keypoint, sizeof(CamKeypoint*) * nbPoints);
    if (fpoints->keypoint == NULL) {
	fpoints->nbPoints = 0;
	fpoints->allocated = 0;
	camError("camKeypointsReallocate", "Memory allocation error");
	return 0;
    }
    fpoints->allocated = nbPoints;
    return 1;
}

/// Keypoints deallocation
int camFreeKeypoints(CamKeypoints *fpoints)
{
    CAM_CHECK_ARGS(camKeypointsDeallocate, fpoints != NULL);
    if (fpoints->keypoint) free(fpoints->keypoint);
    if (fpoints->bag) 
#ifdef __SSE2__
	_mm_free(fpoints->bag);
#else
	free(fpoints->bag);
#endif
    fpoints->keypoint = NULL;
    fpoints->bag = NULL;
    fpoints->nbPoints = 0;
    fpoints->allocated = 0;
    return 1;
}

int camBuildGaussianFilter(CamImage *image, double sigma)
{
    int x,y;
    int width,height;
    CAM_INT16 *imptr,*tmpptr;
    CamInternalROIPolicyStruct iROI;
    DECLARE_MASK_MANAGEMENT;
    double cx, cy;

    CAM_CHECK(camBuildGaussianFilter,camInternalROIPolicy(image, NULL, &iROI, 1));
    CAM_CHECK_ARGS(camBuildGaussianFilter,(image->depth & CAM_DEPTH_MASK) == 16);

    // ROI (Region Of Interest) management
    width = iROI.srcroi.width;
    height = iROI.srcroi.height;
    imptr = (CAM_INT16*)iROI.srcptr;
    cx = iROI.srcroi.xOffset + (iROI.srcroi.width / 2.0) - 0.5;
    cy = iROI.srcroi.yOffset + (iROI.srcroi.height / 2.0) - 0.5;

    INIT_MASK_MANAGEMENT;

    for (y = 0; y < height; y++) {
	tmpptr = imptr; 
        BEGIN_MASK_MANAGEMENT(imptr = tmpptr + startx * iROI.srcinc;)
            for (x = startx; x < endx; x++, imptr += iROI.srcinc) {
                *imptr = (CAM_INT16)(32767 * exp(-((x-cx) * (x-cx) + (y-cy) * (y-cy)) / (2 * sigma * sigma)));
            }
        END_MASK_MANAGEMENT;
        imptr = tmpptr + iROI.srclinc;
    }
    
    camInternalROIPolicyExit(&iROI);
    return 1;    
}

#define CAM_NB_SECTORS 36
#define CAM_ORIENTATION_STAMP_SIZE 30

int camKeypointsInternalsPrepareOrientation()
{
    int i, s, x, y;
    double angle, v1[2], v2[2], v3[2];
    double dp1, dp2; // Dot products
    int xok[CAM_ORIENTATION_STAMP_SIZE * CAM_ORIENTATION_STAMP_SIZE], yok[CAM_ORIENTATION_STAMP_SIZE * CAM_ORIENTATION_STAMP_SIZE], nbok;
    FILE *handle;
    CamImage dummy;
#define OFFSET(x,y) (dummy.widthStep * (y) / 2 + (x)) 

    camAllocateImage(&dummy, CAM_ORIENTATION_STAMP_SIZE, CAM_ORIENTATION_STAMP_SIZE, CAM_DEPTH_16S);
    handle = fopen("src/cam_keypoints_sectors_code.c", "wt");
    fprintf(handle, "const int CamKeypointSectors[CAM_NB_SECTORS][MAX_PIX_PER_SECTOR] = {\n");
    for (s = 0; s < CAM_NB_SECTORS; s++) {
	nbok = 0;
	// Try to find out which points lie within the sector
	angle = 2 * M_PI * s / CAM_NB_SECTORS;
	v1[0] = cos(angle);
	v1[1] = sin(angle);
	angle = 2 * M_PI * (s + 1) / CAM_NB_SECTORS;
	v2[0] = cos(angle);
	v2[1] = sin(angle);
	// Test all the points...
	for (y = 0; y < CAM_ORIENTATION_STAMP_SIZE; y++) {
	    for (x = 0; x < CAM_ORIENTATION_STAMP_SIZE; x++) {
		v3[0] = x - (CAM_ORIENTATION_STAMP_SIZE/2 - 0.5); 
		v3[1] = (CAM_ORIENTATION_STAMP_SIZE/2 - 0.5) - y; // y goes upward
		if (v3[0] * v3[0] + v3[1] * v3[1] < 14 * 14) {
		    // This pixel has a good radius
		    // Does it lie within this sector ?
		    dp1 = - v3[0] * v1[1] + v3[1] * v1[0];
		    dp2 = - v3[0] * v2[1] + v3[1] * v2[0];
		    if (dp1 > 0 && dp2 < 0) {
		        // Yes, it does ! Record it !
			xok[nbok] = x; yok[nbok] = y; nbok++;
		    }	
		}
	    }
	}
	// OK. We have our set of points for the current sector.
	// Record the data
	fprintf(handle, "\t{");
	for (i = 0; i < nbok; i++) {
	    fprintf(handle, "%d, ", OFFSET(xok[i], yok[i]));
	}
	fprintf(handle, "-1}");
	if (s != CAM_NB_SECTORS - 1) fprintf(handle, ",\n");
    }
    fprintf(handle, "\n};\n");
    fclose(handle);
    camDeallocateImage(&dummy);
    return 1;
}

#define MAX_PIX_PER_SECTOR 20 
#include "cam_keypoints_sectors_code.c"

int camSortStrength(const void *p1x, const void *p2x)
{
    int *p1 = (int*)p1x;
    int *p2 = (int*)p2x;
    if (p1[0] < p2[0]) return 1;
    if (p1[0] == p2[0]) return 0;
    return -1;
}

// Returns 1 if two orientations have been selected
int camKeypointOrientation(CamImage *source, CamKeypointShort *point, CamImage *filter, CamKeypointShort *next_point)
{
    CamWarpingParams params;
    CamImage scaled, filtered_h, filtered_v;
    CamArithmParams params2;

    int i, c, d, x, y;
    int scale, angle, angle1, angle2, dangle;
    int vx[CAM_NB_SECTORS], vy[CAM_NB_SECTORS];
    int vxc[CAM_NB_SECTORS], vyc[CAM_NB_SECTORS];
    int strength[CAM_NB_SECTORS], stgx[CAM_NB_SECTORS * 2];
    int maximum;

    const int xp[4] = {-1, 1, 1, -1};
    const int yp[4] = {-1, -1, 1, 1};

    camAllocateImage(&scaled, CAM_ORIENTATION_STAMP_SIZE, CAM_ORIENTATION_STAMP_SIZE, source->depth);

    // Scale the stamp 
    params.perspective=0;
    params.interpolation=1;
    scale = (point->scale * camPatchSizeParam * 3) << 9;
    for (i = 0; i < 4; i++) {
	params.p[i].x = xp[i] * scale;
 	params.p[i].y = yp[i] * scale;
	params.p[i].x += (point->x << 16) + 32768;
	params.p[i].y += (point->y << 16) + 32768;
    }
    camWarping(source, &scaled, &params);
    
#if 0
    static int nbx = 0;
    char strx[256];
    sprintf(strx, "output/keypoints_rotation_%d.pgm", nbx++);
    camSavePGM(&scaled, strx);
#endif

    // We now have the scaled image
    // Let's filter it...
    camAllocateImage(&filtered_h, CAM_ORIENTATION_STAMP_SIZE, CAM_ORIENTATION_STAMP_SIZE, CAM_DEPTH_16S);
    camAllocateImage(&filtered_v, CAM_ORIENTATION_STAMP_SIZE, CAM_ORIENTATION_STAMP_SIZE, CAM_DEPTH_16S);
    camFixedFilter(&scaled, &filtered_h, CAM_SCHARR_H);
    camFixedFilter(&scaled, &filtered_v, CAM_SCHARR_V);
    params2.operation = CAM_ARITHM_MUL;
    params2.c1 = 16;
    camDyadicArithm(&filtered_h, filter, &filtered_h, &params2);
    camDyadicArithm(&filtered_v, filter, &filtered_v, &params2);

    // Great. Let's compute the sum of gradient for all sectors...
    for (i = 0; i < CAM_NB_SECTORS; i++) {
	vx[i] = 0; vy[i] = 0;
	c = 0;
	while (CamKeypointSectors[i][c] != -1) {
	    vx[i] += *(((CAM_INT16*)filtered_v.imageData) + CamKeypointSectors[i][c]);
	    vy[i] -= *(((CAM_INT16*)filtered_h.imageData) + CamKeypointSectors[i][c]);
	    c++;
	}
    }

    // Accumulation into pi/3 sectors
    for (i = 0; i < CAM_NB_SECTORS; i++) {
	vxc[i] = vx[i];
	vyc[i] = vy[i];
	d = i + 1;
	if (d == CAM_NB_SECTORS) d = 0;
	for (c = 1; c < CAM_NB_SECTORS / 6; c++) {
	    vxc[i] += vx[d];
	    vyc[i] += vy[d];
	    d++;
	    if (d == CAM_NB_SECTORS) d = 0;
	}
    }

#if 0
    FILE *handle;
    static int nb = 0;
    char str[256];
    sprintf(str, "output/keypoints_orientation_%d.txt", nb++);
    handle = fopen(str, "wt");
    for (i = 0; i < CAM_NB_SECTORS; i++) {
	fprintf(handle, "%d\n", strength[i]);
    }
    fclose(handle);
#endif

    camDeallocateImage(&scaled);
    camDeallocateImage(&filtered_h);
    camDeallocateImage(&filtered_v);
    
    // What is the best ?
    for (i = 0; i < CAM_NB_SECTORS; i++) {
	x = vxc[i] >> 4;
	y = vyc[i] >> 4;
	strength[i] = stgx[i * 2] = x * x + y * y;
	stgx[i * 2 + 1] = i;
    }
    qsort(stgx, CAM_NB_SECTORS, sizeof(int) * 2, camSortStrength);
    maximum = stgx[1];
    point->angle = (int)floor(atan2(vyc[maximum], vxc[maximum]) * 360 / (2 * M_PI) + 180 + 0.5);

    // Is there another angle that would deserve to be turned into a keypoint ?
    for (i = 1; i < CAM_NB_SECTORS; i++) {
	if (stgx[i * 2] > stgx[0] * 2 / 3) {
	    // Check that it is a local maximum
	    d = stgx[i * 2 + 1]; // Retrieve the index
	    if (d == 0) c = CAM_NB_SECTORS - 1; else c = d - 1;
	    if (stgx[i * 2] > strength[c]) {
		if (d == CAM_NB_SECTORS - 1) c = 0; else c = d + 1;
		if (stgx[i * 2] > strength[c]) {
		    // OK. This is a local maximum
		    angle = (int)floor(atan2(vyc[d], vxc[d]) * 360 / (2 * M_PI) + 180 + 0.5);
		    if (angle > point->angle) {
			angle1 = point->angle;
			angle2 = angle;
		    } else {
			angle1 = angle;
			angle2 = point->angle;
		    }
		    dangle = 360 - angle2 + angle1;
		    if (angle2 - angle1 < dangle) dangle = angle2 - angle1;
		    if (dangle > 60) {
			if (dangle >= 170 && dangle <= 190) {
			    // This is the opposite
			    // Let's keep only one of them
			    if (vyc[maximum] < 0) {
				point->angle = angle;
			    }
			} else {
			    // Add a new keypoint 
			    *next_point = *point;
			    next_point->angle = angle;
			    return 1;
			}
			break;
		    } 
		}
	    }
	    
	} else break;
    }
    return 0;
}

int camSortKeypointsShort(const void *p1x, const void *p2x)
{
    CamKeypointShort *p1 = (CamKeypointShort*)p1x;
    CamKeypointShort *p2 = (CamKeypointShort*)p2x;
    if (p1->value < p2->value) return 1;
    if (p1->value == p2->value) return 0;
    return -1;
}

#define CAM_MAX_SCALES 20
#define CAM_MAX_KEYPOINTS 100000

int camFastHessianDetector(CamImage *source, CamKeypoints *points, int nb_max_keypoints, int options)
{
    CamImage integral;
    CamImage results[sizeof(CamScale) / sizeof(CamScale[0])];
    int width, height;
    int i, j, k, scale;
    CamInternalROIPolicyStruct iROI;
    CamROI *roi, roix;
    const int nbScales = sizeof(CamScale) / sizeof(CamScale[0]);
    CamKeypointShort *keypoints;
    int nb_keypoints = 0, pnb_keypoints;

    int x, y, d, thr;
    int widthClusters, heightClusters;
    int nbClusters, cluster;
    int *nbPerCluster;
    CamKeypointShort ***firstPoint; 
    CamKeypointShort **ptrPoints;
    CamKeypointShort *point1, *point2;
    CamKeypointShort *bigOne, *smallOne;
    int nbClusters2scan, nbPoints2scan[5];
    CamKeypointShort **first2scan[5];

    int x1, x3, y1, y2, y3, p, num, den;
    CAM_UINT16 *ptr;
    CamImage filter;

    int multiplier, shift, coeff[CAM_MAX_SCALES];
    char str[64]; // for debug purpose

    // Parameters checking
    CAM_CHECK(camFastHessianDetector, camInternalROIPolicy(source, NULL, &iROI, 1));
    CAM_CHECK_ARGS(camFastHessianDetector, (source->depth & CAM_DEPTH_MASK) >= 8);
    CAM_CHECK_ARGS(camFastHessianDetector, points->allocated >= 0);
    CAM_CHECK_ARGS(camFastHessianDetector, source->nChannels == 1 || ((source->nChannels == 3) && (source->dataOrder == CAM_DATA_ORDER_PLANE)));
    width = iROI.srcroi.width;
    height = iROI.srcroi.height;

    // Algorithm initialization
    if (options & CAM_APPROX_HESSIAN) {
	for (scale = 0; scale < nbScales; scale++) {
	    coeff[scale] = (1 << (scale + 12)) / ((CamScale[scale]*3 - 2 * CamOffset[scale][0]) * CamScale[scale]);
	    //printf("%d = %d\n", scale, coeff[scale]);
	}
    } else {
	for (scale = 0; scale < nbScales; scale++) {
	    multiplier = 65536 / ((CamScale[scale]*3 - 2 * CamOffset[scale][0]) * CamScale[scale]);  
	    shift = 18 - 2 * scale;
	    if (shift >= 0) 
		coeff[scale] = (multiplier * multiplier) >> shift;
	    else
		coeff[scale] = (multiplier * multiplier) << (-shift); 
	}
    }

    integral.imageData = NULL;
    roi = source->roi;
    camSetMaxROI(&roix, source);
    roix.coi = 1;
    source->roi = &roix;
    camIntegralImage(source, &integral);
    integral.roi = &iROI.srcroi;
    source->roi = roi;

    // Allocate temp memory for keypoints
    keypoints = (CamKeypointShort*)malloc(CAM_MAX_KEYPOINTS * sizeof(CamKeypointShort));

    // Fast Hessian Detector for all scales
    for (scale = 0; scale < nbScales; scale++) {
	results[scale].imageData = NULL; // In order to use automatic allocation
	if (options & CAM_APPROX_HESSIAN) {
	    camFastApproxHessianDetectorFixedScale(&integral, &results[scale], scale);
	} else {
	    camFastHessianDetectorFixedScale(&integral, &results[scale], scale);
	}
	//sprintf(str, "output/features_fast_hessian_%02d.pgm", scale);
	//camSavePGM(&results[scale], str);
	camFindLocalMaximaCircle5(&results[scale], &keypoints[nb_keypoints], &j);
	for (i = nb_keypoints; i < nb_keypoints + j; i++) {
#ifdef POST_SCALING
	    keypoints[i].value = ((keypoints[i].value * coeff[scale]) >> 6);
#endif
	    keypoints[i].scale = scale;
	    keypoints[i].x <<= CamSampling[scale];
	    keypoints[i].y <<= CamSampling[scale];
	    keypoints[i].x += CamSamplingOffset[scale];
	    keypoints[i].y += CamSamplingOffset[scale];
	}
	nb_keypoints += j;
    }	

    // Spatial clustering

    // 1st pass : Count the number of points per cluster
    widthClusters = (((width - 1) >> 5) + 1);
    heightClusters = (((height - 1) >> 5) + 1); 
    nbClusters = widthClusters * heightClusters;
    nbPerCluster = (int*)malloc(nbClusters * sizeof(int));
    for (i = 0; i < nbClusters; i++) nbPerCluster[i] = 0;
    for (i = 0; i < nb_keypoints; i++) {
	cluster = (keypoints[i].y >> 5) * widthClusters + (keypoints[i].x >> 5);
	nbPerCluster[cluster]++;
    }
    firstPoint =  (CamKeypointShort***)malloc(nbClusters * sizeof(CamKeypointShort**));    
    // Cluster ordered pointers to the points in clusters
    ptrPoints = (CamKeypointShort**)malloc(nb_keypoints * sizeof(CamKeypointShort*));
    firstPoint[0] = ptrPoints;
    for (i = 1; i < nbClusters; i++) firstPoint[i] = firstPoint[i - 1] + nbPerCluster[i - 1]; 
    // 2nd pass : Start again and register pointers
    for (i = 0; i < nb_keypoints; i++) {
	cluster = (keypoints[i].y >> 5) * widthClusters + (keypoints[i].x >> 5);
	*firstPoint[cluster] = &keypoints[i];
	firstPoint[cluster]++;	
    }
    firstPoint[0] = ptrPoints;
    for (i = 1; i < nbClusters; i++) firstPoint[i] = firstPoint[i - 1] + nbPerCluster[i - 1]; 
    // OK. clustering is finished.
    
    // Now, scan all the clusters and the points inside
    // for bad points removal
    for (cluster = 0; cluster < nbClusters; cluster++) {
	x = cluster % widthClusters;
	y = cluster / widthClusters;
	for (i = 0; i < nbPerCluster[cluster]; i++) {
	    point1 = *(firstPoint[cluster] + i);
	    
	    // Build the set of points to compare with
	    nbClusters2scan = 0;
	    if (i < nbPerCluster[cluster] - 1) {
		first2scan[0] = firstPoint[cluster] + i + 1;
		nbPoints2scan[0] = nbPerCluster[cluster] - i - 1; 
		nbClusters2scan = 1;
	    }
	    if (x != 0 && y != heightClusters - 1) {
		first2scan[nbClusters2scan] = firstPoint[cluster + widthClusters - 1];
		nbPoints2scan[nbClusters2scan] = nbPerCluster[cluster + widthClusters - 1];
		nbClusters2scan++;
	    }
	    if (y != heightClusters - 1) {
		first2scan[nbClusters2scan] = firstPoint[cluster + widthClusters];
		nbPoints2scan[nbClusters2scan] = nbPerCluster[cluster + widthClusters];
		nbClusters2scan++;
	    }
	    if (x != widthClusters - 1 && y != heightClusters - 1) {
		first2scan[nbClusters2scan] = firstPoint[cluster + widthClusters + 1];
		nbPoints2scan[nbClusters2scan] = nbPerCluster[cluster + widthClusters + 1];
		nbClusters2scan++;
	    }
	    if (x != widthClusters - 1) {
		first2scan[nbClusters2scan] = firstPoint[cluster + 1];
		nbPoints2scan[nbClusters2scan] = nbPerCluster[cluster + 1];
		nbClusters2scan++;
	    }
	    
	    // Now, scan all these points and compare with point1
	    for (j = 0; j < nbClusters2scan; j++) {
		for (k = 0; k < nbPoints2scan[j]; k++) {
		    point2 = *(first2scan[j] + k);
		    // Now I have point1 and point2
		    // They are next to each other...
		    // Let's compare them and possibly keep only one of them

		    // Which one is the bigger ?
		    if (point1->scale > point2->scale) {
			bigOne = point1;
			smallOne = point2;
		    } else {
			bigOne = point2;
			smallOne = point1;
		    }
			
		    // First of all, are these feature points similar in size ?
		    if (bigOne->scale - smallOne->scale != 1) continue;

		    // These two points are interesting
		    d = (point1->x - point2->x) * (point1->x - point2->x) + (point1->y - point2->y) * (point1->y - point2->y);
		    thr = (CamScale[bigOne->scale] - CamScale[smallOne->scale]);
		    thr *= thr;
		    thr <<= 2;
		    if (d < thr) {
			// These are next to each other
			// Are they of similar value ?
			if (abs(point1->value - point2->value) < ((point1->value + point2->value) >> 4)) {
			    // Yes they are. Keep the bigger.
			    smallOne->value = -1;
			} else {
			    // No they are not. Keep the better.
			    if (point1->value > point2->value) {
				point2->value = -1;
			    } else {
				point1->value = -1;
			    }
			}
		    }
		}
	    }	
	}
    }
   
    free(nbPerCluster);
    free(ptrPoints);
    free(firstPoint);

    // Remove all the points that have been marked to be destroyed
    pnb_keypoints = 0;
    i = 0;
    while (i < nb_keypoints && keypoints[i].value != -1) {
	i++;
	pnb_keypoints++;
    }
    for (; i < nb_keypoints; i++) {
	if (keypoints[i].value != -1) {
	    keypoints[pnb_keypoints++] = keypoints[i];
	}
    }
    nb_keypoints = pnb_keypoints;
   
    if (options & CAM_NO_INTERPOLATION) {
	for (i = 0; i < nb_keypoints; i++) {
	    keypoints[i].scale = CamScale[keypoints[i].scale] << 2;
	}
    } else {
    /* Interpolation :
     * Maxima : solve([y1=a*(x1-p)^2+b,y2=a*(-p)^2+b,y3=a*(x3-p)^2+b],[a,b,p])
     * Maxima yields the following formula for parabolic interpolation
                                       2         2     2         2
			             x1  y3 + (x3  - x1 ) y2 - x3  y1
			       p = ------------------------------------
                                   2 x1 y3 + (2 x3 - 2 x1) y2 - 2 x3 y1
    
     */
	for (i = 0; i < nb_keypoints; i++) {
	    if (keypoints[i].scale == 0 || keypoints[i].scale == nbScales - 1) {
		// Remove this point : its scale can't be interpolated
		keypoints[i].value = -1;
	    } else {
		x1 = CamScale[keypoints[i].scale - 1] - CamScale[keypoints[i].scale];
		x3 = CamScale[keypoints[i].scale + 1] - CamScale[keypoints[i].scale];

		scale = keypoints[i].scale - 1;
		ptr = (CAM_UINT16*)(results[scale].imageData + results[scale].widthStep * (keypoints[i].y >> CamSampling[scale])) + (keypoints[i].x >> CamSampling[scale]);
		y1 = *ptr;
		// Scale the point
#ifdef POST_SCALING
		y1 = ((y1 * coeff[scale]) >> 6);
#endif
		scale = keypoints[i].scale;
		ptr = (CAM_UINT16*)(results[scale].imageData + results[scale].widthStep * (keypoints[i].y >> CamSampling[scale])) + (keypoints[i].x >> CamSampling[scale]);
		y2 = *ptr;
		// Scale the point
#ifdef POST_SCALING
		y2 = ((y2 * coeff[scale]) >> 6);
#endif

		scale = keypoints[i].scale + 1;
		ptr = (CAM_UINT16*)(results[scale].imageData + results[scale].widthStep * (keypoints[i].y >> CamSampling[scale])) + (keypoints[i].x >> CamSampling[scale]);
		y3 = *ptr;
		// Scale the point
#ifdef POST_SCALING
		y3 = ((y3 * coeff[scale]) >> 6);
#endif

		if (y3 > y2) {
		    keypoints[i].scale++;
		    keypoints[i].value = -1; // Destroy
		    //i--;
		    continue;
		}
		if (y1 > y2) {
		    keypoints[i].scale--;
		    keypoints[i].value = -1; // Destroy
		    //i--;
		    continue;
		}
		num = (x1 * x1 * y3 + (x3 * x3 - x1 * x1) * y2 - x3 * x3 * y1);
		den = x1 * y3 + (x3 - x1) * y2 - x3 * y1;
		if (den == 0)
		    keypoints[i].value = -1; // Destroy
		else {
		    p = (num << 1) / den; // shift only by 1, because den is shifted by 2	    
		    keypoints[i].scale = p + (CamScale[keypoints[i].scale] << 2);
		}
	    }
	}
    }

    // Memory deallocation 
    for (scale = 0; scale < nbScales; scale++) {
	camDeallocateImage(&results[scale]);
    }
    camDeallocateImage(&integral);

    // Remove again all the points that have been marked to be destroyed
    pnb_keypoints = 0;
    i = 0;
    while (i < nb_keypoints && keypoints[i].value != -1) {
	i++;
	pnb_keypoints++;
    }
    for (; i < nb_keypoints; i++) {
	if (keypoints[i].value != -1) {
	    keypoints[pnb_keypoints++] = keypoints[i];
	}
    }
    nb_keypoints = pnb_keypoints;
 
    // Sort the features according to value
    qsort(keypoints, nb_keypoints, sizeof(CamKeypointShort), camSortKeypointsShort);
    
    // Angle
    if (options & CAM_UPRIGHT) {
	for (i = 0; i < nb_keypoints; i++) {
	    keypoints[i].angle = 0;
	}
    } else {
        camAllocateImage(&filter, CAM_ORIENTATION_STAMP_SIZE, CAM_ORIENTATION_STAMP_SIZE, CAM_DEPTH_16S);
	camBuildGaussianFilter(&filter, camSigmaParam);
	pnb_keypoints = nb_keypoints;
	if (pnb_keypoints > nb_max_keypoints) pnb_keypoints = nb_max_keypoints;
	for (i = 0; i < pnb_keypoints; i++) {
	    nb_keypoints += camKeypointOrientation(source, &keypoints[i], &filter, &keypoints[nb_keypoints]);
	}
	camDeallocateImage(&filter);
    }
    
    // Sort again the features according to value
    qsort(keypoints, nb_keypoints, sizeof(CamKeypointShort), camSortKeypointsShort);
    
    // Keypoints allocation
    pnb_keypoints = nb_keypoints;
    if (pnb_keypoints > nb_max_keypoints) pnb_keypoints = nb_max_keypoints;
    if (points->allocated == 0) {
	camAllocateKeypoints(points, pnb_keypoints);
    } else if (points->allocated < pnb_keypoints) {
	camFreeKeypoints(points);
	camAllocateKeypoints(points, pnb_keypoints);
    }
    if (points->bag == NULL) {
#ifdef __SSE2__
	points->bag = (CamKeypoint*)_mm_malloc(sizeof(CamKeypoint) * pnb_keypoints, 16);
#else
	points->bag = (CamKeypoint*)malloc(sizeof(CamKeypoint) * pnb_keypoints);
#endif
    }
    for (i = 0; i < pnb_keypoints; i++) {
        points->keypoint[i] = &points->bag[i];
	points->keypoint[i]->x = keypoints[i].x;
	points->keypoint[i]->y = keypoints[i].y;
	points->keypoint[i]->scale = keypoints[i].scale;
	points->keypoint[i]->value = keypoints[i].value;
	points->keypoint[i]->angle = keypoints[i].angle;
    }
    points->nbPoints = pnb_keypoints;
    free(keypoints);

    camKeypointsDescriptor(points, source, options);

    // Finally, set the points' set 
    for (i = 0; i < points->nbPoints; i++) {
	points->keypoint[i]->set = points;
    }

    camInternalROIPolicyExit(&iROI);
    return 1;
}

int camDrawKeypoints(CamKeypoints *points, CamImage *dest, int color)
{
    int i;
    for (i = 0; i < points->nbPoints; i++) {
	camDrawKeypoint(points->keypoint[i], dest, color);
    }
    return 1;
}

int camDrawKeypoint(CamKeypoint *point, CamImage *dest, int color)
{
    int x, y;
    double costheta, sintheta;

    if (dest->roi) {
	x = dest->roi->xOffset;
	y = dest->roi->yOffset;
    } else {
	x = 0;
	y = 0;
    }

    camDrawCircle(dest, x + point->x, y + point->y, (point->scale >> 2), color);
    if (point->angle != 0) {
	costheta = cos(point->angle * 2 * M_PI / 360);
	sintheta = sin(point->angle * 2 * M_PI / 360);
	camDrawLine(dest, (int)floor(x + point->x + costheta * (point->scale >> 2) + 0.5), (int)floor(y + point->y - sintheta * (point->scale >> 2) + 0.5), x + point->x, y + point->y, color);
    }
    return 1;
}

