/* CamKeypoints alternative implementation using gradient search
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

#define CAM_MAX_SCALE 100
#define CAM_MAX_KEYPOINTS 10000
#define CAM_ORIENTATION_STAMP_SIZE 30

typedef struct {
    CamImage *integral;
    int offset[CAM_MAX_SCALE][8 * 4];
    int preshift[CAM_MAX_SCALE];
    int coeff[CAM_MAX_SCALE];
} CamHessianEstimateData;

void camHessianEstimateDataBuild(CamHessianEstimateData *data, CamImage *integral)
{
    int sizep, widthp, linc;
    int left, top, right, bottom;
    double v;
    int scale;

    data->integral = integral;
    linc = integral->widthStep / 4;

    for (scale = 0; scale < CAM_MAX_SCALE; scale++) {
	
	// Size of patch is odd and linear wrt scale : sizep = 1 + scale * 2
	// width of rectangular patch is the odd value next to widthp = 5 * size / 3
	sizep = 1 + scale * 2;
	v = (5.0 * sizep) / 3.0;
	widthp = (int)v;
	if (!(widthp & 1)) { // if is even
	    widthp++;
	}

// 0 : top-left, 1 : top-right, 2 : bottom-left, 3 : bottom-right
#define SET_OFFSET(i) \
	data->offset[scale][i * 4] = left - 1 + (top - 1) * linc; \
	data->offset[scale][i * 4 + 1] = right + (top - 1) * linc; \
	data->offset[scale][i * 4 + 2] = left - 1 + bottom * linc; \
	data->offset[scale][i * 4 + 3] = right + bottom * linc; 

	left = -widthp / 2; top = -sizep / 2; right = -left; bottom = -top;
	SET_OFFSET(0);
	top = -sizep * 3/2; bottom = -top;
	SET_OFFSET(1);
	left = -sizep / 2; top = -widthp / 2; right = -left; bottom = -top;
	SET_OFFSET(2);
	left = -sizep * 3/2; right = -left;
	SET_OFFSET(3);
	// totalwidth = sizep * 3
	// occupied space = sizep * 2
	// free space = sizep
	// free space on left = sizep / 4
	left = -sizep * 3/2 + sizep / 4; top = left; right = left + sizep - 1; bottom = top + sizep - 1;
	SET_OFFSET(4);
	left = sizep / 4; top = sizep / 4; right = left + sizep - 1; bottom = top + sizep - 1;
	SET_OFFSET(5);
	left = -sizep * 3/2 + sizep / 4; top = sizep / 4; right = left + sizep - 1; bottom = top + sizep - 1;
	SET_OFFSET(6);
	left = sizep / 4; top = -sizep * 3/2 + sizep / 4; right = left + sizep - 1; bottom = top + sizep - 1;
	SET_OFFSET(7);

	data->preshift[scale] = (int)(2 * log(sizep) / log(2));
	data->preshift[scale] -= 5; // So that when it should have shifted back to 8 bits, it shifts back only to 13...
	if (data->preshift[scale] < 1) data->preshift[scale] = 1;

	// Coeff : for sizep, divider should be sizep^4
	// With data in 8 bits range, 2D integral is in (8 bits * sizep) ^ 2
	// In order to stay in 8 bits range, one should multiply by 1 / (8 * sizep^4)
	// which is equivalent to (2^16 / (sizep ^ 4)) >> 24
	// We need to stay in 32 bits range, so we can multiply by (2^16 / (sizep ^ 4)) only
	// This is without considering preshifting
	// determinant will have already been shifted by data->preshift[scale] * 2 bits
	// that is divided by 2^(data->preshift[scale] * 2). 
	// The multiplier should be augmented accordingly
	data->coeff[scale] = (int)floor(pow(2, 16) / pow(sizep, 4) * pow(2, data->preshift[scale] * 2) + 0.5);
	// All data->coeff[scale] should be in the range from 1 to 2^6 = 64, turning a 26 bits max determinant into a 32 bits max determinant
	printf("scale = %d sizep = %d widthp = %d preshift = %d coeff = %d\n", scale, sizep, widthp, data->preshift[scale], data->coeff[scale]); 
    }
}

// All 2D integrals are proportional to sizep * sizep
// so the determinant is proportional to sizep^4.
// For scale = 100, sizep = 201 and sizep^4 = 1.6b (compared to 81 width scale = 1). We need prescaling before multiply (preshifting indeed).

// This code involves 8 * 4 = 32 memory accesses to the integral image
// 4 right shifts (1 with constant), and 5 32-bits multiplications (not overflowing - no need for 64 bits result)
// (and approximately 8 * (4 + 3) + 2 * 3 + 3 + 1 = 66 additions) 
int camHessianEstimate(CamHessianEstimateData *data, int x, int y, int scale)
{
    int det, Dxx, Dyy, Dxy, tmp;
    unsigned long *imptr = ((unsigned long*)(data->integral + y * data->integral->widthStep)) + x;
    int *offset_ptr = data->offset[scale];

#define CAM_INTEGRAL(i) \
    ( *(imptr + offset_ptr[i * 4]) - *(imptr + offset_ptr[i * 4 + 1]) - *(imptr + offset_ptr[i * 4 + 2]) + *(imptr + offset_ptr[i * 4 + 3]) )
		    
    tmp = CAM_INTEGRAL(0);
    Dxx = CAM_INTEGRAL(1) - (tmp + tmp + tmp);
    Dxx >>= data->preshift[scale]; 
    tmp = CAM_INTEGRAL(2); 
    Dyy = CAM_INTEGRAL(3) - (tmp + tmp + tmp);
    Dyy >>= data->preshift[scale];
    Dxy = CAM_INTEGRAL(4) + CAM_INTEGRAL(5) - CAM_INTEGRAL(6) - CAM_INTEGRAL(7);
    Dxy >>= data->preshift[scale];
    // With this preshifting, Dxx, Dyy and Dxy are 13 bits wide max
    // det = Dxx * Dyy - 0.81 * Dxy^2
    det = Dxx * Dyy - ((13 * Dxy * Dxy) >> 4);
    // det certainly stands in 26 bits 
    det *= data->coeff[scale];
    // it should stay in 32 bits range
    return det;
}

extern int camSigmaParam;
int camKeypointOrientation(CamImage *source, CamKeypointShort *point, CamImage *filter, CamKeypointShort *next_point);
int camSortKeypointsShort(const void *p1x, const void *p2x);
int camBuildGaussianFilter(CamImage *image, double sigma);
void camKeypointsInternalsPrepareDescriptor();
int camKeypointsDescriptor(CamImage *source, CamKeypoint *point, CamImage *filter, int option);

int camKeypointsDetector(CamImage *source, CamKeypoints *points, int nb_max_keypoints, int options)
{
    CamImage integral, filter;
    CamInternalROIPolicyStruct iROI;
    CamHessianEstimateData data;
    CamROI *roi, roix;
    int width, height;
    int i;
    int nb_keypoints = 0, pnb_keypoints;
    CamKeypointShort *keypoints;

    // Parameters checking
    CAM_CHECK(camKeypointsDetector, camInternalROIPolicy(source, NULL, &iROI, 1));
    CAM_CHECK_ARGS(camKeypointsDetector, (source->depth & CAM_DEPTH_MASK) >= 8);
    CAM_CHECK_ARGS(camKeypointsDetector, points->allocated >= 0);
    CAM_CHECK_ARGS(camKeypointsDetector, source->nChannels == 1 || ((source->nChannels == 3) && (source->dataOrder == CAM_DATA_ORDER_PLANE)));
    width = iROI.srcroi.width;
    height = iROI.srcroi.height;
    points->width = width;
    points->height = height;

    // Compute integral image
    integral.imageData = NULL;
    roi = source->roi;
    camSetMaxROI(&roix, source);
    roix.coi = 1;
    source->roi = &roix;
    camIntegralImage(source, &integral);
    points->nbPoints = 0;
    integral.roi = &iROI.srcroi;
    source->roi = roi;

    // Allocate temp memory for keypoints
    keypoints = (CamKeypointShort*)malloc(CAM_MAX_KEYPOINTS * sizeof(CamKeypointShort));

    // Integral Image is now useless
    camDeallocateImage(&integral);

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

    camAllocateImage(&filter, 20, 20, CAM_DEPTH_16S);
    camBuildGaussianFilter(&filter, camSigmaParam);
    
    camKeypointsInternalsPrepareDescriptor();
    // Get the signature from the selected feature points
    for (i = 0; i < points->nbPoints; i++) {
	camKeypointsDescriptor(source, points->keypoint[i], &filter, options);
    }

    // Finally, set the points' set 
    for (i = 0; i < points->nbPoints; i++) {
	points->keypoint[i]->set = points;
    }

    camDeallocateImage(&filter);
    camInternalROIPolicyExit(&iROI);
    return 1;
}

void test_camKeypointsAlt()
{
    CamHessianEstimateData data;
    CamImage image, integral, dest;
    CamKeypoints points;
    int i;
    int c, t1, t2;
    const int x = 8;
    char str[256];
    FILE *handle;
    int angle;
    double costheta;
    double sintheta;
    CamROI roi;

    const int xp[4] = {-1, 1, 1, -1};
    const int yp[4] = {-1, -1, 1, 1};
    CamWarpingParams params;

    printf("Alternative keypoints detection :\n");
    camAllocateImage(&image, 256, 256, CAM_DEPTH_8U);
    camAllocateKeypoints(&points, 1000);

    camSet(&image, 0);
    camDrawRectangle(&image, 102, 120, 156, 152, 255);
    camFillColor(&image, 103, 121, 255, -1);
    camDrawRectangle(&image, 122, 35, 176, 67, 128);
    camFillColor(&image, 123, 36, 128, -1);

    camDrawCircle(&image, 50, 50, 10, 255);
    camFillColor(&image, 50, 50, 255, -1);
    camDrawCircle(&image, 80, 50, 5, 255);
    camFillColor(&image, 80, 50, 255, -1);
    camDrawCircle(&image, 100, 50, 3, 255);
    camFillColor(&image, 100, 50, 255, -1);

#if 1
    angle = 20;
    costheta = cos(angle * 2 * M_PI / 360);
    sintheta = sin(angle * 2 * M_PI / 360);
    for (i = 0; i < 4; i++) {
	params.p[i].x = (int)floor((costheta * xp[i] - sintheta * yp[i]) * 15 + 0.5);
	params.p[i].y = (int)floor((sintheta * xp[i] + costheta * yp[i]) * 15 + 0.5);
	params.p[i].x += 192;
	params.p[i].y += 192;
    }
    for (i = 0; i < 4; i++) {
	camDrawLine(&image, params.p[i].x, params.p[i].y, params.p[(i+1)%4].x, params.p[(i+1)%4].y, 255);
    }
    camFillColor(&image, 192, 192, 255, -1);
 
    angle = 30;
    costheta = cos(angle * 2 * M_PI / 360);
    sintheta = sin(angle * 2 * M_PI / 360);
    for (i = 0; i < 4; i++) {
	params.p[i].x = (int)floor((costheta * xp[i] - sintheta * yp[i]) * 10 + 0.5);
	params.p[i].y = (int)floor((sintheta * xp[i] + costheta * yp[i]) * 10 + 0.5);
	params.p[i].x += 50;
	params.p[i].y += 192;
    }

    for (i = 0; i < 4; i++) {
	camDrawLine(&image, params.p[i].x, params.p[i].y, params.p[(i+1)%4].x, params.p[(i+1)%4].y, 255);
    }
    camFillColor(&image, 50, 192, 255, -1);

#endif
    integral.imageData = NULL; /* in order to use automatic allocation */
    dest.imageData = NULL; 
	
    camIntegralImage(&image, &integral);
    camHessianEstimateDataBuild(&data, &integral);

    camDeallocateImage(&image);
    camDeallocateImage(&integral);
    camDeallocateImage(&dest);
    camFreeKeypoints(&points);
}


