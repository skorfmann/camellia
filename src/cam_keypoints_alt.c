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
#define CAM_MAX_KEYPOINTS 100000
#define CAM_MAX_SEEDS 20000
#define CAM_ORIENTATION_STAMP_SIZE 30

typedef struct {
    CamImage *integral;
    int count;
    int min_distance_to_border[CAM_MAX_SCALE];
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
    data->count = 0;
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
	data->min_distance_to_border[scale] = sizep + scale + 1;

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
	//printf("scale = %d sizep = %d widthp = %d preshift = %d coeff = %d\n", scale, sizep, widthp, data->preshift[scale], data->coeff[scale]); 
    }
}

// All 2D integrals are proportional to sizep * sizep
// so the determinant is proportional to sizep^4.
// For scale = 100, sizep = 201 and sizep^4 = 1.6b (compared to 81 width scale = 1). We need prescaling before multiply (preshifting indeed).

typedef struct {
    int c[3];   
} CamKeypointLocation;

// This code involves (exscluding boundaries checking) 8 * 4 = 32 memory accesses to the integral image
// 4 right shifts (1 with constant), and 5 32-bits multiplications (not overflowing - no need for 64 bits result)
// (and approximately 8 * (4 + 3) + 2 * 3 + 3 + 1 = 66 additions) 
int camHessianEstimate(CamHessianEstimateData *data, CamKeypointLocation *keypoint)
{
    int det, Dxx, Dyy, Dxy, tmp;
    int x = keypoint->c[0], y = keypoint->c[1], scale = keypoint->c[2];
    int *offset_ptr = data->offset[scale];
    unsigned long *imptr;

    data->count++;
    imptr = ((unsigned long*)(data->integral->imageData + y * data->integral->widthStep)) + x;

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
    // det certainly stands within 26 bits 
    det *= data->coeff[scale];
    // it should stay within 32 bits range
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
    int c, i, x, y, ok, tmp;
    int nb_keypoints = 0, pnb_keypoints;
    int nb_seeds = 0, nb_tests, found_better, move, counter_not_found;
    CamKeypointShort *keypoints;
    CamKeypointLocation *seeds, neighbour, best_keypoint, current_keypoint, gradient, min_gradient;
    int best_keypoint_value, current_keypoint_value, neighbour_value;
    int stat_upscales, stat_moves, stat_good;

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

    camHessianEstimateDataBuild(&data, &integral);

    // Allocate temp memory for keypoints and seeds
    keypoints = (CamKeypointShort*)malloc(CAM_MAX_KEYPOINTS * sizeof(CamKeypointShort));

    // Ready to go
    // Seeding
    // Should take ROI into consideration
    seeds = (CamKeypointLocation*)malloc(CAM_MAX_SEEDS * sizeof(CamKeypointLocation));
#define INTERVAL 20
    for (y = INTERVAL * 2; y < height - INTERVAL * 2; y += INTERVAL) {
	for (x = INTERVAL * 2; x < width - INTERVAL * 2; x += INTERVAL) {
	    seeds[nb_seeds].c[0] = x + iROI.srcroi.xOffset;
	    seeds[nb_seeds].c[1] = y + iROI.srcroi.yOffset;
	    seeds[nb_seeds].c[2] = INTERVAL;
	    nb_seeds++;
	}
    }	

    // Go !
#define TEST_NEIGHBOUR \
    if ((neighbour.c[2] < 1) || (neighbour.c[2] >= CAM_MAX_SCALE)) break; \
    tmp = data.min_distance_to_border[neighbour.c[2]]; \
    if ((neighbour.c[0] < tmp) || (neighbour.c[0] > integral.width - tmp) || (neighbour.c[1] < tmp) || (neighbour.c[1] > integral.height - tmp)) break; \
    neighbour_value = camHessianEstimate(&data, &neighbour); \
    nb_tests++; \
    if (neighbour_value > best_keypoint_value) { \
	found_better = 1; counter_not_found = 0; best_keypoint = neighbour; best_keypoint_value = neighbour_value; stat_good++; \
    } 
    
    stat_upscales = 0; stat_moves = 0; stat_good = 0;
    for (c = 0; c < nb_seeds; c++) {
#define MAX_NOT_FOUND 10
#define INV_GAIN_MOVE 2
#define INV_GAIN_GRADIENT 24
	current_keypoint = seeds[c];
        current_keypoint_value = best_keypoint_value = camHessianEstimate(&data, &current_keypoint);
	nb_tests = 1;
	best_keypoint = current_keypoint;
	counter_not_found = 0;
	do {
	    found_better = 0;
	    move = current_keypoint.c[2] >> INV_GAIN_MOVE;
	    if (move == 0) move = 1;
	    // Compute gradient
	    for (i = 0; i != 3; i++) {
		neighbour = current_keypoint;
		neighbour.c[i] += move;
	        TEST_NEIGHBOUR;	
		gradient.c[i] = neighbour_value - current_keypoint_value;
	    }
	    // Search at gradient ascent location
	    move = 0;
	    for (i = 1; i != 3; i++) {
		if (abs(gradient.c[i]) > abs(gradient.c[move])) move = i;
	    }
	    for (i = 0; i != 3; i++) min_gradient.c[i] = 0;
	    min_gradient.c[move] = (gradient.c[move] > 0) ? 1 : -1;
	    for (ok = 0, i = 0; i != 3; i++) {
		gradient.c[i] = (gradient.c[i] * current_keypoint.c[2]) >> INV_GAIN_GRADIENT;
		if (gradient.c[i] != 0) ok = 1;
	    }
	    if (!ok) gradient = min_gradient;
	    neighbour = current_keypoint;
	    for (i = 0; i != 3; i++) neighbour.c[i] += gradient.c[i];
	    TEST_NEIGHBOUR;
	    if (!found_better) {
		counter_not_found++;
		if (counter_not_found >= MAX_NOT_FOUND) {
		    // Try upscaling
		    stat_upscales++;
		    neighbour = best_keypoint;
		    neighbour.c[2] *= 2;
		    TEST_NEIGHBOUR;
		    current_keypoint = neighbour;		    
		} else {
		    // Take the direction of gradient, even if it's not better
		    current_keypoint = neighbour;
		    current_keypoint_value = neighbour_value;
		    stat_moves++;
		}
	    } else {
		// OK. Good ! Let's move on !
		current_keypoint = best_keypoint;
		current_keypoint_value = best_keypoint_value;
	    }
	} while (1);
	
	// Let's record this keypoint
	if (best_keypoint_value > 0) {
	    keypoints[nb_keypoints].x = best_keypoint.c[0];
	    keypoints[nb_keypoints].y = best_keypoint.c[1];
	    keypoints[nb_keypoints].scale = (best_keypoint.c[2] * 2 + 1) * 4;
	    keypoints[nb_keypoints].value = best_keypoint_value;
	    nb_keypoints++;
	}
    }

    free(seeds);
    
    printf("Good = %d, Moves = %d, Upscales = %d, ", stat_good, stat_moves, stat_upscales);
    printf("Count = %d\n", data.count);

    // Integral Image is now useless
    camDeallocateImage(&integral);

    // Sort the features according to value
    qsort(keypoints, nb_keypoints, sizeof(CamKeypointShort), camSortKeypointsShort);

    // Remove duplicates

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
    CamKeypointLocation key;

    const int xp[4] = {-1, 1, 1, -1};
    const int yp[4] = {-1, -1, 1, 1};
    CamWarpingParams params;

    printf("Alternative keypoints detection :\n");
    camAllocateImage(&image, 256, 256, CAM_DEPTH_8U);
    camAllocateKeypoints(&points, 100);

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
    dest.imageData = NULL; 

    integral.imageData = NULL; // in order to use automatic allocation
    camIntegralImage(&image, &integral);
    camHessianEstimateDataBuild(&data, &integral);
    camDeallocateImage(&integral);
    key.c[0] = 50; key.c[1] = 50; key.c[2] = 6;
    printf("Value = %d\n", camHessianEstimate(&data, &key));

    /*
    camKeypointsDetector(&image, &points, 100, 0);
    for (i = 0; i < points.nbPoints; i++) {
	printf("x=%d y=%d value=%d scale=%d size=%d angle=%d\n", points.keypoint[i]->x, points.keypoint[i]->y, points.keypoint[i]->value, points.keypoint[i]->scale, points.keypoint[i]->size, points.keypoint[i]->angle);
    }
    camDrawKeypoints(&points, &image, 128);
    */camSavePGM(&image, "output/features_reference.pgm");
    
    camDeallocateImage(&image);
    camDeallocateImage(&dest);
    camFreeKeypoints(&points);
}

