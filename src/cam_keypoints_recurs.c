/* CamKeypoints alternative implementation using recursive search 
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

#define CAM_MAX_KEYPOINTS 100000
#define CAM_ORIENTATION_STAMP_SIZE 30
#define CAM_MIN_SCALE 3
#define CAM_SCALE_MARGIN 0

extern int camPatchSizeParam;
extern double camSigmaParam;
int camKeypointOrientation(CamImage *source, CamKeypointShort *point, CamImage *filter, CamKeypointShort *next_point);
int camSortKeypointsShort(const void *p1x, const void *p2x);
int camBuildGaussianFilter(CamImage *image, double sigma);

int camKeypointsRecursiveDetector(CamImage *source, CamKeypoints *points, int nb_max_keypoints, int options)
{
    CamImage integral, filter;
    CamInternalROIPolicyStruct iROI;
    CamROI *roi, roix;
    int i, width, height;
    CamKeypointShort *keypoints;
    int pnb_keypoints;
    int nb_keypoints = 0;
    unsigned int *abs_value_lines[2];
    int *value_lines[2];
    unsigned char *scale_lines[2];
    unsigned char *lmax_lines[2];
    unsigned int *lmax[2], nblmax[2];
    int widthStep;

    camPatchSizeParam = 16; //32 * 2 / 3; // Equivalent optimal value wrt SURF (the descriptor can partially lie outside screen)
    
    // Parameters checking
    CAM_CHECK(camKeypointsRecursiveDetector, camInternalROIPolicy(source, NULL, &iROI, 1));
    CAM_CHECK_ARGS(camKeypointsRecursiveDetector, (source->depth & CAM_DEPTH_MASK) >= 8);
    CAM_CHECK_ARGS(camKeypointsRecursiveDetector, points->allocated >= 0);
    CAM_CHECK_ARGS(camKeypointsRecursiveDetector, source->nChannels == 1 || ((source->nChannels == 3) && (source->dataOrder == CAM_DATA_ORDER_PLANE)));
    width = iROI.srcroi.width;
    height = iROI.srcroi.height;

    // Compute integral image
    integral.imageData = NULL;
    roi = source->roi;
    camSetMaxROI(&roix, source);
    roix.coi = 1;
    source->roi = &roix;
    camIntegralImage(source, &integral);
    integral.roi = &iROI.srcroi;
    source->roi = roi;
    widthStep = integral.widthStep / 4;

    // Allocate temp memory for keypoints
    keypoints = (CamKeypointShort*)malloc(CAM_MAX_KEYPOINTS * sizeof(CamKeypointShort));
    points->nbPoints = 0;
    // Allocate value and scale lines
    for (i = 0; i < 2; i++) {
        scale_lines[i] = (unsigned char*)malloc(width * sizeof(unsigned char));
        value_lines[i] = (int*)malloc(width * sizeof(int));
        abs_value_lines[i] = (unsigned int*)malloc(width * sizeof(unsigned int));
        lmax_lines[i] = (unsigned char*)malloc(width * sizeof(unsigned char));
        lmax[i] = (unsigned int*)malloc(width * sizeof(unsigned int));
        nblmax[i] = 0;
    }

    // Go !
    {
        int y;
        int max_scale_y;
        unsigned int *abs_current_value_line, *abs_previous_value_line;
        int *current_value_line, *previous_value_line;
        unsigned char *current_scale_line, *previous_scale_line;
        unsigned char *current_lmax_line, *previous_lmax_line;
        unsigned int *current_lmax, *previous_lmax;
        unsigned int current_nblmax, previous_nblmax;

        current_value_line = value_lines[0]; 
        previous_value_line = value_lines[1];
        abs_current_value_line = abs_value_lines[0]; 
        abs_previous_value_line = abs_value_lines[1];
        current_scale_line = scale_lines[0]; 
        previous_scale_line = scale_lines[1];
        current_lmax_line = lmax_lines[0]; 
        previous_lmax_line = lmax_lines[1];
        current_lmax = lmax[0]; 
        previous_lmax = lmax[1];
        // Fill the previous value line with a dummy value
        for (i = 0; i < width; i++) {
            abs_previous_value_line[i] = 0;
            previous_value_line[i] = 0;
            previous_scale_line[i] = 1;
            previous_lmax_line[i] = 0;
        }
        current_nblmax = 0;
        previous_nblmax = 0;

        for (y = 0; y < height; y++) {
            int x, v;
            // max_scale is the biggest scale we may apply. It must be more than 0, otherwise we can't do any image processing
            max_scale_y = (y + iROI.srcroi.yOffset) >> 1;
            // The preceding formula comes from :
            // - When iROI.yOffset is 0, max_scale should be 1 when y is 2 (1 pixel margin due to integral image)
            // - When iROI.yOffset is 0, max_scale should be 1 when y is 3 (1 pixel margin due to integral image)
            // - When iROI.yOffset is 0, max_scale should be 2 when y is 4 (1 pixel margin due to integral image)
            // Check max_scale with the bottom line
            v = (source->height - 1 - (y + iROI.srcroi.yOffset)) >> 1;
            // The preceding formula comes from "When iROI.yOffset is 0, max_scale should be 1 when y is equal to source->height - 3
            if (v < max_scale_y) max_scale_y = v;

            if (max_scale_y < 1) {
                // We can't compute anything for that line
                for (i = 0; i < width; i++) {
                    current_value_line[i] = 0;
                    current_scale_line[i] = 1;
                }
            } else {
                int scale_up = 1;
                int current_value = 0;
                unsigned int abs_current_value = 0;
                int current_scale = 1;
                unsigned int *ptr = ((unsigned int*)(integral.imageData + (y + iROI.srcroi.yOffset) * integral.widthStep)) + iROI.srcroi.xOffset;
                int s1, s2;

                for (x = 0; x < width; x++, ptr++) {
                    int first_scale = current_scale;
                    int max_scale = max_scale_y;
                    // Check max_scale with the left side of the frame
                    v = (x + iROI.srcroi.xOffset) >> 1;
                    if (v < max_scale) max_scale = v;
                    // Check max_scale with the right side of the frame
                    v = (source->width -1 - (x + iROI.srcroi.xOffset)) >> 1;
                    if (v < max_scale) max_scale = v;
                    // We need a scale margin in order to compute the descriptor
                    max_scale -= CAM_SCALE_MARGIN;

                    // Choose the first scale to evaluate, depending on upper and left results
                    if (abs_previous_value_line[x] > abs_current_value) 
                        first_scale = previous_scale_line[x];

                    if (first_scale <= max_scale) {
                        // Evaluate at first scale
                        int scale = first_scale; // alias
                        int yoffset = scale * widthStep;
                        int value;
#define CAM_RECURSIVE_PATTERN \
                        { \
                            unsigned int *ptrAi = ptr - (yoffset + scale); \
                            unsigned int *ptrDi = ptr + (yoffset + scale); \
                            unsigned int *ptrBi = ptr - (yoffset - scale); \
                            unsigned int *ptrCi = ptr + (yoffset - scale); \
                            unsigned int valin = *ptrDi - *ptrBi - *ptrCi + *ptrAi; \
                            unsigned int *ptrAo = ptr - ((yoffset + scale) << 1); \
                            unsigned int *ptrDo = ptr + ((yoffset + scale) << 1); \
                            unsigned int *ptrBo = ptr - ((yoffset - scale) << 1); \
                            unsigned int *ptrCo = ptr + ((yoffset - scale) << 1); \
                            unsigned int valout = *ptrDo - *ptrBo - *ptrCo + *ptrAo; \
                            value = valout - (valin << 2); \
			    if (abs(value) > 256 * height * width) printf("value : %i width : %i height : %i\n xOffset : %i yOffset : %i \n", value, width, height, iROI.srcroi.xOffset, iROI.srcroi.yOffset); \
                            assert(abs(value) < 256 * height * width);	\
                        }
                        CAM_RECURSIVE_PATTERN;
                        current_scale = scale;
                        current_value = value;
                        abs_current_value = abs(value);

                        // Let's see if we go upscale or not...
                        if (scale_up && scale < max_scale) {
                            // OK. Let's test the upscale
                            scale++;
                            yoffset += widthStep;
                            CAM_RECURSIVE_PATTERN;
                            s1 = current_scale * current_scale; s2 = scale * scale;
                            if (((CAM_INT64)abs(value)) * s1 > (CAM_INT64)abs_current_value * s2) {
                                // Yes. That was a good hint
                                current_scale = scale;
                                current_value = value;
                                abs_current_value = abs(value);
                                // Let's test again at higher scale
                                if (scale < max_scale) { 
                                    scale++;
                                    yoffset += widthStep;
                                    CAM_RECURSIVE_PATTERN;
                                    s1 = current_scale * current_scale; s2 = scale * scale;
                                    if (((CAM_INT64)abs(value)) * s1 > (CAM_INT64)abs_current_value * s2) {
                                        // Even better
                                        current_scale = scale;
                                        current_value = value;
                                        abs_current_value = abs(value);
                                        // And we stop there, satisfied...
                                    }
                                }
                            } else {
                                // Oh, oh... Apparently, that was not a good hint to go upscale...
                                // Let's try to downscale then...
                                scale -= 2;
                                if (scale > 0 && scale <= max_scale) {
                                    yoffset -= widthStep << 1;
                                    CAM_RECURSIVE_PATTERN;
                                    s1 = current_scale * current_scale; s2 = scale * scale;
                                    if (((CAM_INT64)abs(value)) * s1 > (CAM_INT64)abs_current_value * s2) {
                                        // Yes. It's better
                                        scale_up = 0; // This is the hint for next time
                                        current_scale = scale;
                                        current_value = value;
                                        abs_current_value = abs(value);
                                    }
                                    // And we stop there, satisfied...
                                }
                            }
                        } else {
                            // OK. Let's test the downscale
                            scale--;
                            if (scale > 0 && scale <= max_scale) {
                                yoffset -= widthStep;
                                CAM_RECURSIVE_PATTERN;
                                if (((CAM_INT64)abs(value)) * current_scale * current_scale > (CAM_INT64)abs_current_value * scale * scale) {
                                    // Yes. That was a good hint.
                                    current_scale = scale;
                                    current_value = value;
                                    abs_current_value = abs(value);
                                    // Let's test again at lower scale
                                    scale--;
                                    if (scale > 0 && scale <= max_scale) { 
                                        yoffset -= widthStep;
                                        CAM_RECURSIVE_PATTERN;
                                        s1 = current_scale * current_scale; s2 = scale * scale;
                                        if (((CAM_INT64)abs(value)) * s1 > (CAM_INT64)abs_current_value * s2) {
                                            // Even better
                                            current_scale = scale;
                                            current_value = value;
                                            abs_current_value = abs(value);
                                            // And we stop there, satisfied...
                                        }
                                    }
                                } else {
                                    // Oh, oh... Apparently, that was not a good hint to go downscale...
                                    // Let's try to upscale then...
                                    scale += 2;
                                    if (scale > 0 && scale <= max_scale) {
                                        yoffset += widthStep << 1;
                                        CAM_RECURSIVE_PATTERN;
                                        s1 = current_scale * current_scale; s2 = scale * scale;
                                        if (((CAM_INT64)abs(value)) * s1 > (CAM_INT64)abs_current_value * s2) {
                                            // Yes. It's better
                                            scale_up = 1; // This is the hint for next time
                                            current_scale = scale;
                                            current_value = value;
                                            abs_current_value = abs(value);
                                        }
                                        // And we stop there, satisfied...
                                    }
                                }
                            } else {
                                // Oh, oh... Apparently, that was not a good hint to go downscale...
                                // Let's try to upscale then...
                                scale += 2;
                                if (scale > 0 && scale <= max_scale) {
                                    yoffset += widthStep;
                                    CAM_RECURSIVE_PATTERN;
                                    s1 = current_scale * current_scale; s2 = scale * scale;
                                    if (((CAM_INT64)abs(value)) * s1 > (CAM_INT64)abs_current_value * s2) {
                                        // Yes. It's better
                                        scale_up = 1; // This is the hint for next time
                                        current_scale = scale;
                                        current_value = value;
                                        abs_current_value = abs(value);
                                    }
                                    // And we stop there, satisfied...
                                }
                            }
                        }
                    } else {
                        // Skip this pixel and go back to small scale
                        current_value = 0;
                        abs_current_value = 0;
                        current_scale = 1;
                        scale_up = 1;
                    }

                    current_value = (current_value << 4) / (current_scale * current_scale);
                    abs_current_value = abs(current_value);

                    // Now, we do have the current_value, current_scale and abs_current_value
                    // Let's check whether it is a local maximum or not
                    {
                        // If the current value is strictly greater than all the other values, then this is a local maximum and
                        // the other pixels are marked as NOT being local maxima
                        // If the current value is greater or equal to all neighbours, then this ia marked to be a local maximum
                        // but the other pixels are kept as local maxima
                        // If the current value is greater than any neighbour, this neighbour is marked as not being a local maximum
                        int local_maximum = 1;
                        if (x > 0) {
                            // Compare to the left and upper-left pixel
                            if (abs_current_value_line[x - 1] > abs_current_value)
                                local_maximum = 0;
                            else current_lmax_line[x - 1] = 0;
                            if (abs_previous_value_line[x - 1] > abs_current_value)
                                local_maximum = 0;
                            else previous_lmax_line[x - 1] = 0;
                        }   
                        // Compare to pixel above
                        if (abs_previous_value_line[x] > abs_current_value)
                            local_maximum = 0;
                        else previous_lmax_line[x] = 0;
                        if (x < width - 1) {
                            // Compare to the upper right pixel
                            if (abs_previous_value_line[x + 1] > abs_current_value)
                                local_maximum = 0;
                            else previous_lmax_line[x + 1] = 0;
                        }
                        if (local_maximum && abs_current_value > 0) {
                            current_lmax[current_nblmax++] = x;
                            current_lmax_line[x] = 1;
                        }
                    }

                    // Record the data for next line evaluation
                    current_value_line[x] = current_value;
                    abs_current_value_line[x] = abs_current_value;
                    current_scale_line[x] = current_scale;

                    //if (y == 100) {
                    //printf("y=%d x=%d value=%d abs=%d scale=%d\n", y, x, current_value, abs_current_value, current_scale);
                    //}
                }
            }

            // Check the local maxima on the previous line, and record the keypoints IF they are local maxima
            for (i = 0; i != previous_nblmax; i++) {
                if (previous_lmax_line[previous_lmax[i]]) {
                    int scale = previous_scale_line[previous_lmax[i]];
                    if (scale >= CAM_MIN_SCALE) {
                        // Yes, we have definitely a local maximum
                        // Record the keypoint
                        keypoints[nb_keypoints].x = previous_lmax[i];
                        keypoints[nb_keypoints].y = y - 1;
                        keypoints[nb_keypoints].scale = scale;
                        keypoints[nb_keypoints].value = previous_value_line[previous_lmax[i]];
                        nb_keypoints++;
                    }
                }
            }

            // Switch previous and current lines
            {
                unsigned char *tmp;
                tmp = current_scale_line;
                current_scale_line = previous_scale_line;
                previous_scale_line = tmp;
            }{
                unsigned int *tmp;
                tmp = abs_current_value_line;
                abs_current_value_line = abs_previous_value_line;
                abs_previous_value_line = tmp;
            }{
                int *tmp;
                tmp = current_value_line;
                current_value_line = previous_value_line;
                previous_value_line = tmp;
            }{
                unsigned char *tmp;
                tmp = current_lmax_line;
                current_lmax_line = previous_lmax_line;
                previous_lmax_line = tmp;
            }{
                unsigned int *tmp;
                tmp = current_lmax;
                current_lmax = previous_lmax;
                previous_lmax = tmp;
            }
            previous_nblmax = current_nblmax;
            current_nblmax = 0;
        }    
    }    

    for (i = 0; i < 2; i++) {
        free(value_lines[i]);
        free(abs_value_lines[i]);
        free(scale_lines[i]);
        free(lmax_lines[i]);
        free(lmax[i]);
    }

    // Post-processing : remove keypoints in a given neighborhood (proportionnal to scale)
    for (i = 1; i < nb_keypoints; i++) {
        CamKeypointShort *keypoint = &keypoints[i];
        int j = i - 1;
        int neighborhood = keypoint->scale >> 2;
        while (j >= 0 && keypoints[j].y >= keypoint->y - neighborhood) {
            CamKeypointShort *keypoint2 = &keypoints[j];
            if (keypoint2->x >= keypoint->x - neighborhood && keypoint2->x <= keypoint->x + neighborhood) {
               // They are next to each other...
               if (abs(keypoint2->value) > abs(keypoint->value))
                  keypoint->scale = 0; // Let's get rid of this one... It has a lower value 
               else
                  keypoint2->scale = 0;
            }
            j--;
        }
    }
    for (i = 0; i < nb_keypoints; i++) 
        if (!keypoints[i].scale) keypoints[i].value = 0;

    if (options & CAM_UPRIGHT) {
        // Remove keypoints the descriptor of which would be outside the frame boundaries
        for (i = 0; i < nb_keypoints; i++) {
            if (keypoints[i].scale) {
                int scale = keypoints[i].scale; // Save the scale
                keypoints[i].scale = (scale << 2) + 1;
                if (!camKeypointDescriptorCheckBounds(&keypoints[i], &integral)) 
                    keypoints[i].value = 0;
                keypoints[i].scale = scale;
            }
        }
    }

    // Sort the features according to value
    qsort(keypoints, nb_keypoints, sizeof(CamKeypointShort), camSortKeypointsShort);
    for (i = 0; i < nb_keypoints; i++) 
        if (keypoints[i].value == 0) break; 
    nb_keypoints = i;
    
    if (nb_keypoints > nb_max_keypoints) nb_keypoints = nb_max_keypoints;

    /* Interpolation :
     * Maxima : solve([y1=a*(x1-p)^2+b,y2=a*(-p)^2+b,y3=a*(x3-p)^2+b],[a,b,p])
     * Maxima yields the following formula for parabolic interpolation
                                       2         2     2         2
			             x1  y3 + (x3  - x1 ) y2 - x3  y1
			       p = ------------------------------------
                                   2 x1 y3 + (2 x3 - 2 x1) y2 - 2 x3 y1
    
     */
    // Scale super-resolution
    for (i = 0; i < nb_keypoints; i++) {
        CamKeypointShort *keypoint = &keypoints[i];
#if 1
        int x = keypoint->x, y = keypoint->y;
        int v, max_scale;    
        
        max_scale = (y + iROI.srcroi.yOffset) >> 1;
        v = (source->height - 1 - (y + iROI.srcroi.yOffset)) >> 1;
        if (v < max_scale) max_scale = v;
        v = (x + iROI.srcroi.xOffset) >> 1;
        if (v < max_scale) max_scale = v;
        v = (source->width -1 - (x + iROI.srcroi.xOffset)) >> 1;
        if (v < max_scale) max_scale = v;
        max_scale -= CAM_SCALE_MARGIN;

        if (keypoint->scale < max_scale) {
            unsigned int *ptr = ((unsigned int*)(integral.imageData + (y + iROI.srcroi.yOffset) * integral.widthStep)) + iROI.srcroi.xOffset + x;
            int scale = keypoint->scale - 1; 
            int yoffset = scale * widthStep;
            int value, num, den, x1, x3, y1, y2, y3;
            x1 = -1;
            x3 = 1;
            y2 = keypoint->value;
            CAM_RECURSIVE_PATTERN;
            y1 = (value << 4) / (scale * scale);
            if (y1 > y2) keypoint->scale <<= 2; // This is a problem : this is not a local maximum in scale...
            else {
                int found = 0;
                scale = keypoint->scale + 1;
                do {
                    yoffset = scale * widthStep;
                    CAM_RECURSIVE_PATTERN;
                    y3 = (value << 4) / (scale * scale);
                    if (y3 <= y2) {
                        num = (x1 * x1 * y3 + (x3 * x3 - x1 * x1) * y2 - x3 * x3 * y1);
                        den = x1 * y3 + (x3 - x1) * y2 - x3 * y1;
                        if (den == 0)
                            keypoint->scale <<= 2;
                        else {
                            int p = (num << 1) / den; // shift only by 1, because den is shifted by 2	    
                            keypoint->scale = p + ((scale - 1) << 2);
                            found = 1;
                        }
                        break;
                    }
                    y1 = y2;
                    y2 = y3;
                    scale++;
                } while (scale <= max_scale);
                if (!found) keypoint->scale <<= 2;
            }
        }
        else 
#endif
            keypoint->scale <<= 2;
    }

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
	points->keypoint[i]->x = keypoints[i].x + iROI.srcroi.xOffset;
	points->keypoint[i]->y = keypoints[i].y + iROI.srcroi.yOffset;
	points->keypoint[i]->scale = keypoints[i].scale;
	points->keypoint[i]->value = keypoints[i].value;
	points->keypoint[i]->angle = keypoints[i].angle;
    }
    points->nbPoints = pnb_keypoints;
    free(keypoints);

    if (options & CAM_UPRIGHT) {
        camKeypointsDescriptor(points, &integral, options);
    } else { 
        camKeypointsDescriptor(points, source, options);
    }

    // Finally, set the points' set 
    for (i = 0; i < points->nbPoints; i++) {
	points->keypoint[i]->set = points;
    }
    
    // Integral Image is now useless
    camDeallocateImage(&integral);
    camInternalROIPolicyExit(&iROI);
    return 1;
}

void test_camRecursiveKeypoints()
{
    CamImage image, dest;
    CamKeypoints points;
    int i;
    const int x = 8;
    int angle;
    double costheta;
    double sintheta;

    const int xp[4] = {-1, 1, 1, -1};
    const int yp[4] = {-1, -1, 1, 1};
    CamWarpingParams params;

    printf("Recursive keypoints detection :\n");
    camAllocateImage(&image, 256, 256, CAM_DEPTH_8U);
    camAllocateKeypoints(&points, 100);

    camSet(&image, 0);
    camDrawRectangle(&image, 102, 120, 156, 152, 128);
    camFillColor(&image, 103, 121, 128, -1);
    camDrawRectangle(&image, 100, 35, 150, 67, 128);
    camFillColor(&image, 123, 36, 128, -1);

    camDrawCircle(&image, 50, 50, 10, 128);
    camFillColor(&image, 50, 50, 128, -1);
    camDrawCircle(&image, 80, 50, 5, 128);
    camFillColor(&image, 80, 50, 128, -1);
    camDrawCircle(&image, 100, 50, 3, 128);
    camFillColor(&image, 100, 50, 128, -1);
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
	camDrawLine(&image, params.p[i].x, params.p[i].y, params.p[(i+1)%4].x, params.p[(i+1)%4].y, 128);
    }
    camFillColor(&image, 192, 192, 128, -1);
 
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
	camDrawLine(&image, params.p[i].x, params.p[i].y, params.p[(i+1)%4].x, params.p[(i+1)%4].y, 128);
    }
    camFillColor(&image, 50, 192, 128, -1);

#endif
    dest.imageData = NULL; 
    
    camKeypointsRecursiveDetector(&image, &points, 20, 0);
    for (i = 0; i < points.nbPoints; i++) {
	printf("x=%d y=%d value=%d scale=%d size=%d angle=%d\n", points.keypoint[i]->x, points.keypoint[i]->y, points.keypoint[i]->value, points.keypoint[i]->scale, points.keypoint[i]->size, points.keypoint[i]->angle);
    }
    camDrawKeypoints(&points, &image, 192);
    camSavePGM(&image, "output/keypoints_recursive.pgm");
    
    camDeallocateImage(&image);
    camDeallocateImage(&dest);
    camFreeKeypoints(&points);
}

