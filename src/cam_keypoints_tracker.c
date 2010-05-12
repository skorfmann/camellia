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

    Copyright (c) 2002-2006, Ecole des Mines de Paris - Centre de Robotique
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

        * Redistributions of source code must retain the above copyright
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

#include <stdlib.h>	// malloc, calloc, qsort
#include <string.h>	// memcpy
#include <stdio.h>	// printf
#ifdef __SSE2__
#include <emmintrin.h>	// _mm_malloc
#endif
#include "camellia.h"

extern double camSigmaParam;
#define CAM_ORIENTATION_STAMP_SIZE 30

typedef struct
{
  int width;	// half width reseach size
  int height;	// half height reseach size
  int scale;	// half scale reseach size
  int ds;	// scale amplification, positive to start search at bigger scale than previous one and reciprocally
} researchVolume;

typedef struct
{
  CamKeypoints		*previousFeatures;
  int			nbFrames;
  int			nbFeatures;
  int			nbSeeds;
  researchVolume	rv;		
} CamTrackingContext;

#define max(a, b) (a > b ? a : b)
#define min(a, b) (a > b ? b : a)

inline int	cam_keypoints_tracking_compute_detector(CamImage *integralImage, CamKeypointShort *keypoint)
{
  unsigned int	*ptr;
  unsigned int	*ptrAi;
  unsigned int	*ptrDi;
  unsigned int	*ptrBi;
  unsigned int	*ptrCi;
  unsigned int	valin;
  unsigned int	*ptrAo;
  unsigned int	*ptrDo;
  unsigned int	*ptrBo;
  unsigned int	*ptrCo;
  unsigned int	valout;
  int		res;
  int		yoffset;
  int		widthStep;
  int		scale;

  scale = keypoint->scale;
  ptr = ((unsigned int*)(integralImage->imageData + keypoint->y * integralImage->widthStep)) + keypoint->x;
  widthStep = integralImage->widthStep / 4;
  yoffset = scale * widthStep;
  ptrAi = ptr - (yoffset + scale);
  ptrDi = ptr + (yoffset + scale);
  ptrBi = ptr - (yoffset - scale);
  ptrCi = ptr + (yoffset - scale);
  valin = *ptrDi - *ptrBi - *ptrCi + *ptrAi;
  ptrAo = ptr - ((yoffset + scale) << 1);
  ptrDo = ptr + ((yoffset + scale) << 1);
  ptrBo = ptr - ((yoffset - scale) << 1);
  ptrCo = ptr + ((yoffset - scale) << 1);
  valout = *ptrDo - *ptrBo - *ptrCo + *ptrAo;
  res = valout - (valin << 2);
  return (res);
}

inline unsigned int	*cam_keypoints_tracking_extract_seeds(CamTrackingContext *tc)
{
  unsigned int		*seedsIndexes;
  unsigned int		currentIndex;
  int			i;
  int			j;

  seedsIndexes = (unsigned int*)malloc(tc->nbSeeds * sizeof(unsigned int));
  for (i = 0 ; i < tc->nbSeeds ; ++i)
    {
    _seedSelection:
      currentIndex = rand() % tc->nbFeatures;
      for (j = 0 ; j < i ; ++j)
	{
	  if (seedsIndexes[j] == currentIndex)
	    goto _seedSelection;
	}
      seedsIndexes[i] = currentIndex;
    }
  return (seedsIndexes);
}

int	camKeypointDescriptor(CamKeypoint *point, CamImage *source, CamImage *filter, int options);
int	camSortKeypointsShort(const void *p1x, const void *p2x);
int	camKeypointsDistance(CamKeypoint *point1, CamKeypoint *point2);
int	camKeypointOrientation(CamImage *source, CamKeypointShort *point, CamImage *filter, CamKeypointShort *next_point);

inline CamKeypointShort	*cam_keypoints_tracking_extract_max_on_each_scale(CamTrackingContext *tc, int (*detectorValue)(CamImage *, CamKeypointShort*), int i, CamImage *integralImage, unsigned int *seedsIndexes)
{
  int			j;
  int			k;
  int			seedX;
  int			seedY;
  int			seedScale;
  CamKeypointShort	currentPoint;
  CamKeypointShort	*maxOnEachScale;
  
  seedX = tc->previousFeatures->keypoint[seedsIndexes[i]]->x;
  seedY = tc->previousFeatures->keypoint[seedsIndexes[i]]->y;
  seedScale = tc->previousFeatures->keypoint[seedsIndexes[i]]->scale >> 2;
  maxOnEachScale = (CamKeypointShort*)calloc(((tc->rv.scale << 1) + 1), sizeof(CamKeypointShort));
  // on each scale, get the highest detector's value
  for (currentPoint.scale = max(seedScale - tc->rv.scale, 1), j = 0 ; currentPoint.scale <= seedScale + tc->rv.scale ; ++currentPoint.scale, ++j)
    {
      for (currentPoint.y = max(seedY - tc->rv.height, 0) ; currentPoint.y <= min(seedY + tc->rv.height, integralImage->height - 1) ; ++currentPoint.y)
	{
	  if ((currentPoint.y - (currentPoint.scale << 1)) <= 0)
	    break;
	  if ((currentPoint.y + (currentPoint.scale << 1)) >= integralImage->height)
	    break;
	  for (currentPoint.x = max(seedX - tc->rv.width, 0) ; currentPoint.x <= min(seedX + tc->rv.width, integralImage->width - 1) ; ++currentPoint.x)
	    {
	      if ((currentPoint.x - (currentPoint.scale << 1)) <= 0)
		break;
	      if ((currentPoint.x + (currentPoint.scale << 1)) >= integralImage->width)
		break;
	      currentPoint.value = abs(((*detectorValue)(integralImage, &currentPoint)));
	      if (currentPoint.value > maxOnEachScale[j].value)
		memcpy(&maxOnEachScale[j], &currentPoint, sizeof(CamKeypointShort));
	    }
	}
    }
  for (k = 0 ; k < j && maxOnEachScale[k].scale ; ++k)
	maxOnEachScale[k].value = (maxOnEachScale[k].value << 4) / (maxOnEachScale[k].scale * maxOnEachScale[k].scale);
  /* sur un bord ? => TO BE DONE */
  qsort(maxOnEachScale, j, sizeof(CamKeypointShort), camSortKeypointsShort);
  return (maxOnEachScale);
}

CamKeypointsMatches	*cam_keypoints_tracking_extract_seed_matches(CamTrackingContext *tc, CamImage *image, CamImage *integralImage, int (*detectorValue)(CamImage *, CamKeypointShort*),   CamKeypoints *currentFeatures, int options)
{
  int			i;
  int			l;
  int			distance1;
  int			distance2;
  CamKeypoint		possibleSeedMatch1;
  CamKeypoint		possibleSeedMatch2;
  CamROI		roi;
  CamImage		filter;
  CamKeypointShort	*maxOnEachScale;
  unsigned int		*seedsIndexes;
  CamKeypointsMatches	*seedsMatches;

  seedsIndexes = cam_keypoints_tracking_extract_seeds(tc);  

  seedsMatches = (CamKeypointsMatches*)malloc(sizeof(*seedsMatches));
  camAllocateKeypointsMatches(seedsMatches, tc->nbSeeds);
  
  // in the specified research volume, finds the more coherent point corresponding to the current seed <=> highest detector value
  l = 0;
  for (i = 0 ; i < tc->nbSeeds ; ++i)
    {
      maxOnEachScale = cam_keypoints_tracking_extract_max_on_each_scale(tc, detectorValue, i, integralImage, seedsIndexes);
      
      /* begin of angle computation of the 2 best */
      camAllocateImage(&filter, CAM_ORIENTATION_STAMP_SIZE, CAM_ORIENTATION_STAMP_SIZE, CAM_DEPTH_16S);
      camBuildGaussianFilter(&filter, camSigmaParam);
      camKeypointOrientation(image, &maxOnEachScale[0], &filter, &maxOnEachScale[1]); // check with bruno
      camKeypointOrientation(image, &maxOnEachScale[1], &filter, &maxOnEachScale[2]);
      camDeallocateImage(&filter);
      /* end of angle computation of the 2 best */
 
      /* begin check descriptor coherency */
      memcpy(&possibleSeedMatch1.x, &maxOnEachScale[0], sizeof(CamKeypointShort));
      memcpy(&possibleSeedMatch2.x, &maxOnEachScale[1], sizeof(CamKeypointShort));
      free(maxOnEachScale);
      possibleSeedMatch1.scale = possibleSeedMatch1.scale << 2;
      possibleSeedMatch2.scale = possibleSeedMatch2.scale << 2;
      camAllocateImage(&filter, 20, 20, CAM_DEPTH_16S);
      camBuildGaussianFilter(&filter, camSigmaParam);
      camKeypointsInternalsPrepareDescriptor();
      camKeypointDescriptor(&possibleSeedMatch1, image, &filter, options);
      camKeypointDescriptor(&possibleSeedMatch2, image, &filter, options);
      distance1 = camKeypointsDistance(tc->previousFeatures->keypoint[seedsIndexes[i]], &possibleSeedMatch1);
      distance2 = camKeypointsDistance(tc->previousFeatures->keypoint[seedsIndexes[i]], &possibleSeedMatch2);
      camDeallocateImage(&filter);
      
      if (distance1 < distance2)
	{
	  currentFeatures->keypoint[l] = &currentFeatures->bag[l];
	  memcpy(currentFeatures->keypoint[l], &possibleSeedMatch1, sizeof(possibleSeedMatch1));
	  seedsMatches->pairs[l].p1 = tc->previousFeatures->keypoint[seedsIndexes[i]];
	  seedsMatches->pairs[l].p2 = currentFeatures->keypoint[l];
	  seedsMatches->pairs[l].mark = distance1;
	  seedsMatches->nbMatches++;
	  ++l;
	}
      else
	seedsMatches->nbOutliers++;

      /* end check descriptor coherency */
    }
  free (seedsIndexes);
  return (seedsMatches);
}

CamKeypointsMatches	*cam_keypoints_tracking(CamTrackingContext *tc, CamImage *image, int (*roiDetector)(CamImage*, CamKeypoints*, int, int), int (*detectorValue)(CamImage *, CamKeypointShort*), int options)
{
  CamImage		integralImage;
  CamKeypointsMatches	*seedsMatches;
  CamKeypoints		*currentFeatures;

  // initialisation of the tracker
  if (!(tc->previousFeatures->bag))
    {
      //srandom(time(NULL)); // used to extract seeds
      (*roiDetector)(image, tc->previousFeatures, tc->nbFeatures, options);
      return (NULL);
    }
  // initialisation already done
  else
    {
      currentFeatures = (CamKeypoints*)malloc(sizeof(CamKeypoints));
      camAllocateKeypoints(currentFeatures, tc->nbFeatures);
      if (currentFeatures->bag == NULL)
	{
#ifdef __SSE2__
	  currentFeatures->bag = (CamKeypoint*)_mm_malloc(sizeof(CamKeypoint) * tc->nbFeatures, 16);
#else
	  currentFeatures->bag = (CamKeypoint*)malloc(sizeof(CamKeypoint) * tc->nbFeatures);
#endif
      }
      
      /* begin integral image computing */
      integralImage.imageData = NULL;
      camIntegralImage(image, &integralImage);
      /* end integral image computing */

      seedsMatches = cam_keypoints_tracking_extract_seed_matches(tc, image, &integralImage, detectorValue, currentFeatures, options); // need to check the descriptor computation
      

      camFreeKeypoints(tc->previousFeatures);
      free(tc->previousFeatures);
      tc->previousFeatures = currentFeatures;
      free (seedsMatches);
      camDeallocateImage(&integralImage);
      return (NULL);
    }
  return (NULL);
}

int camKeypointsDetector(CamImage *source, CamKeypoints *points, int nb_max_keypoints, int options);
int camKeypointsRecursiveDetector(CamImage *source, CamKeypoints *points, int nb_max_keypoints, int options);

void			test_cam_keypoints_tracking()
{
  CamTrackingContext	tc;
  CamImage		modelImage;
  CamImage		firstImage;
  CamImage		secondImage;
  int			t1;
  int			t2;
  int			t3;

  /* begin images building */
  modelImage.imageData = NULL;
  camLoadBMP(&modelImage, "./resources/photos/scene1.bmp");
  camAllocateYUVImage(&firstImage, modelImage.width, modelImage.height);
  camRGB2YUV(&modelImage, &firstImage);
  camLoadBMP(&modelImage, "./resources/photos/scene1.bmp");
  camAllocateYUVImage(&secondImage, modelImage.width, modelImage.height);
  camRGB2YUV(&modelImage, &secondImage);
  /* end images building */

  /* begin tracking context */
  tc.nbFeatures = 100;
  tc.nbFrames = 3;
  tc.nbSeeds = 10;
  tc.rv.width = 1;
  tc.rv.height = 1;
  tc.rv.scale = 1;
  tc.rv.ds = 0;
  tc.previousFeatures = (CamKeypoints*)malloc(sizeof(CamKeypoints));
  camAllocateKeypoints(tc.previousFeatures, 100);
  /* end tracking context */

  t1 = camGetTimeMs();
  cam_keypoints_tracking(&tc, &firstImage, camKeypointsRecursiveDetector, cam_keypoints_tracking_compute_detector, 0);
  t2 = camGetTimeMs();
  cam_keypoints_tracking(&tc, &secondImage, camKeypointsRecursiveDetector, cam_keypoints_tracking_compute_detector, 0);
  t3 = camGetTimeMs();
  printf("initial detection timing : %ims\ntracking timing : %ims\n", t2 - t1, t3 - t2);
}
