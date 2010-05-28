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
#include <sys/time.h>	// gettimeofday
#include <limits.h>	// INT_MAX
#include <math.h>	// sqrt
#ifdef __SSE2__
#include <emmintrin.h>	// _mm_malloc
#endif
#include "camellia.h"

extern double camSigmaParam;

#define CAM_ORIENTATION_STAMP_SIZE 30
//#define CAM_TRACKING_SUBTIMINGS
#define CAM_TRACKING_TIMINGS
#define CAM_TRACKING_DEBUG_1

typedef	struct
{
  int	width;	// half width reseach size
  int	height;	// half height reseach size
  int	scale;	// half scale reseach size
  int	ds;	// scale amplification, positive to start search at bigger scale than previous one and reciprocally, not yet used
} researchVolume;

typedef	struct
{
  int	height;
  int	width;
} researchWindow;

typedef struct
{
  CamKeypoints		*previousFeatures;
  int			nbFrames;		//not yet used
  int			nbFeatures;
  int			nbSeeds;
  researchVolume	rv;			// volume around feature point
  researchWindow	rw;			// window for the seeding
} CamTrackingContext;

typedef struct		s_pointList
{
  int			x;
  int			y;
  int			scale;
  int			value;
  struct s_pointList	*next;
} CamPointsList;

typedef enum
  {
    TRUE,
    FALSE
  } BOOL;

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

void		cam_keypoints_tracking_print_matches(CamImage *img1, CamImage *img2, char *outfile, CamKeypointsMatches *matches)
{
  CamImage	res;
  CamROI	roi;
  char		filename[256];
  register int	i;
  int		x1;
  int		x2;
  int		y1;
  int		y2;

  camAllocateRGBImage(&res, img1->width, img1->height * 2);
  camSetROI(&roi, 0, 0, img1->height, img1->width, img1->height);
  res.roi = NULL;
  camCopy(img1, &res);
  res.roi = &roi;
  camCopy(img2, &res);
  res.roi = NULL;

  for (i = 0 ; i < matches->nbMatches ; ++i)
    {
      camDrawKeypoint(matches->pairs[i].p1, &res, CAM_RGB(255, 0, 0));
      x1 = matches->pairs[i].p1->x;
      y1 = matches->pairs[i].p1->y;
      x2 = matches->pairs[i].p2->x;
      y2 = matches->pairs[i].p2->y;
      y2 += img1->height;
      camDrawLine(&res, x1, y1, x2, y2, CAM_RGB(0, 255, 0));
    }
  for (i = 0 ; i < matches->nbMatches ; ++i)
    {
      matches->pairs[i].p2->y += img1->height;
      camDrawKeypoint(matches->pairs[i].p2, &res, 128);
    }

  sprintf(filename, "output/%s.bmp", outfile);
  camSaveBMP(&res, filename);
  camDeallocateImage(&res);
}

inline unsigned int	*cam_keypoints_tracking_extract_seeds(CamTrackingContext *tc)
{
  unsigned int		*seedsIndexes;
  register unsigned int	i;
  register unsigned int	j;
  register unsigned int	k;
  unsigned int			step;
  
  step = sqrt(tc->nbFeatures) / sqrt(tc->nbSeeds);
  seedsIndexes = (unsigned int*)malloc(tc->nbSeeds * sizeof(unsigned int));
  k = 0;
  for (i = step - 1 ; i < sqrt(tc->nbFeatures) ; i += step)
    {
      for (j = step - 1 ; j < sqrt(tc->nbFeatures) && k < tc->nbSeeds ; j += step)
	seedsIndexes[k++] = i * sqrt(tc->nbFeatures) + j;
    }

  return (seedsIndexes);
}

int	camKeypointDescriptor(CamKeypoint *point, CamImage *source, CamImage *filter, int options);
int	camSortKeypointsShort(const void *p1x, const void *p2x);
int	camKeypointsDistance(CamKeypoint *point1, CamKeypoint *point2);
int	camKeypointOrientation(CamImage *source, CamKeypointShort *point, CamImage *filter, CamKeypointShort *next_point);

inline CamKeypointShort	*cam_keypoints_tracking_extract_max_on_each_scale(CamTrackingContext *tc, int (*detectorValue)(CamImage *, CamKeypointShort*), CamImage *integralImage, int seedX, int seedY, int seedScale)
{
  int			j;
  int			k;
  CamKeypointShort	currentPoint;
  CamKeypointShort	*maxOnEachScale;
  
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
  qsort(maxOnEachScale, j, sizeof(CamKeypointShort), camSortKeypointsShort);
  return (maxOnEachScale);
}

inline int	cam_keypoints_tracking_max_on_border(CamTrackingContext *tc, CamKeypointShort *localMax, int seedX, int seedY, int seedScale)
{
  if ( abs(localMax->x - seedX) == tc->rv.width || abs(localMax->y - seedY) == tc->rv.height ||
       abs(localMax->scale - seedScale) == tc->rv.scale )
    return (1);
  return (0);
}

inline void	cam_keypoints_tracking_extract_new_research_box(CamTrackingContext *tc, CamKeypointShort *localMax, int *seedX, int *seedY, int *seedScale)
{
  if (*seedX - localMax->x == -(tc->rv.width))
    *seedX += tc->rv.width << 1;
  else if (*seedX - localMax->x == (tc->rv.width))
    *seedX -= tc->rv.width << 1;
  if (*seedY - localMax->y == -(tc->rv.height))
    *seedY += tc->rv.height << 1;
  else if (*seedY - localMax->y == (tc->rv.height))
    *seedY -= tc->rv.height << 1;
  if (*seedScale - localMax->scale == -(tc->rv.scale))
    *seedScale += tc->rv.scale << 1;
  else if (*seedScale - localMax->scale == (tc->rv.scale))
    *seedScale -= tc->rv.scale << 1;
}

void		cam_keypoints_tracking_free_linked_list(CamPointsList *l)
{
  CamPointsList	*ptr;

  ptr = l;
  while (ptr)
    {
      l = l->next;
      free (ptr);
      ptr = l;
    }
}

inline CamPointsList	*cam_keypoints_tracking_add_to_linked_list(CamPointsList *l, int x, int y, int scale, int value)
{
  CamPointsList		*head;

  head = (CamPointsList*)malloc(sizeof(CamPointsList));
  head->x = x;
  head->y = y;
  head->scale = scale;
  head->value = value;
  head->next = l;
  return (head);
}

inline BOOL	cam_keypoints_tracking_point_not_visited(CamPointsList *l)
{
  CamPointsList	*point;

  point = l->next;
  while (point)
    {
      if (l->x == point->x && l->y == point->y && l->scale == point->scale)
	return (FALSE);
      point = point->next;
    }
  return (TRUE);
}


inline CamKeypointShort	*cam_keypoints_tracking_extract_overall_max_on_each_scale(CamTrackingContext *tc, CamImage *integralImage, int (*detectorValue)(CamImage *, CamKeypointShort*), int index, int shiftX, int shiftY, int shiftScale)
{
  int			featureX;
  int			featureY;
  int			featureScale;
  CamKeypointShort	*maxOnEachScale;
  CamPointsList		*pointsList;
#ifdef CAM_TRACKING_SUBTIMINGS
  int			deltaTimers;
  struct timeval	tv1;
  struct timeval	tv2;
#endif
  
  featureX = tc->previousFeatures->keypoint[index]->x + shiftX;
  featureY = tc->previousFeatures->keypoint[index]->y + shiftY;
  featureScale = max((tc->previousFeatures->keypoint[index]->scale + shiftScale) >> 2, 1);
  
#ifdef CAM_TRACKING_SUBTIMINGS
  gettimeofday(&tv1, NULL);
#endif
  maxOnEachScale = cam_keypoints_tracking_extract_max_on_each_scale(tc, detectorValue, integralImage, featureX, featureY, featureScale);
  pointsList = NULL;
  pointsList = cam_keypoints_tracking_add_to_linked_list(pointsList, maxOnEachScale[0].x, maxOnEachScale[0].y, maxOnEachScale[0].scale, maxOnEachScale[0].value);
  /* max on a border of the research space => shift the research cube */
  while (cam_keypoints_tracking_max_on_border(tc, &maxOnEachScale[0], featureX, featureY, featureScale) && cam_keypoints_tracking_point_not_visited(pointsList))
    {
      cam_keypoints_tracking_extract_new_research_box(tc, &maxOnEachScale[0], &featureX, &featureY, &featureScale);
      free(maxOnEachScale);
      maxOnEachScale = cam_keypoints_tracking_extract_max_on_each_scale(tc, detectorValue, integralImage, featureX, featureY, featureScale);
      pointsList = cam_keypoints_tracking_add_to_linked_list(pointsList, maxOnEachScale[0].x, maxOnEachScale[0].y, maxOnEachScale[0].scale, maxOnEachScale[0].value);
    }
#ifdef CAM_TRACKING_SUBTIMINGS
  gettimeofday(&tv2, NULL);
  deltaTimers = tv2.tv_usec - tv1.tv_usec;
  if (deltaTimers < 0)
    deltaTimers = 1000000 + deltaTimers;
  printf("Max on each scale : %ius\n", deltaTimers);
#endif
  return (maxOnEachScale);
}

void			cam_keypoints_tracking_compute_feature_description(CamKeypoint *keypoint, CamKeypointShort *feature, CamKeypointShort *nextFeature, CamImage *image, CamImage *integralImage, int options)
{
  CamImage		filter;  

  if (!(options & CAM_UPRIGHT))
    {
      /* begin of angle computation of the best */
      camAllocateImage(&filter, CAM_ORIENTATION_STAMP_SIZE, CAM_ORIENTATION_STAMP_SIZE, CAM_DEPTH_16S);
      camBuildGaussianFilter(&filter, camSigmaParam);
      camKeypointOrientation(image, feature, &filter, nextFeature);
      camDeallocateImage(&filter);
      /* end of angle computation of the best */
    }
  camAllocateImage(&filter, 20, 20, CAM_DEPTH_16S);
  camBuildGaussianFilter(&filter, camSigmaParam);
  camKeypointsInternalsPrepareDescriptor();
  if (options & CAM_UPRIGHT)
    {
      camKeypointDescriptor(keypoint, integralImage, &filter, options);
      keypoint->angle = 0;
    }
  else
    camKeypointDescriptor(keypoint, image, &filter, options);
  camDeallocateImage(&filter);
}

CamKeypointsMatches	*cam_keypoints_tracking_extract_seed_matches(CamTrackingContext *tc, CamImage *image, CamImage *integralImage, int (*detectorValue)(CamImage *, CamKeypointShort*),   CamKeypoints *currentFeatures, unsigned int *seedsIndexes, int options)
{
  register int		i;
  register int		l;
  int			distance1;
  int			distance2;
  CamKeypoint		seedMatch;
  CamROI		roi;
  CamImage		filter;
  CamKeypointShort	*maxOnEachScale;
  CamKeypointsMatches	*seedsMatches;
#ifdef CAM_TRACKING_SUBTIMINGS
  int			deltaTimers;
  struct timeval	tv1;
  struct timeval	tv2;
#endif

  seedsMatches = (CamKeypointsMatches*)malloc(sizeof(*seedsMatches));
  camAllocateKeypointsMatches(seedsMatches, tc->nbSeeds);
  
  // in the specified research volume, finds the more coherent point corresponding to the current seed <=> highest detector value
  l = 0;
  for (i = 0 ; i < tc->nbSeeds ; ++i)
    {
      maxOnEachScale = cam_keypoints_tracking_extract_overall_max_on_each_scale(tc, integralImage, detectorValue, seedsIndexes[i], 0, 0, 0);
#ifdef CAM_TRACKING_SUBTIMINGS
      gettimeofday(&tv1, NULL);      
#endif
      memcpy(&seedMatch.x, &maxOnEachScale[0], sizeof(CamKeypointShort));
      cam_keypoints_tracking_compute_feature_description(&seedMatch, &maxOnEachScale[0], &maxOnEachScale[1], image, integralImage, options);
      seedMatch.scale = seedMatch.scale << 2;
#ifdef CAM_TRACKING_SUBTIMINGS
      gettimeofday(&tv2, NULL);      
#endif
      currentFeatures->keypoint[l] = &currentFeatures->bag[l];
      memcpy(currentFeatures->keypoint[l], &seedMatch, sizeof(seedMatch));
      ++currentFeatures->nbPoints;
      seedsMatches->pairs[l].p1 = tc->previousFeatures->keypoint[seedsIndexes[i]];
      seedsMatches->pairs[l].p2 = currentFeatures->keypoint[l];
      seedsMatches->pairs[l].error = camKeypointsDistance(tc->previousFeatures->keypoint[seedsIndexes[i]], currentFeatures->keypoint[l]);
      seedsMatches->pairs[l].mark = 1;
      seedsMatches->nbMatches++;
      ++l;
#ifdef CAM_TRACKING_SUBTIMINGS
      deltaTimers = tv2.tv_usec - tv1.tv_usec;
      if (deltaTimers < 0)
	deltaTimers = 1000000 + deltaTimers;
      printf("Feature description : %ius\n", deltaTimers);
#endif
      /* end computation of descriptor */
    }
  return (seedsMatches);
}

BOOL			cam_keypoints_tracking_is_seed(CamTrackingContext *tc, unsigned int *seedsIndexes, unsigned int featureIndex)
{
  register unsigned int	k;
  
  for (k = 0 ; k < tc->nbSeeds ; ++k)
    {
      if (seedsIndexes[k] == featureIndex)
	return (TRUE);
    }  
  return (FALSE);
}

inline unsigned int	cam_keypoints_tracking_distance_square(CamKeypoint *p1, CamKeypoint *p2)
{
  return ((p1->x - p2->x) * (p1->x - p2->x) + (p1->y - p2->y) * (p1->y - p2->y));
}

unsigned int		cam_keypoints_tracking_extract_closest_seed_index(CamTrackingContext *tc, CamKeypointsMatches *seedsMatches, unsigned int featureIndex)
{
  register unsigned int	i;
  unsigned int		res;
  float			minDistanceSquare;
  float			currentDistance;

  minDistanceSquare = INT_MAX;
  for (i = 0 ; i < seedsMatches->nbMatches ; ++i)
    {
      currentDistance = cam_keypoints_tracking_distance_square(seedsMatches->pairs[i].p1, tc->previousFeatures->keypoint[featureIndex]);
      if (currentDistance < minDistanceSquare)
	{
	  minDistanceSquare = currentDistance;
	  res = i;
	}
    }
  return (res);
}

CamKeypointsMatches	*cam_keypoints_tracking_extract_points_matching(CamTrackingContext *tc, CamImage *image, CamImage *integralImage, int (*detectorValue)(CamImage *, CamKeypointShort*), CamKeypointsMatches *seedsMatches, unsigned int *seedsIndexes, CamKeypoints *currentFeatures, int options)
{

  register unsigned int	i;
  register unsigned int	j;
  register unsigned int	k;
  unsigned int		closestSeedIndex;
  CamKeypointShort	*maxOnEachScale;
  CamKeypoint		keypoint;
  CamKeypointsMatches	*featuresMatches;
  CamKeypointsMatches	*res;
  int			shiftX;
  int			shiftY;
  int			shiftScale;
  int			supposedMinDistance;
  int			currentDistance;
  int			*goodMatchesTab;
  int			nbGoodMatches;

  featuresMatches = (CamKeypointsMatches*)malloc(sizeof(*featuresMatches));
  camAllocateKeypointsMatches(featuresMatches, tc->nbFeatures - tc->nbSeeds);
  goodMatchesTab = (int*)calloc(tc->nbFeatures - tc->nbSeeds, sizeof(int));

  j = 0;
  for (i = 0 ; i < tc->nbFeatures ; ++i)
    {
      if (cam_keypoints_tracking_is_seed(tc, seedsIndexes, i) == TRUE)
	continue ;
      closestSeedIndex = cam_keypoints_tracking_extract_closest_seed_index(tc, seedsMatches, i);
      shiftX = seedsMatches->pairs[closestSeedIndex].p2->x - seedsMatches->pairs[closestSeedIndex].p1->x;
      shiftY = seedsMatches->pairs[closestSeedIndex].p2->y - seedsMatches->pairs[closestSeedIndex].p1->y;
      shiftScale = seedsMatches->pairs[closestSeedIndex].p2->scale - seedsMatches->pairs[closestSeedIndex].p1->scale;
      maxOnEachScale = cam_keypoints_tracking_extract_overall_max_on_each_scale(tc, integralImage, detectorValue, i, shiftX, shiftY, shiftScale);
      memcpy(&keypoint.x, &maxOnEachScale[0], sizeof(CamKeypointShort));
      keypoint.scale = keypoint.scale << 2;
      cam_keypoints_tracking_compute_feature_description(&keypoint, &maxOnEachScale[0], &maxOnEachScale[1], image, integralImage, options);
      free(maxOnEachScale);
      featuresMatches->pairs[j].p1 = tc->previousFeatures->keypoint[i];
      currentFeatures->keypoint[currentFeatures->nbPoints] = &currentFeatures->bag[currentFeatures->nbPoints];
      memcpy(currentFeatures->keypoint[currentFeatures->nbPoints++], &keypoint, sizeof(keypoint));
      featuresMatches->pairs[j].p2 = currentFeatures->keypoint[j + seedsMatches->nbMatches];
      ++j;
    }
  /* begin of checking for each point if the found one is the good one */
  for (i = 0 ; i < j ; ++i)
    {
      supposedMinDistance = camKeypointsDistance(featuresMatches->pairs[i].p1, featuresMatches->pairs[i].p2);
      for (k = 0 ; k < j ; ++k)
	{
	  currentDistance = camKeypointsDistance(featuresMatches->pairs[i].p1, featuresMatches->pairs[k].p2);
	  if (supposedMinDistance > currentDistance)
	    break;
	}
      if (k == j)
	{
	  featuresMatches->pairs[i].mark = 1;
	  featuresMatches->pairs[i].error = supposedMinDistance;
	  ++featuresMatches->nbMatches;
	}
      else
	{
	  featuresMatches->pairs[i].mark = 0;
	  featuresMatches->pairs[i].error = currentDistance; // the first better found, maybe others
	  ++featuresMatches->nbOutliers;
	}
    }
  /* end of checking for each point if the found one is the good one */
  return (featuresMatches);
}

void	cam_keypoints_tracking_fill_empty_area()
{

}

CamKeypointsMatches	*cam_keypoints_tracking(CamTrackingContext *tc, CamImage *image, int (*roiDetector)(CamImage *,CamImage *, CamKeypoints *, int, int), int (*detectorValue)(CamImage *, CamKeypointShort*), int options)
{
  CamImage		integralImage;
  CamKeypointsMatches	*seedsMatches;
  CamKeypointsMatches	*keypointsMatches;
  CamKeypoints		*currentFeatures;
  unsigned int		*seedsIndexes;
  CamKeypointsMatches	*res;
  CamKeypoints		tmpPoint;
  CamROI		roi;
  CamROI		roix;
  register int		i;
  register int		j;
  
#ifdef CAM_TRACKING_TIMINGS
  int			timeInROIDetector;
  int			t1;
  int			t2;
  int			t3;
  int			t4;
  int			t5;
#endif

  // initialisation of the tracker
  if (!(tc->previousFeatures))
    {
      tc->previousFeatures = (CamKeypoints*)malloc(sizeof(CamKeypoints));
      camAllocateKeypoints(tc->previousFeatures, tc->nbFeatures);
#ifdef __SSE2__
      tc->previousFeatures->bag = (CamKeypoint*)_mm_malloc(sizeof(CamKeypoint) * tc->nbFeatures, 16);
#else
      tc->previousFeatures->bag = (CamKeypoint*)malloc(sizeof(CamKeypoint) * tc->nbFeatures);
#endif

      camAllocateKeypoints(&tmpPoint, 1);

      /* computation of integral image */
      integralImage.imageData = NULL;
      camSetMaxROI(&roix, image);
      roix.coi = 1;
      image->roi = &roix;
      camIntegralImage(image, &integralImage);
      integralImage.roi = &roi;

      // configuration of the window related to the seeding
      image->roi = &roi;
      image->roi->width = tc->rw.width;
      image->roi->height = tc->rw.height;
      i = 0;
      timeInROIDetector = 0;
      for (image->roi->xOffset = 0 ; image->roi->xOffset < image->width ; image->roi->xOffset += tc->rw.width)
	{
	  for (image->roi->yOffset = 0 ; image->roi->yOffset < image->height ; image->roi->yOffset += tc->rw.height)
	    {	
#ifdef CAM_TRACKING_SUBTIMINGS
	      t1 = camGetTimeMs();
#endif
	      (*roiDetector)(image, &integralImage, &tmpPoint, 1, options);
#ifdef CAM_TRACKING_SUBTIMINGS
	      t2 = camGetTimeMs();
	      timeInROIDetector += t2 - t1;
#endif
	      tc->previousFeatures->keypoint[i] = &tc->previousFeatures->bag[i];
	      memcpy(tc->previousFeatures->keypoint[i], tmpPoint.keypoint[0], sizeof(CamKeypoint));
	      ++i;
	    }
	}
#ifdef CAM_TRACKING_SUBTIMINGS
      printf("Time spent in roiDetector : %ims\n", timeInROIDetector);
#endif
      camFreeKeypoints(&tmpPoint);
      camDeallocateImage(&integralImage);
      return (NULL);
    }
  // initialisation already done
  else
    {
      seedsIndexes = cam_keypoints_tracking_extract_seeds(tc);  
      currentFeatures = (CamKeypoints*)malloc(sizeof(CamKeypoints));
      res = (CamKeypointsMatches*)malloc(sizeof(CamKeypointsMatches));
      camAllocateKeypoints(currentFeatures, tc->nbFeatures);
      if (currentFeatures->bag == NULL)
	{
#ifdef __SSE2__
	  currentFeatures->bag = (CamKeypoint*)_mm_malloc(sizeof(CamKeypoint) * tc->nbFeatures, 16);
#else
	  currentFeatures->bag = (CamKeypoint*)malloc(sizeof(CamKeypoint) * tc->nbFeatures);
#endif
	}

#ifdef CAM_TRACKING_TIMINGS
      t4 = camGetTimeMs();
#endif
      /* begin integral image computing */ 
      integralImage.imageData = NULL;
      camIntegralImage(image, &integralImage); 
#ifdef CAM_TRACKING_TIMINGS
      t5 = camGetTimeMs();
#endif
     /* end integral image computing */
#ifdef CAM_TRACKING_TIMINGS
      t1 = camGetTimeMs();
#endif
      seedsMatches = cam_keypoints_tracking_extract_seed_matches(tc, image, &integralImage, detectorValue, currentFeatures, seedsIndexes, options); // need to check the descriptor computation
#ifdef CAM_TRACKING_TIMINGS
      t2 = camGetTimeMs();
#endif
      keypointsMatches = cam_keypoints_tracking_extract_points_matching(tc, image, &integralImage, detectorValue, seedsMatches, seedsIndexes, currentFeatures, options);
#ifdef CAM_TRACKING_TIMINGS
      t3 = camGetTimeMs();
      printf("Seed matches : %ims, features matches : %ims, integral image : %ims\n", t2 - t1, t3 - t2, t5 - t4);
#endif

      camAllocateKeypointsMatches(res, seedsMatches->nbMatches + keypointsMatches->nbMatches);

      j = 0;
      for (i = 0 ; i < tc->nbSeeds ; ++i)
	{
	  if (seedsMatches->pairs[i].mark)
	    {
	      res->pairs[j].p1 = (CamKeypoint*)malloc(sizeof(CamKeypoint));
	      res->pairs[j].p2 = (CamKeypoint*)malloc(sizeof(CamKeypoint));
	      memcpy(res->pairs[j].p1, seedsMatches->pairs[i].p1, sizeof(CamKeypoint));
	      memcpy(res->pairs[j].p2, seedsMatches->pairs[i].p2, sizeof(CamKeypoint));
	      ++j;
	    }
	}
      for (i = 0 ; i < tc->nbFeatures - tc->nbSeeds ; ++i)
	{
	  if (keypointsMatches->pairs[i].mark)
	    {
	      res->pairs[j].p1 = (CamKeypoint*)malloc(sizeof(CamKeypoint));
	      res->pairs[j].p2 = (CamKeypoint*)malloc(sizeof(CamKeypoint));
	      memcpy(res->pairs[j].p1, keypointsMatches->pairs[i].p1, sizeof(CamKeypoint));
	      memcpy(res->pairs[j].p2, keypointsMatches->pairs[i].p2, sizeof(CamKeypoint));
	      ++j;
	    }
	}
      res->nbMatches = seedsMatches->nbMatches + keypointsMatches->nbMatches;
      res->nbOutliers = tc->nbFeatures - seedsMatches->nbMatches + keypointsMatches->nbMatches;

      // cam_keypoints_tracking_clean_keypoints(currentFeatures, keypointsMatches, seedsMatches->nbMatches); => to be done
      cam_keypoints_tracking_fill_empty_area();
      
      free (seedsIndexes);
      camFreeKeypoints(tc->previousFeatures);
      free(tc->previousFeatures);
      camFreeKeypointsMatches(seedsMatches);
      free (seedsMatches);

      tc->previousFeatures = currentFeatures;
      camFreeKeypointsMatches(keypointsMatches);
      free (keypointsMatches);
      camDeallocateImage(&integralImage);
      return (res);
    }
}

//int camKeypointsDetector(CamImage *source, CamKeypoints *points, int nb_max_keypoints, int options);
int camKeypointsRecursiveDetector(CamImage *source, CamImage *integral, CamKeypoints *points, int nb_max_keypoints, int options);

void			test_cam_keypoints_tracking()
{
  CamTrackingContext	tc;
  CamImage		modelImage;
  CamImage		firstImage;
  CamImage		secondImage;
  int			t1;
  int			t2;
  int			t3;
  CamKeypointsMatches	*track;
  char			img1[] = "./resources/rover/rotation1.bmp";
  char			img2[] = "./resources/rover/rotation2.bmp";

  /* begin images building */
  modelImage.imageData = NULL;
  //camLoadBMP(&modelImage, "./resources/photos/scene1.bmp");
  camLoadBMP(&modelImage, img1);
  camAllocateYUVImage(&firstImage, modelImage.width, modelImage.height);
  camRGB2YUV(&modelImage, &firstImage);
  //camLoadBMP(&modelImage, "./resources/photos/scene1.bmp");
  camLoadBMP(&modelImage, img2);
  camAllocateYUVImage(&secondImage, modelImage.width, modelImage.height);
  camRGB2YUV(&modelImage, &secondImage);
  /* end images building */

  /* begin tracking context */
  tc.nbFeatures = 100;
  tc.nbFrames = 3;
  tc.nbSeeds = sqrt(tc.nbFeatures) - 1;
  tc.rv.width = 2;
  tc.rv.height = 2;
  tc.rv.scale = 2;
  tc.rv.ds = 0;
  tc.rw.height = firstImage.height / sqrt(tc.nbFeatures);
  tc.rw.width = firstImage.width / sqrt(tc.nbFeatures);
  tc.previousFeatures = NULL;
  /* end tracking context */

#ifdef CAM_TRACKING_TIMINGS
  t1 = camGetTimeMs();
#endif
  
  cam_keypoints_tracking(&tc, &firstImage, camKeypointsRecursiveDetector, cam_keypoints_tracking_compute_detector, CAM_UPRIGHT);

#ifdef CAM_TRACKING_TIMINGS
  t2 = camGetTimeMs();
#endif

  track = cam_keypoints_tracking(&tc, &secondImage, camKeypointsRecursiveDetector, cam_keypoints_tracking_compute_detector, CAM_UPRIGHT);
  camDeallocateImage(&firstImage);
  camDeallocateImage(&secondImage);

#ifdef CAM_TRACKING_TIMINGS
  t3 = camGetTimeMs();
  printf("initial detection timing : %ims\ntracking timing : %ims\n", t2 - t1, t3 - t2);
#endif

  camLoadBMP(&firstImage, img1);
  camLoadBMP(&secondImage, img2);
  cam_keypoints_tracking_print_matches(&firstImage, &secondImage, "track", track);
  camDeallocateImage(&firstImage);
  camDeallocateImage(&secondImage);

  // camfullFREE(track)
  camFreeKeypointsMatches(track);
  free(track);
}
