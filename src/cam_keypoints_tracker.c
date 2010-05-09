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

#include <stdlib.h>
#include <string.h> // memcpy
#include "camellia.h"

typedef struct
{
  int width;	// half width reseach size
  int height;	// half height reseach size
  int scale;	// half scale reseach size
  int ds;	// scale amplification, positive to start search at bigger scale than previous one and reciprocally
} researchVolume;

typedef struct
{
  CamKeypoints		previousFeatures;
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

  ptr = (unsigned int*)(integralImage->imageData + (keypoint->y * integralImage->widthStep) + keypoint->x);
  widthStep = integralImage->widthStep / 4;
  yoffset = keypoint->scale * widthStep;
  ptrAi = ptr - (yoffset + keypoint->scale);
  ptrDi = ptr + (yoffset + keypoint->scale);
  ptrBi = ptr - (yoffset - keypoint->scale);
  ptrCi = ptr + (yoffset - keypoint->scale);
  valin = *ptrDi - *ptrBi - *ptrCi + *ptrAi;
  ptrAo = ptr - ((yoffset + keypoint->scale) << 1);
  ptrDo = ptr + ((yoffset + keypoint->scale) << 1);
  ptrBo = ptr - ((yoffset - keypoint->scale) << 1);
  ptrCo = ptr + ((yoffset - keypoint->scale) << 1);
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

CamKeypointsMatch	*cam_keypoints_tracking(CamTrackingContext *tc, CamImage *image, int (*roiDetector)(CamImage*, CamKeypoints*, int, int), int (*detectorValue)(CamImage *, CamKeypointShort*))
{
  unsigned int		*seedsIndexes;
  CamKeypointsMatch	*matches;
  int			i;
  int			j;
  int			k;
  int			seedX;
  int			seedY;
  int			seedScale;
  int			fullScale;
  CamKeypoints		seedMatch;
  CamROI		roi;
  CamImage		integralImage;
  CamKeypointShort	*maxOnEachScale;
  CamKeypointShort	currentPoint;

   // initialisation of the tracker
  if (!(tc->previousFeatures.bag))
    {
      srandom(time(NULL)); // used to extract seeds
      (*roiDetector)(image, &tc->previousFeatures, tc->nbFeatures, 0);
      return (NULL);
    }
  // initialisation already done
  else
    {
      seedsIndexes = cam_keypoints_tracking_extract_seeds(tc);

      /* begin integral image computing */
      integralImage.imageData = NULL;
      camIntegralImage(image, &integralImage);
      /* end integral image computing */

      // in the specified research volume, finds the more coherent point corresponding to the current seed <=> highest detector value
      fullScale = (tc->rv.scale << 1) + 1;
      for (i = 0 ; i < tc->nbSeeds ; ++i)
	{
	  seedX = tc->previousFeatures.keypoint[seedsIndexes[i]]->x;
	  seedY = tc->previousFeatures.keypoint[seedsIndexes[i]]->y;
	  seedScale = tc->previousFeatures.keypoint[seedsIndexes[i]]->scale;
	  maxOnEachScale = (CamKeypointShort*)calloc(fullScale, sizeof(CamKeypointShort));
	  // on each scale, get the highest detector's value
	  for (currentPoint.scale = max(seedScale - tc->rv.scale,1), j = 0 ; currentPoint.scale <= seedScale + tc->rv.scale ; ++currentPoint.scale, ++j)
	    {
	      for (currentPoint.y = max(seedY - tc->rv.height, 0) ; currentPoint.y <= min(seedY + tc->rv.height, integralImage.height - 1) ; ++currentPoint.y)
		{
		  for (currentPoint.x = max(seedX - tc->rv.width, 0) ; currentPoint.x <= min(seedX + tc->rv.width, integralImage.width - 1) ; ++currentPoint.x)
		    {
		      currentPoint.value = abs(((*detectorValue)(&integralImage, &currentPoint)));
		      if (currentPoint.value > maxOnEachScale[j].value)
			memcpy(&maxOnEachScale[j], &currentPoint, sizeof(CamKeypointShort));
		    }
		}
	    }
	  for (k = 0 ; k < j ; ++k)
	    {
	      maxOnEachScale[k].value = (maxOnEachScale[k].value << 4) / (maxOnEachScale[k].scale * maxOnEachScale[k].scale);
	    }
	  free (maxOnEachScale);
	}
      free (seedsIndexes);
    }
}

int camKeypointsDetector(CamImage *source, CamKeypoints *points, int nb_max_keypoints, int options);
int camKeypointsRecursiveDetector(CamImage *source, CamKeypoints *points, int nb_max_keypoints, int options);

void			test_cam_keypoints_tracking()
{
  CamTrackingContext	tc;
  CamImage		modelImage;
  CamImage		firstImage;
  CamImage		secondImage;

  /* begin images building */
  modelImage.imageData = NULL;
  camLoadBMP(&modelImage, "resources/photos/scene1.bmp");
  camAllocateYUVImage(&firstImage, modelImage.width, modelImage.height);
  camRGB2YUV(&modelImage, &firstImage);
  camLoadBMP(&modelImage, "resources/photos/scene1.bmp");
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
  camAllocateKeypoints(&tc.previousFeatures, 100);
  /* end tracking context */

  cam_keypoints_tracking(&tc, &firstImage, camKeypointsRecursiveDetector, cam_keypoints_tracking_compute_detector);
  cam_keypoints_tracking(&tc, &secondImage, camKeypointsRecursiveDetector, cam_keypoints_tracking_compute_detector);
}
