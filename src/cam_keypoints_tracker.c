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
#include "camellia.h"

typedef struct
{
  int x;	// dx reseach size
  int y;	// dy reseach size
  int s;	// scale reseach size
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

CamKeypointsMatch	*cam_keypoints_tracking(CamTrackingContext *tc, CamImage *image, int (*detector)(CamImage*, CamKeypoints*, int, int))
{
  unsigned int		*seedsIndexes;
  CamKeypointsMatch	*matches;

   // initialisation of the tracker
  if (!(tc->previousFeatures.bag))
    {
      srandom(time(NULL)); // used to extract seeds
      (*detector)(image, &tc->previousFeatures, tc->nbFeatures, 0);
      return (NULL);
    }
  // initialisation already done
  else
    {
      seedsIndexes = cam_keypoints_tracking_extract_seeds(tc);
      
      free (seedsIndexes);
    }
}

int camKeypointsDetector(CamImage *source, CamKeypoints *points, int nb_max_keypoints, int options);

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
  camLoadBMP(&modelImage, "resources/photos/scene4.bmp");
  camAllocateYUVImage(&secondImage, modelImage.width, modelImage.height);
  camRGB2YUV(&modelImage, &secondImage);
  /* end images building */

  /* begin tracking context */
  tc.nbFeatures = 100;
  tc.nbFrames = 3;
  tc.nbSeeds = 10;
  tc.rv.x = 3;
  tc.rv.y = 3;
  tc.rv.s = 3;
  tc.rv.ds = 0;
  camAllocateKeypoints(&tc.previousFeatures, 100);
  /* end tracking context */

  cam_keypoints_tracking(&tc, &firstImage, camKeypointsDetector);
  cam_keypoints_tracking(&tc, &secondImage, camKeypointsDetector);
}
