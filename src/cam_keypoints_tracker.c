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
#define	CAM_TRACKING_DEBUG_3
#define	CAM_TRACKING_DEBUG_2
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
  CamImage		*previousImage;
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
  register int		k;
  register int		l;
  register unsigned int	m;
  int			valueMax;
  unsigned int		step;
  int			width;

  width = sqrt(tc->nbFeatures);
  step = width / sqrt(tc->nbSeeds);
  seedsIndexes = (unsigned int*)calloc(tc->nbSeeds, sizeof(unsigned int));
  m = 0;
  for (i = step - 1 ; i < width ; i += step)
    {
      for (j = step - 1 ; j < width && m < tc->nbSeeds ; j += step)
	{
	  valueMax = 0;
	  for (k = -1 ; k < 2 ; ++k)
	    {
	      for (l = -1 ; l < 2 ; ++l)
		{
		  if (tc->previousFeatures->keypoint[(i + k) * width + (j + l)] &&
		      abs(tc->previousFeatures->keypoint[(i + k) * width + (j + l)]->value) > valueMax)
		    {
		      valueMax = abs(tc->previousFeatures->keypoint[(i + k) * width + (j + l)]->value);
		      seedsIndexes[m] = (i + k) * width + (j + l);
		    }
		}
	    }
#ifdef CAM_TRACKING_DEBUG_3
	  printf("Seed %i is better in %i with value %i\n", i * width + j, seedsIndexes[m], tc->previousFeatures->keypoint[seedsIndexes[m]]->value);
#endif
	  ++m;
	}
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
	      currentPoint.value = ((*detectorValue)(integralImage, &currentPoint));
	      if (abs(currentPoint.value) > abs(maxOnEachScale[j].value))
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

void		cam_keypoints_tracking_compute_intensity_difference(CamTrackingContext *tc, CamImage *image, int index, CamImage *intensityDifference)
{
  register int	i;
  register int	j;
  int		x;
  int		y;
  int		hh;
  int		hw;

  x = tc->previousFeatures->keypoint[index]->x;
  y = tc->previousFeatures->keypoint[index]->y;
  hh = intensityDifference->height / 2;
  hw = intensityDifference->width / 2;
  for (j = -hh ; j <= hh ; ++j)
    {
      for (i = -hw ; i <= hw ; ++i)
	{
	  intensityDifference->imageData[(j + hh) * intensityDifference->width + (i + hw)] = tc->previousImage->imageData[(y + j) * tc->previousImage->widthStep + (x + i)] - image->imageData[(y + j) * image->widthStep + (x + i)];
	}
    }
}

void	cam_keypoints_tracking_solve_displacement(float gxx, float gxy, float gyy, float ex, float ey, float *dx, float *dy)
{
  float det;

  det = gxx * gyy - gxy * gxy;
  *dx = (gyy * ex - gxy * ey) / det;
  *dy = (gxx * ey - gxy * ex) / det;
}

void		cam_keypoints_tracking_convolve_horiz(CamImage *imgIn, float *kernel, CamImage *imgOut)
{
  register int	i;
  register int	j;
  register int	k;
  float		sum;
  int		ncols;
  int		nrows;
  unsigned char	*ptrrow;
  unsigned char	*ptrout;
  unsigned char	*ppp;

  ptrrow = imgIn->imageData;
  ncols = imgIn->width;
  nrows = imgIn->height;
  ptrout = imgOut->imageData;
  for  (j = 0 ; j < nrows ; ++j)
    {
      for (i = 0 ; i < ncols ; ++i)
	{
	  ppp = ptrrow + i;
	  sum = 0.0f;
	  for (k = 0 ; k < imgOut->width ; ++k)	// size of the kernel is equal to width and height
	    {
	      sum += *ppp ++ * kernel[k];
	    }
	  *ptrout++ = (unsigned char)sum;
	}
    }
}

void		cam_keypoints_tracking_convolve_vert(CamImage *imgIn, float *kernel, CamImage *imgOut)
{
  register int	i;
  register int	j;
  register int	k;
  int		ncols;
  int		nrows;
  unsigned char	*ptrcol;
  unsigned char	*ptrout;
  unsigned char	*ppp;
  float		sum;

  ncols = imgIn->width;
  nrows = imgIn->height;
  ptrcol = imgIn->imageData;
  ptrout = imgOut->imageData;
  for (i = 0 ; i < ncols ; ++i)
    {
      for (j = 0 ; j < nrows ; ++j)
	{
	  ppp = ptrcol + ncols * j;
	  sum = 0.0f;
	  for (k = 0 ; k < imgOut->width ; ++k)	// size of the kernel is equal to width and height
	    {
	      sum += *ppp * kernel[k];
	      ppp += ncols;
	    }
	  *ptrout = (unsigned char)sum;
	  ptrout += ncols;
	}
      ptrcol++;
      ptrout -= nrows * ncols - 1;
    }
}

void		cam_keypoints_tracking_convolve_separate(CamImage *image, float *horizKernel, float *vertKernel, CamImage *outImg)
{
  CamImage	tmpImage;

  camAllocateYUVImage(&tmpImage, image->width, image->height);
  cam_keypoints_tracking_convolve_horiz(image, horizKernel, &tmpImage);
  cam_keypoints_tracking_convolve_vert(&tmpImage, vertKernel, outImg);
  camDeallocateImage(&tmpImage);
}

void	cam_keypoints_tracking_compute_gradients(CamImage *img, CamImage *gradX, CamImage *gradY, float *gaussKernel, float *gaussDerivKernel)
{
  cam_keypoints_tracking_convolve_separate(img, gaussDerivKernel, gaussKernel, gradX);
  cam_keypoints_tracking_convolve_separate(img, gaussKernel, gaussDerivKernel, gradY);
}

void		cam_keypoints_tracking_compute_kernels(float sigma, float *gaussKernel, float *gaussDerivKernel, int windowSize)
{
  register int	i;
  int		hw;

  i = 0;
  hw = windowSize / 2;
  for (i = -hw ; i < hw ; ++i)
    {
      gaussKernel[i + hw] = exp(-i*i / (2 * sigma * sigma));
      gaussDerivKernel[i + hw] = -i / (sigma * sigma) * gaussKernel[i + hw];
    }
}

void				cam_keypoints_tracking_compute_gradients_sum(CamImage *gradX1, CamImage *gradY1, CamImage *gradX2, CamImage *gradY2,
						     int x1, int y1, int x2, int y2, int windowSize, CamImage *gradX, CamImage *gradY)
{
  register int			i;
  register int			j;
  register unsigned char	*pGradX;
  register unsigned char	*pGradY;
  unsigned char			g1;
  unsigned char			g2;

  pGradX = gradX->imageData;
  pGradY = gradY->imageData;
  for (j = -windowSize / 2 ; j <= windowSize / 2 ; ++j)
    {
      for (i = -windowSize / 2 ; i <= windowSize / 2 ; ++i)
	{
	  g1 = gradX1->imageData[(y1 + j) * gradX1->width + (x1 + i)];
	  g2 = gradX2->imageData[(y2 + j) * gradX2->width + (x2 + i)];
	  *pGradX++ = g2 - g1;
	  g1 = gradY1->imageData[(y1 + j) * gradY1->width + (x1 + i)];
	  g2 = gradY2->imageData[(y2 + j) * gradY2->width + (x2 + i)];
	  *pGradY++ = g2 - g1;
	}
    }
}

void				cam_keypoints_tracking_compute_2by2_gradient_matrix(CamImage *gradX, CamImage *gradY, int windowSize, int *gxx, int *gxy, int *gyy)
{
  register unsigned char	*pGradX;
  register unsigned char	*pGradY;
  register int			gx;
  register int			gy;
  register int			i;

  pGradX = gradX->imageData;
  pGradY = gradY->imageData;
  *gxx = 0;
  *gxy = 0;
  *gyy = 0;
  for (i = 0 ; i < windowSize * windowSize ; ++i)
    {
      gx = *pGradX++;
      gy = *pGradY++;
      *gxx += gx * gx;
      *gyy += gy * gy;
      *gxy += gx * gy;
    }
}

void		cam_keypoints_tracking_compute_2by1_error_vector(CamImage *intensityImage, CamImage *gradX, CamImage *gradY, int *ex, int *ey)
{
  register int	i;
  unsigned char	*ppp;
  unsigned char	*gradx;
  unsigned char	*grady;
  int		diff;

  *ex = 0;
  *ey = 0;
  ppp = intensityImage->imageData;
  gradx = gradX->imageData;
  grady = gradY->imageData;
  for (i = 0 ; i < intensityImage->width * intensityImage->height ; ++i)
    {
      diff = *ppp++;
      *ex += diff * (*gradx++);
      *ey += diff * (*grady++); 
    }
}

void	cam_keypoints_tracking_solve_displacement_equation(int gxx, int gxy, int gyy, int ex, int ey, int *shiftX, int *shiftY)
{
  int	det;

  det = gxx * gyy - gxy * gxy;
  if (det * det < 0.01)
    {
      *shiftX = 0;
      *shiftY = 0;
    }
  else
    {
      *shiftX = (gyy * ex - gxy *ey) / det;
      *shiftY = (gxx * ey - gxy *ey) / det;
    }
}

void		cam_keypoints_tracking_compute_local_image_displacement(CamTrackingContext *tc, CamImage *image, int index, int *shiftX, int *shiftY,
									int windowSize, float *gaussKernel, float *gaussDerivKernel)
{
  CamImage	intensityDifference;
  CamImage	gradX1;
  CamImage	gradY1;
  CamImage	gradX2;
  CamImage	gradY2;
  CamImage	gradX;
  CamImage	gradY;
  CamImage	img;
  int		x1;
  int		y1;
  int		x2;
  int		y2;
  int		gxx;
  int		gxy;
  int		gyy;
  int		ex;
  int		ey;

  camAllocateYUVImage(&intensityDifference, windowSize, windowSize);
  camAllocateYUVImage(&gradX, windowSize, windowSize);
  camAllocateYUVImage(&gradY, windowSize, windowSize);
  camAllocateYUVImage(&gradX1, windowSize, windowSize);
  camAllocateYUVImage(&gradY1, windowSize, windowSize);
  camAllocateYUVImage(&gradX2, windowSize, windowSize);
  camAllocateYUVImage(&gradY2, windowSize, windowSize);
  camAllocateYUVImage(&img, windowSize, windowSize);
  x1 = (tc->previousFeatures->keypoint[index]->x - windowSize / 2);
  y1 = (tc->previousFeatures->keypoint[index]->y - windowSize / 2);

  memcpy(img.imageData, &tc->previousImage->imageData[y1 * image->widthStep + x1], windowSize * windowSize * sizeof(unsigned char)); // fix me to be not necessary
  cam_keypoints_tracking_compute_gradients(&img, &gradX1, &gradY1, gaussKernel, gaussDerivKernel);
  memcpy(img.imageData, &image->imageData[y1 * image->widthStep + x1], windowSize * windowSize * sizeof(unsigned char)); // fix me to be not necessary
  cam_keypoints_tracking_compute_gradients(&img, &gradX2, &gradY2, gaussKernel, gaussDerivKernel);

  x1 = windowSize / 2;
  y1 = windowSize / 2;
  x2 = x1;
  y2 = y1;
  cam_keypoints_tracking_compute_intensity_difference(tc, image, index, &intensityDifference);
  cam_keypoints_tracking_compute_gradients_sum(&gradX1, &gradY1, &gradX2, &gradY2, x1, y1, x2, y2, windowSize, &gradX, &gradY);
  cam_keypoints_tracking_compute_2by2_gradient_matrix(&gradX, &gradY, windowSize, &gxx, &gxy, &gyy);
  cam_keypoints_tracking_compute_2by1_error_vector(&intensityDifference, &gradX, &gradY, &ex, &ey);
  cam_keypoints_tracking_solve_displacement_equation(gxx, gxy, gyy, ex, ey, shiftX, shiftY);

#ifdef CAM_TRACKING_DEBUG_1
  printf("shiftX : %i shiftY : %i\n", *shiftX, *shiftY);
#endif

  camDeallocateImage(&intensityDifference);
  camDeallocateImage(&gradX);
  camDeallocateImage(&gradY);
  camDeallocateImage(&gradX1);
  camDeallocateImage(&gradY1);
  camDeallocateImage(&gradX2);
  camDeallocateImage(&gradY2);
  camDeallocateImage(&img);
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
  int			shiftX;
  int			shiftY;
  int			windowSize;
  float			*gaussKernel;
  float			*gaussDerivKernel;
#ifdef CAM_TRACKING_SUBTIMINGS
  int			deltaTimers1;
  struct timeval	tv1;
  struct timeval	tv2;
  int			deltaTimers2;
  struct timeval	tv3;
  struct timeval	tv4;
#endif

  seedsMatches = (CamKeypointsMatches*)malloc(sizeof(*seedsMatches));
  camAllocateKeypointsMatches(seedsMatches, tc->nbSeeds);
  
  // in the specified research volume, finds the more coherent point corresponding to the current seed <=> highest detector value
  l = 0;
  windowSize = min(tc->rw.height, tc->rw.width) / 3;
  if (!(windowSize % 2))
    windowSize--;
  gaussKernel = (float *)malloc(windowSize * sizeof(float));
  gaussDerivKernel = (float *)malloc(windowSize * sizeof(float));
#ifdef CAM_TRACKING_SUBTIMINGS
  gettimeofday(&tv3, NULL);      
#endif
  cam_keypoints_tracking_compute_kernels(1, gaussKernel, gaussDerivKernel, windowSize);
#ifdef CAM_TRACKING_SUBTIMINGS
  gettimeofday(&tv4, NULL);      
  deltaTimers2 = tv4.tv_usec - tv3.tv_usec;
  if (deltaTimers2 < 0)
    deltaTimers2 = 1000000 + deltaTimers2;
  printf("Kernels computation : %ius\n", deltaTimers2);
#endif
  for (i = 0 ; i < tc->nbSeeds ; ++i)
    {
      cam_keypoints_tracking_compute_local_image_displacement(tc, image, seedsIndexes[i], &shiftX, &shiftY, windowSize, gaussKernel, gaussDerivKernel);
      maxOnEachScale = cam_keypoints_tracking_extract_overall_max_on_each_scale(tc, integralImage, detectorValue, seedsIndexes[i], shiftX, shiftY, 0);
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
      deltaTimers1 = tv2.tv_usec - tv1.tv_usec;
      if (deltaTimers1 < 0)
	deltaTimers1 = 1000000 + deltaTimers1;
      printf("Feature description : %ius\n", deltaTimers1);
#endif
      /* end computation of descriptor */
    }
  free(gaussKernel);
  free(gaussDerivKernel);
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
      if (!tc->previousFeatures->keypoint[i])
	continue ;
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
#ifdef	CAM_TRACKING_DEBUG_2
	  printf("Missmatch : %i %i %i %i, %i %i %i %i\n", featuresMatches->pairs[i].p1->x, featuresMatches->pairs[i].p1->y, featuresMatches->pairs[i].p1->scale, featuresMatches->pairs[i].p1->value, featuresMatches->pairs[i].p2->x,featuresMatches->pairs[i].p2->y,featuresMatches->pairs[i].p2->scale, featuresMatches->pairs[i].p2->value);
#endif
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

#ifdef CAM_TRACKING_DEBUG_3
void		printMatchings(CamKeypointsMatches *matches)
{
  register int	i;
  
  for (i = 0 ; i < matches->nbMatches ; ++i)
    {
      printf("x1: %i\ty1: %i\tscale1: %i\tvalue: %i\t// x2: %i\ty2: %i\tscale2: %i\tvalue: %i\n", matches->pairs[i].p1->x, matches->pairs[i].p1->y, matches->pairs[i].p1->scale, matches->pairs[i].p1->value, matches->pairs[i].p2->x, matches->pairs[i].p2->y, matches->pairs[i].p2->scale, matches->pairs[i].p2->value);
    }
}
#endif

CamKeypointsMatches	*cam_keypoints_tracking(CamTrackingContext *tc, CamImage *image, int (*roiDetector)(CamImage *,CamImage *, CamKeypoints *, int, int), int (*detectorValue)(CamImage *, CamKeypointShort*), int options)
{
  CamImage		*integralImage;
  CamImage		*previousIntegral;
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
      integralImage = (CamImage*)malloc(sizeof(CamImage));
      integralImage->imageData = NULL;
      camSetMaxROI(&roix, image);
      roix.coi = 1;
      image->roi = &roix;
      camIntegralImage(image, integralImage);
      integralImage->roi = &roi;

      // configuration of the window related to the seeding
      image->roi = &roi;
      image->roi->width = tc->rw.width;
      image->roi->height = tc->rw.height;
      i = 0;
#ifdef CAM_TRACKING_SUBTIMINGS
      timeInROIDetector = 0;
#endif
      for (image->roi->xOffset = 0 ; image->roi->xOffset < image->width ; image->roi->xOffset += tc->rw.width)
	{
	  for (image->roi->yOffset = 0 ; image->roi->yOffset < image->height ; image->roi->yOffset += tc->rw.height)
	    {	
#ifdef CAM_TRACKING_SUBTIMINGS
	      t1 = camGetTimeMs();
#endif
	      (*roiDetector)(image, integralImage, &tmpPoint, 1, options);
#ifdef CAM_TRACKING_SUBTIMINGS
	      t2 = camGetTimeMs();
	      timeInROIDetector += t2 - t1;
#endif
	      if (tmpPoint.nbPoints)
		{
		  tc->previousFeatures->keypoint[i] = &tc->previousFeatures->bag[i];
		  memcpy(tc->previousFeatures->keypoint[i], tmpPoint.keypoint[0], sizeof(CamKeypoint));
		  tc->previousFeatures->nbPoints++;
		}
	      else
		{
#ifdef CAM_TRACKING_DEBUG_3
		  printf("Index %i has no keypoint\n", i);
#endif
		  tc->previousFeatures->keypoint[i] = NULL;
		}
	      ++i;
	    }
	}
#ifdef CAM_TRACKING_SUBTIMINGS
      printf("Time spent in roiDetector : %ims\n", timeInROIDetector);
#endif
      /*      camDeallocateImage(integralImage);
      free(integralImage);
      */
      previousIntegral = integralImage;
      camFreeKeypoints(&tmpPoint);
      tc->previousImage = image;
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

      /* begin integral image computing */ 
#ifdef CAM_TRACKING_TIMINGS
      t4 = camGetTimeMs();
#endif
      integralImage = (CamImage*)malloc(sizeof(CamImage));
      integralImage->imageData = NULL;

      camSetMaxROI(&roix, image);
      roix.coi = 1;
      image->roi = &roix;
      camIntegralImage(image, integralImage);
      integralImage->roi = &roi;

#ifdef CAM_TRACKING_TIMINGS
      t5 = camGetTimeMs();
#endif
     /* end integral image computing */
#ifdef CAM_TRACKING_TIMINGS
      t1 = camGetTimeMs();
#endif
      seedsMatches = cam_keypoints_tracking_extract_seed_matches(tc, image, integralImage, detectorValue, currentFeatures, seedsIndexes, options); // need to check the descriptor computation
#ifdef CAM_TRACKING_TIMINGS
      t2 = camGetTimeMs();
#endif
      keypointsMatches = cam_keypoints_tracking_extract_points_matching(tc, image, integralImage, detectorValue, seedsMatches, seedsIndexes, currentFeatures, options);
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

#ifdef CAM_TRACKING_DEBUG_1
      printf("Seeds matches inliers : %i outliers : %i\n", seedsMatches->nbMatches, seedsMatches->nbOutliers);
#endif
#ifdef CAM_TRACKING_DEBUG_3
      printf("Matches of the seeds :\n");
      printMatchings(seedsMatches);
#endif
#ifdef CAM_TRACKING_DEBUG_1
      printf("Keypoints matches inliers : %i outliers : %i\n", keypointsMatches->nbMatches, keypointsMatches->nbOutliers);
#endif
#ifdef CAM_TRACKING_DEBUG_3
      printf("Matches of the keypoints :\n");
      printMatchings(keypointsMatches);
#endif

      res->nbMatches = seedsMatches->nbMatches + keypointsMatches->nbMatches;
      res->nbOutliers = tc->nbFeatures - seedsMatches->nbMatches + keypointsMatches->nbMatches;

      // cam_keypoints_tracking_clean_keypoints(currentFeatures, keypointsMatches, seedsMatches->nbMatches); => to be done
      cam_keypoints_tracking_fill_empty_area();
      
      camDeallocateImage(tc->previousImage);
      tc->previousImage = image;

      camDeallocateImage(integralImage);
      free(integralImage);

      free (seedsIndexes);
      camFreeKeypoints(tc->previousFeatures);
      free(tc->previousFeatures);
      camFreeKeypointsMatches(seedsMatches);
      free (seedsMatches);

      tc->previousFeatures = currentFeatures;
      camFreeKeypointsMatches(keypointsMatches);
      free (keypointsMatches);
      return (res);
    }
}

void	cam_keypoint_tracking_configure_context(CamTrackingContext *tc, int nbFeatures, int height, int width)
{
  tc->nbFeatures = nbFeatures;
  tc->nbFrames = 3;
  tc->nbSeeds = sqrt(tc->nbFeatures) - 1;
  tc->rv.width = 2;
  tc->rv.height = 2;
  tc->rv.scale = 2;
  tc->rv.ds = 0;
  tc->rw.height = height / sqrt(tc->nbFeatures);
  tc->rw.width = width / sqrt(tc->nbFeatures);
  tc->previousFeatures = NULL;
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
  char			img2[] = "./resources/rover/rotation1.bmp";

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

  cam_keypoint_tracking_configure_context(&tc, 100, firstImage.height, firstImage.width);

#ifdef CAM_TRACKING_TIMINGS
  t1 = camGetTimeMs();
#endif
  
  cam_keypoints_tracking(&tc, &firstImage, camKeypointsRecursiveDetector, cam_keypoints_tracking_compute_detector, CAM_UPRIGHT);

#ifdef CAM_TRACKING_TIMINGS
  t2 = camGetTimeMs();
#endif

  track = cam_keypoints_tracking(&tc, &secondImage, camKeypointsRecursiveDetector, cam_keypoints_tracking_compute_detector, CAM_UPRIGHT);

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
