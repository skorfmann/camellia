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

#ifndef LINUX
#define inline
#endif

#include <stdio.h>	// printf
#include <math.h>	// fabs, exp
#include <emmintrin.h>	// _mm_malloc
#include <stdlib.h>	// malloc, qsort
#include <limits.h>	// INT_MAX
#include <string.h>	// memcpy
#include "camellia.h"

#define max(a, b) (a > b ? a : b)
#define min(a, b) (a > b ? b : a)

#define MAX_KERNEL_WIDTH 71

#define CAM_TRACKING2_TIMINGS

typedef enum
  {
    TRUE,
    FALSE
  } BOOL;

typedef enum
  {
    TRACKED,
    MAX_ITERATIONS,
    OOB,
    SMALL_DET,
    LARGE_RESIDUE,
    NOT_YET_EVALUATED
  } TRACKING_STATUS;

typedef struct
{
  int	width;
  float	data[MAX_KERNEL_WIDTH];
} CamConvolutionKernel;

typedef	struct
{
  int	width;	// half width reseach size
  int	height;	// half height reseach size
  int	scale;	// half scale reseach size
}	researchVolume;

typedef struct		s_camList
{
  void			*data;
  struct s_camList	*next;
  int			index;
}			CamList;

typedef	struct
{
  int	height;
  int	width;
}	researchWindow;

typedef struct
{
  int	ncols;
  int	nrows;
  float	*data;
}	CamFloatImage;

typedef struct
{
  CamFloatImage	gradX;
  CamFloatImage	gradY;
  CamFloatImage	image;
  float		x;
  float		y;
}		pyramidLevel;

typedef struct
{
  pyramidLevel	*img1;
  pyramidLevel	*img2;
  int		scale;
}		pyramidLevels;

typedef struct
{
  int		nbLevels;
  pyramidLevels *levels;
}		pyramid;

typedef struct
{
  int	borderX;
  int	borderY;
}	borders;

typedef struct
{
  float	x;
  float	y;
}	vector;

typedef struct
{
  CamKeypoints		*previousFeatures;
  CamKeypoints		*previousCorners;
  CamImage		*previousImage;
  CamImage		*previousIntegralImage;
  int			nbFeatures;
  int			nbDetectedFeatures;
  researchVolume	rv;
  researchWindow	rw;
  borders		b;
  pyramid		*pyramidImages;
  CamConvolutionKernel	gaussianKernel;
  CamConvolutionKernel	gaussianDerivedKernel;
}			CamTrackingContext;

inline int	cam_keypoints_tracking2_compute_detector(CamImage *integralImage, CamKeypointShort *keypoint)
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

void		cam_keypoints_tracking_free_linked_list(CamList *l)
{
  CamList	*ptr;

  ptr = l;
  while (ptr)
    {
      l = l->next;
      free (ptr);
      ptr = l;
    }
}

inline CamList	*cam_keypoints_tracking_add_to_linked_list(CamList *l, void *data)
{
  CamList	*head;

  head = (CamList*)malloc(sizeof(CamList));
  head->data = data;
  head->next = l;
  if (l)
    head->index = l->index + 1;
  else
    head->index = 1;
  return (head);
}

void	CamAllocateFloatImage(CamFloatImage *res, int ncols, int nrows)
{
  res->nrows = nrows;
  res->ncols = ncols;
  res->data = (float *)malloc(ncols * nrows * sizeof(float));
}

void	CamDisallocateFloatImage(CamFloatImage *res)
{
  res->ncols = 0;
  res->nrows = 0;
  if (res->data)
    {
      free (res->data);
      res->data = NULL;
    }
}

void		cam_keypoints_tracking2_compute_kernels(float sigma, CamConvolutionKernel *gauss, CamConvolutionKernel *gaussderiv)
{
  const float	factor = 0.01f;
  int i;

  {
    const int hw = MAX_KERNEL_WIDTH / 2;
    float max_gauss = 1.0f, max_gaussderiv = (float) (sigma*exp(-0.5f));
    
    for (i = -hw ; i <= hw ; i++)
      {
	gauss->data[i+hw] = (float) (exp(-i*i / (2*sigma*sigma)));
	gaussderiv->data[i+hw] = -i * gauss->data[i+hw];
      }
    gauss->width = MAX_KERNEL_WIDTH;
    for (i = -hw ; fabs(gauss->data[i+hw] / max_gauss) < factor ; i++, gauss->width -= 2);
    gaussderiv->width = MAX_KERNEL_WIDTH;
    for (i = -hw ; fabs(gaussderiv->data[i+hw] / max_gaussderiv) < factor ; i++, gaussderiv->width -= 2);
    if (gauss->width == MAX_KERNEL_WIDTH || 
        gaussderiv->width == MAX_KERNEL_WIDTH)
      camError("cam_keypoints_tracking_compute_kernels", "Max kernel width too small\n");
  }

  for (i = 0 ; i < gauss->width ; i++)
    gauss->data[i] = gauss->data[i+(MAX_KERNEL_WIDTH-gauss->width)/2];
  for (i = 0 ; i < gaussderiv->width ; i++)
    gaussderiv->data[i] = gaussderiv->data[i+(MAX_KERNEL_WIDTH-gaussderiv->width)/2];
  {
    const int hw = gaussderiv->width / 2;
    float den;
    
    den = 0.0;
    for (i = 0 ; i < gauss->width ; i++)  den += gauss->data[i];
    for (i = 0 ; i < gauss->width ; i++)  gauss->data[i] /= den;
    den = 0.0;
    for (i = -hw ; i <= hw ; i++)  den -= i*gaussderiv->data[i+hw];
    for (i = -hw ; i <= hw ; i++)  gaussderiv->data[i+hw] /= den;
  }
}

void		cam_keypoints_tracking_free_context(CamTrackingContext *tc)
{
  register int	i;

  if (tc->previousCorners)
    {
      camFreeKeypoints(tc->previousCorners);
      free(tc->previousCorners);
      tc->previousCorners = NULL;
    }
  if (tc->previousFeatures)
    {
      camFreeKeypoints(tc->previousFeatures);
      free(tc->previousFeatures);
      tc->previousFeatures = NULL;
    }
  if (tc->previousImage)
    {
      camDeallocateImage(tc->previousImage);
      tc->previousImage = NULL;
    }
  if (tc->previousIntegralImage)
    {
      camDeallocateImage(tc->previousIntegralImage);
      free(tc->previousIntegralImage);
      tc->previousIntegralImage = NULL;
    }
  if (tc->pyramidImages)
    {
      for (i = 0 ; i < tc->pyramidImages->nbLevels ; ++i)
	{
	  CamDisallocateFloatImage(&tc->pyramidImages->levels[i].img1->gradX);
	  CamDisallocateFloatImage(&tc->pyramidImages->levels[i].img1->gradY);
	  CamDisallocateFloatImage(&tc->pyramidImages->levels[i].img1->image);
	  free(tc->pyramidImages->levels[i].img1);
	  CamDisallocateFloatImage(&tc->pyramidImages->levels[i].img2->gradX);
	  CamDisallocateFloatImage(&tc->pyramidImages->levels[i].img2->gradY);
	  CamDisallocateFloatImage(&tc->pyramidImages->levels[i].img2->image);
	  free(tc->pyramidImages->levels[i].img2);
	}
      free(tc->pyramidImages->levels);
      free(tc->pyramidImages);
      tc->pyramidImages = NULL;
    }
  tc->nbFeatures = 0;
  tc->nbDetectedFeatures = 0;
  tc->rv.width = 0;
  tc->rv.height = 0;
  tc->rv.scale = 0;
  tc->rw.height = 0;
  tc->rw.width = 0;
}

void			cam_keypoint_tracking2_configure_context(CamTrackingContext *tc, int nbFeatures, int rwHeight, int rwWidth, int rvHeight, int rvWidth, int rvScale, CamList *scales, CamImage *image)
{
  register CamList	*ptr;
  register int		i;

  if (!tc)
    return ;
  cam_keypoints_tracking2_compute_kernels(1.0f, &tc->gaussianKernel, &tc->gaussianDerivedKernel);
  tc->nbFeatures = nbFeatures;
  tc->nbDetectedFeatures = 0;
  tc->rv.width = rvWidth;
  tc->rv.height = rvHeight;
  tc->rv.scale = rvScale;
  tc->rw.height = rwHeight;
  tc->rw.width = rwWidth;
  tc->b.borderX = max(max(tc->rw.height, tc->rw.width), min(image->height / 20, image->width / 20));
  tc->b.borderY = max(max(tc->rw.height, tc->rw.width), min(image->height / 20, image->width / 20));
  tc->previousCorners = NULL;
  tc->previousFeatures = NULL;
  tc->previousImage = NULL;
  tc->previousIntegralImage = NULL;
  tc->pyramidImages = NULL;
  tc->previousCorners = (CamKeypoints*)malloc(sizeof(CamKeypoints));
  tc->previousFeatures = (CamKeypoints*)malloc(sizeof(CamKeypoints));
  camAllocateKeypoints(tc->previousCorners, tc->nbFeatures);
  camAllocateKeypoints(tc->previousFeatures, tc->nbFeatures);
#ifdef __SSE2__
  tc->previousCorners->bag = (CamKeypoint*)_mm_malloc(sizeof(CamKeypoint) * tc->nbFeatures, 16);
  tc->previousFeatures->bag = (CamKeypoint*)_mm_malloc(sizeof(CamKeypoint) * tc->nbFeatures, 16);
#else
  tc->previousCorners->bag = (CamKeypoint*)malloc(sizeof(CamKeypoint) * tc->nbFeatures);
  tc->previousFeatures->bag = (CamKeypoint*)malloc(sizeof(CamKeypoint) * tc->nbFeatures);
#endif
  if (scales)
    {
      ptr = scales;
      tc->pyramidImages = (pyramid*)malloc(sizeof(pyramid));
      tc->pyramidImages->nbLevels = scales->index;
      tc->pyramidImages->levels = (pyramidLevels*)malloc(scales->index * sizeof(pyramidLevels));
      for (ptr = scales, i = 0 ; ptr ; ptr = ptr->next, ++i)
	{
	  tc->pyramidImages->levels[i].scale = (int)ptr->data;
	  tc->pyramidImages->levels[i].img1 = (pyramidLevel*)malloc(sizeof(pyramidLevel));
	  CamAllocateFloatImage(&tc->pyramidImages->levels[i].img1->gradX, image->width / tc->pyramidImages->levels[i].scale, image->height / tc->pyramidImages->levels[i].scale);
	  CamAllocateFloatImage(&tc->pyramidImages->levels[i].img1->gradY, image->width / tc->pyramidImages->levels[i].scale, image->height / tc->pyramidImages->levels[i].scale);
	  CamAllocateFloatImage(&tc->pyramidImages->levels[i].img1->image, image->width / tc->pyramidImages->levels[i].scale, image->height / tc->pyramidImages->levels[i].scale);
	  tc->pyramidImages->levels[i].img2 = (pyramidLevel*)malloc(sizeof(pyramidLevel));
	  CamAllocateFloatImage(&tc->pyramidImages->levels[i].img2->gradX, image->width / tc->pyramidImages->levels[i].scale, image->height / tc->pyramidImages->levels[i].scale);
	  CamAllocateFloatImage(&tc->pyramidImages->levels[i].img2->gradY, image->width / tc->pyramidImages->levels[i].scale, image->height / tc->pyramidImages->levels[i].scale);
	  CamAllocateFloatImage(&tc->pyramidImages->levels[i].img2->image, image->width / tc->pyramidImages->levels[i].scale, image->height / tc->pyramidImages->levels[i].scale);
	}
    }
}

void			cam_keypoints_tracking2_convolve_horiz(CamFloatImage *imgin, CamConvolutionKernel *kernel, CamFloatImage *imgout)
{
  float			*ptrrow;
  register float	*ptrout;
  register float	*ppp;
  register float	sum;
  register int		radius;
  register int		ncols;
  register int		nrows;
  register int		i, j, k;

  ptrrow = imgin->data;
  ptrout = imgout->data;
  radius = kernel->width / 2;
  ncols = imgin->ncols;
  nrows = imgin->nrows;
  for (j = 0 ; j < nrows ; j++)
    {
      
      for (i = 0 ; i < radius ; i++)
	*ptrout++ = 0.0;
      
      for ( ; i < ncols - radius ; i++)
	{
	  ppp = ptrrow + i - radius;
	  sum = 0.0;
	  for (k = kernel->width-1 ; k >= 0 ; k--)
	    sum += *ppp++ * kernel->data[k];
	  *ptrout++ = sum;
	}

      for ( ; i < ncols ; i++)
	*ptrout++ = 0.0;
      
      ptrrow += ncols;
    }
}

void			cam_keypoints_tracking2_convolve_vert(CamFloatImage *imgin, CamConvolutionKernel *kernel, CamFloatImage *imgout)
{
  float			*ptrcol;
  register float	*ptrout;
  register float	*ppp;
  register float	sum;
  register int		radius;
  register int		ncols;
  register int		nrows;
  register int		i, j, k;

  ptrcol = imgin->data;
  ptrout = imgout->data;
  radius = kernel->width / 2;
  ncols = imgin->ncols;
  nrows = imgin->nrows;
  for (i = 0 ; i < ncols ; i++)
    {
      
      for (j = 0 ; j < radius ; j++)
	{
	  *ptrout = 0.0;
	  ptrout += ncols;
	}
      
      for ( ; j < nrows - radius ; j++)
	{
	  ppp = ptrcol + ncols * (j - radius);
	  sum = 0.0;
	  for (k = kernel->width-1 ; k >= 0 ; k--)
	    {
	      sum += *ppp * kernel->data[k];
	      ppp += ncols;
	    }
	  *ptrout = sum;
	  ptrout += ncols;
	}
      
      for ( ; j < nrows ; j++)
	{
	  *ptrout = 0.0;
	  ptrout += ncols;
	}

      ptrcol++;
      ptrout -= nrows * ncols - 1;
    }
}

void		cam_keypoints_tracking2_convolve_separate(CamFloatImage *imgin, CamConvolutionKernel *horiz_kernel, CamConvolutionKernel *vert_kernel, CamFloatImage *imgout)
{
  CamFloatImage	tmpimg;

  CamAllocateFloatImage(&tmpimg, imgin->ncols, imgin->nrows);
  cam_keypoints_tracking2_convolve_horiz(imgin, horiz_kernel, &tmpimg);
  cam_keypoints_tracking2_convolve_vert(&tmpimg, vert_kernel, imgout);
  CamDisallocateFloatImage(&tmpimg);
}

void	cam_keypoints_tracking2_compute_gradients(CamFloatImage *img, CamFloatImage *gradx, CamFloatImage *grady, CamConvolutionKernel *gauss_kernel, CamConvolutionKernel *gaussderiv_kernel)
{
  cam_keypoints_tracking2_convolve_separate(img, gaussderiv_kernel, gauss_kernel, gradx);
  cam_keypoints_tracking2_convolve_separate(img, gauss_kernel, gaussderiv_kernel, grady);
}

void			cam_keypoints_tracking2_copy_image_to_float_image(CamFloatImage *dst, CamImage *src, int scale)
{
  register unsigned int	i;

  for (i = 0 ; i < dst->ncols * dst->nrows ; ++i)
    dst->data[i] = (float)src->imageData[i];
}

float	cam_keypoints_tracking2_min_eigen_value(float gxx, float gxy, float gyy)
{
  return ((gxx + gyy - sqrt((gxx - gyy) * (gxx - gyy) + 4 * gxy* gxy)) / 2.0f);
}

float	cam_keypoints_tracking2_max_eigen_value(float gxx, float gxy, float gyy)
{
  return ((gxx + gyy + sqrt((gxx - gyy) * (gxx - gyy) + 4 * gxy* gxy)) / 2.0f);
}

float	cam_keypoints_tracking2_compute_main_eigen_vector(float gxx, float gxy, float gyy)
{
  float	eigenValue;
  float	res;

  eigenValue = cam_keypoints_tracking2_max_eigen_value(gxx, gxy, gyy);
  res = (eigenValue - gxx) / gxy;
  return (res);
}

int		cam_keypoints_tracking2_position_filter(CamTrackingContext *tc, CamKeypointShort *pointsList, CamKeypointShort *sortedPointsList, int nbPoints, int minDistance)
{
  register int	x;
  register int	y;
  register int	i;
  int		ncols;
  int		nrows;
  int		found;
  int		borderX;
  int		borderY;

  ncols = tc->pyramidImages->levels[tc->pyramidImages->nbLevels - 1].img1->image.ncols ;
  nrows = tc->pyramidImages->levels[tc->pyramidImages->nbLevels - 1].img1->image.nrows;
  found = 0;
  borderX = tc->b.borderX / tc->pyramidImages->levels[tc->pyramidImages->nbLevels - 1].scale;
  borderY = tc->b.borderY / tc->pyramidImages->levels[tc->pyramidImages->nbLevels - 1].scale;
  minDistance /= tc->pyramidImages->levels[tc->pyramidImages->nbLevels - 1].scale;
  for (i = 0 ; i < (ncols - 2 * borderX) * (nrows - 2 * borderY) && found < tc->nbFeatures ; ++i)
    {
      if (!pointsList[(sortedPointsList[i].y - borderY) * (ncols- 2 * borderX) + (sortedPointsList[i].x - borderX)].value)
	{
	  sortedPointsList[i].value = 0;
	  continue ;
	}
      for (y = max(borderY, sortedPointsList[i].y - minDistance / 2) ; y < min(nrows - borderY, sortedPointsList[i].y + 1 + minDistance / 2) ; ++y)
	{
	  for (x = max(borderX, sortedPointsList[i].x - minDistance / 2) ; x < min(ncols - borderX, sortedPointsList[i].x + 1 + minDistance / 2) ; ++x)
	    {
	      pointsList[(y - borderY) * (ncols - 2 * borderX) + (x - borderX)].value = 0;
	    }
	}
      ++found;
    }
  return (found);
}

BOOL	cam_keypoints_tracking2_is_in_range(CamImage *image, CamKeypointShort *point)
{
  if (point->x - 2 * point->scale >= 0 && point->y - 2 * point->scale && point->x + 2 * point->scale < image->width && point->y + 2 * point->scale < image->height)
    return (TRUE);
  return (FALSE);
}

void			cam_keypoints_tracking2_locate_keypoints(CamTrackingContext *tc, CamKeypoints *corners, CamKeypoints *features)
{
  register int		i;
  int			detectorValue1;
  int			detectorValue2;
  int			previousDetectorValue;
  CamKeypointShort	keypoint1;
  CamKeypointShort	keypoint2;
  CamKeypointShort	keypointMax;
  BOOL			increasingValue;
  BOOL			previousIncreasingValue;
  float			shiftX;
  float			shiftY;
  float			incX;
  float			incY;
  int			scale;

  for (i = 0 ; i < tc->nbDetectedFeatures  ; ++i)
    {
      previousDetectorValue = 0;
      if (abs(corners->keypoint[i]->angle) > 1)
	{
	  incX = 100.0f / (float)corners->keypoint[i]->angle;
	  incY = 1.0f;
	}
      else
	{
	  incX = 1.0f;
	  incY = (float)corners->keypoint[i]->angle / 100.0f;
	}
      scale = 1;
      increasingValue = TRUE;
      previousIncreasingValue = TRUE;
      memcpy(&keypointMax, &corners->keypoint[i]->x, sizeof(CamKeypointShort));
      //printf("%f %f\n", incX, incY);
      while (increasingValue == TRUE || previousIncreasingValue == TRUE)
	{
	  keypoint1.x = corners->keypoint[i]->x + incX * scale;
	  keypoint1.y = corners->keypoint[i]->y + incY * scale;
	  keypoint1.scale = scale;
	  keypoint2.x = corners->keypoint[i]->x - incX * scale;
	  keypoint2.y = corners->keypoint[i]->y - incY * scale;
	  keypoint2.scale = scale;
	  //printf("%i %i %i\n", keypoint1.x, keypoint1.y, keypoint1.scale);
	  //printf("%i %i %i\n", keypoint2.x, keypoint2.y, keypoint2.scale);
	  if (cam_keypoints_tracking2_is_in_range(tc->previousImage, &keypoint1) == TRUE)
	    detectorValue1 = cam_keypoints_tracking2_compute_detector(tc->previousIntegralImage, &keypoint1);
	  else
	    detectorValue1 = 0;
	  if (cam_keypoints_tracking2_is_in_range(tc->previousImage, &keypoint2) == TRUE)
	    detectorValue2 = cam_keypoints_tracking2_compute_detector(tc->previousIntegralImage, &keypoint2);
	  else
	    detectorValue2 = 0;
	  previousIncreasingValue = increasingValue;
	  if (max(abs(detectorValue1), abs(detectorValue2)) / (scale * scale) > abs(previousDetectorValue))
	    increasingValue = TRUE;
	  else
	    increasingValue = FALSE;
	  if (increasingValue == TRUE)
	    {
	      if (abs(detectorValue1) > abs(detectorValue2))
		{
		  previousDetectorValue = detectorValue1 / (scale * scale);
		  memcpy(&keypointMax, &keypoint1, sizeof(CamKeypointShort));
		}
	      else
		{
		  previousDetectorValue = detectorValue2 / (scale * scale);
		  memcpy(&keypointMax, &keypoint2, sizeof(CamKeypointShort));
		}
	    }
	  ++scale;
	}
      memcpy(&features->keypoint[i]->x, &keypointMax, sizeof(CamKeypointShort));
      printf("finished at scale : %i\n", scale);
    }
}

int	camSortKeypointsShort(const void *p1x, const void *p2x);

void				cam_keypoints_tracking2_select_good_features(CamTrackingContext *tc, CamImage *image)
{
  float				val;
  int				x;
  int				y;
  CamKeypointShort		*pointsList;
  CamKeypointShort		*sortedPointsList;
  int				window_hh;
  int				window_hw;
  int				nbPoints;
  int				nrows;
  int				ncols;
  int				nbFound;
  int				borderX;
  int				borderY;
  register CamKeypointShort	*ptr;
  register int			xx;
  register int			yy;
  register float		gx;
  register float		gy;
  register float		gxx;
  register float		gxy;
  register float		gyy;
  register int			i;
  register int			j;

  borderX = tc->b.borderX / tc->pyramidImages->levels[tc->pyramidImages->nbLevels - 1].scale;
  borderY = tc->b.borderY / tc->pyramidImages->levels[tc->pyramidImages->nbLevels - 1].scale;
  window_hh = tc->rw.height / 2;
  window_hw = tc->rw.width / 2;
  ncols = tc->pyramidImages->levels[tc->pyramidImages->nbLevels - 1].img1->image.ncols;
  nrows = tc->pyramidImages->levels[tc->pyramidImages->nbLevels - 1].img1->image.nrows;
  pointsList = (CamKeypointShort *)malloc((nrows - 2 * borderY) * (ncols - 2 * borderX) * sizeof(CamKeypointShort));
  sortedPointsList = (CamKeypointShort *)malloc((nrows - 2 * borderY) * (ncols - 2 * borderX) * sizeof(CamKeypointShort));
  ptr = pointsList;
  nbPoints = 0;
  for (y = borderY ; y < nrows - borderY ; ++y)
    {
      for (x = borderX ; x < ncols - borderX ; ++x)
	{
	  gxx = 0;
	  gxy = 0;
	  gyy = 0;
	  for (yy = y - window_hh ; yy <= y + window_hh ; ++yy)
	    {
	      for (xx = x - window_hw ; xx <= x + window_hw ; ++xx)
		{
		  gx = *(tc->pyramidImages->levels[tc->pyramidImages->nbLevels - 1].img1->gradX.data + ncols * yy + xx);
		  gy = *(tc->pyramidImages->levels[tc->pyramidImages->nbLevels - 1].img1->gradY.data + ncols * yy + xx);
		  gxx += gx * gx;
		  gxy += gx * gy;
		  gyy += gy * gy;
		}
	    }
	  ptr->x = x;
	  ptr->y = y;
	  val = cam_keypoints_tracking2_min_eigen_value(gxx, gxy, gyy);
	  if (val > (float)INT_MAX)
	    val = (float)INT_MAX;
	  ptr->value = (int)val;
	  ptr->scale = 4;
	  ptr->angle = cam_keypoints_tracking2_compute_main_eigen_vector(gxx, gxy, gyy) * 100;
	  ++ptr;
	  nbPoints++;
	}
    }
  memcpy(sortedPointsList, pointsList, nbPoints * sizeof(CamKeypointShort));
  qsort(sortedPointsList, nbPoints, sizeof(CamKeypointShort), camSortKeypointsShort);
  nbFound = cam_keypoints_tracking2_position_filter(tc, pointsList, sortedPointsList, nbPoints, 29);
  for (i = 0 , j = 0 ; j < min(nbFound, tc->nbFeatures) ; ++i)
    {
      if (sortedPointsList[i].value > 1)
	{
	  tc->previousCorners->keypoint[j] = &tc->previousCorners->bag[j];
	  tc->previousFeatures->keypoint[j] = &tc->previousFeatures->bag[j];
	  sortedPointsList[i].x *= tc->pyramidImages->levels[tc->pyramidImages->nbLevels - 1].scale;
	  sortedPointsList[i].y *= tc->pyramidImages->levels[tc->pyramidImages->nbLevels - 1].scale;
	  memcpy(&tc->previousCorners->keypoint[j]->x, &sortedPointsList[i], sizeof(CamKeypointShort));
	  ++j;
	}
    }
  tc->nbDetectedFeatures = j;
  cam_keypoints_tracking2_locate_keypoints(tc, tc->previousCorners, tc->previousFeatures);
  free(sortedPointsList);
  free(pointsList);
}

CamKeypointsMatches	*cam_keypoints_tracking2(CamTrackingContext *tc, CamImage *image, int options)
{
  register int		i;
  CamKeypointsMatches	*res;
  CamROI		roix;

  if (!tc->previousImage)
    {
      for (i = 0 ; i < tc->pyramidImages->nbLevels ; ++i)
	{
	  cam_keypoints_tracking2_copy_image_to_float_image(&tc->pyramidImages->levels[i].img1->image, image, tc->pyramidImages->levels[i].scale);
	  cam_keypoints_tracking2_compute_gradients(&tc->pyramidImages->levels[i].img1->image, &tc->pyramidImages->levels[i].img1->gradX, &tc->pyramidImages->levels[i].img1->gradY, &tc->gaussianKernel, &tc->gaussianDerivedKernel);
	}

      tc->previousImage = image;
      tc->previousIntegralImage = (CamImage*)malloc(sizeof(CamImage));
      tc->previousIntegralImage->imageData = NULL;
      camSetMaxROI(&roix, image);
      roix.coi = 1;
      image->roi = &roix;
      camIntegralImage(image, tc->previousIntegralImage);
      tc->previousIntegralImage->roi = &roix;

      cam_keypoints_tracking2_select_good_features(tc, image);
      return (NULL);
    }
  else
    {
      for (i = 0 ; i < tc->pyramidImages->nbLevels ; ++i)
	{
	  cam_keypoints_tracking2_copy_image_to_float_image(&tc->pyramidImages->levels[i].img2->image, image, tc->pyramidImages->levels[i].scale);
	  cam_keypoints_tracking2_compute_gradients(&tc->pyramidImages->levels[i].img2->image, &tc->pyramidImages->levels[i].img2->gradX, &tc->pyramidImages->levels[i].img2->gradY, &tc->gaussianKernel, &tc->gaussianDerivedKernel);
	}
      res = (CamKeypointsMatches*)malloc(sizeof(CamKeypointsMatches));
      camAllocateKeypointsMatches(res, tc->nbFeatures);      
      for (i = 0 ; i < tc->nbDetectedFeatures ; ++i)
	{
	  res->pairs[i].p1 = (CamKeypoint*)malloc(sizeof(CamKeypoint));
	  res->pairs[i].p2 = (CamKeypoint*)malloc(sizeof(CamKeypoint));
	  res->pairs[i].mark = 1;
	  memcpy(res->pairs[i].p1, tc->previousCorners->keypoint[i], sizeof(CamKeypoint));
	  memcpy(res->pairs[i].p2, tc->previousCorners->keypoint[i], sizeof(CamKeypoint));
	}
      res->nbMatches = tc->nbDetectedFeatures;
      return (res);
    }
}
void		cam_keypoints_tracking2_print_matches(CamImage *img1, CamImage *img2, char *outfile, CamKeypointsMatches *matches)
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
      if (!matches->pairs[i].mark)
	continue ;
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
      if (!matches->pairs[i].mark)
	continue ;
      matches->pairs[i].p2->y += img1->height;
      camDrawKeypoint(matches->pairs[i].p2, &res, 128);
    }
  sprintf(filename, "output/%s.bmp", outfile);
  camSaveBMP(&res, filename);
  camDeallocateImage(&res);
}

void		cam_keypoints_tracking2_release_matches(CamKeypointsMatches *track)
{
  register int	i;

  for (i = 0 ; i < track->nbMatches ; ++i)
    {
      free(track->pairs[i].p1);
      free(track->pairs[i].p2);
      track->pairs[i].mark = 1;
    }
}

void			test_cam_keypoints_tracking2()
{
  CamImage		modelImage;
  CamImage		firstImage;
  CamImage		secondImage;
  CamList		*scales;
  CamTrackingContext	tc;
  CamKeypointsMatches	*track;
  char			img1[] = "./resources/klt/img0.bmp";
  char			img2[] = "./resources/klt/img0.bmp";
#ifdef CAM_TRACKING2_TIMINGS
  int			t1;
  int			t2;
  int			t3;
#endif

  modelImage.imageData = NULL;
  camLoadBMP(&modelImage, img1);
  camAllocateYUVImage(&firstImage, modelImage.width, modelImage.height);
  camRGB2YUV(&modelImage, &firstImage);
  camDeallocateImage(&modelImage);
  camLoadBMP(&modelImage, img2);
  camAllocateYUVImage(&secondImage, modelImage.width, modelImage.height);
  camRGB2YUV(&modelImage, &secondImage);
  camDeallocateImage(&modelImage);

  scales = NULL;
  scales = cam_keypoints_tracking_add_to_linked_list(scales, (void*)1);
  scales = cam_keypoints_tracking_add_to_linked_list(scales, (void*)4);
  
  cam_keypoint_tracking2_configure_context(&tc, 100, 7, 7, 3, 3, 3, scales, &firstImage);
  cam_keypoints_tracking_free_linked_list(scales);

#ifdef CAM_TRACKING2_TIMINGS
  t1 = camGetTimeMs();
#endif
  track = cam_keypoints_tracking2(&tc, &firstImage, CAM_UPRIGHT);
#ifdef CAM_TRACKING2_TIMINGS
  t2 = camGetTimeMs();
  printf("initial detection : %ims\n", t2 - t1);
#endif
  track = cam_keypoints_tracking2(&tc, &secondImage, CAM_UPRIGHT);
#ifdef CAM_TRACKING2_TIMINGS
  t3 = camGetTimeMs();
  printf("tracking : %ims\n", t3 - t2);
#endif
  cam_keypoints_tracking_free_context(&tc);
  camDeallocateImage(&secondImage);

  camLoadBMP(&firstImage, img1);
  camLoadBMP(&secondImage, img2);
  cam_keypoints_tracking2_print_matches(&firstImage, &secondImage, "track", track);
  camDeallocateImage(&firstImage);
  camDeallocateImage(&secondImage);
  cam_keypoints_tracking2_release_matches(track);
  camFreeKeypointsMatches(track);
  free(track);
}
