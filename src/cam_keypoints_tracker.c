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
#define MAX_KERNEL_WIDTH 71

//#define CAM_TRACKING_SUBTIMINGS
#define CAM_TRACKING_TIMINGS
//#define	CAM_TRACKING_DEBUG_4
//#define	CAM_TRACKING_DEBUG_3
//#define	CAM_TRACKING_DEBUG_2
//#define CAM_TRACKING_DEBUG_1

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

typedef struct
{
  int	width;
  float	data[MAX_KERNEL_WIDTH];
} CamConvolutionKernel;

typedef struct
{
  int	ncols;
  int	nrows;
  float	*data;
} CamFloatImage;

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
  free (res->data);
}

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
		      abs(tc->previousFeatures->keypoint[(i + k) * width + (j + l)]->value) > valueMax &&
		      (i + k) && (j + l) && (i + k) < width - 1 && (j + l) < width - 1 )
		    {
		      valueMax = abs(tc->previousFeatures->keypoint[(i + k) * width + (j + l)]->value);
		      seedsIndexes[m] = (i + k) * width + (j + l);
		    }
		}
	    }
#ifdef CAM_TRACKING_DEBUG_3
	  if (tc->previousFeatures->keypoint[i * width + j])
	    printf("Seed %i with value %i is better in %i with value %i\n", i * width + j, tc->previousFeatures->keypoint[i * width + j]->value, seedsIndexes[m], tc->previousFeatures->keypoint[seedsIndexes[m]]->value);
	  else
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

void		cam_keypoints_tracking_compute_kernels(float sigma, CamConvolutionKernel *gauss, CamConvolutionKernel *gaussderiv)
{
  const float	factor = 0.01f;
  int i;

  {
    const int hw = MAX_KERNEL_WIDTH / 2;
    float max_gauss = 1.0f, max_gaussderiv = (float) (sigma*exp(-0.5f));
    
    for (i = -hw ; i <= hw ; i++)  {
      gauss->data[i+hw]      = (float) exp(-i*i / (2*sigma*sigma));
      gaussderiv->data[i+hw] = -i * gauss->data[i+hw];
    }

    gauss->width = MAX_KERNEL_WIDTH;
    for (i = -hw ; fabs(gauss->data[i+hw] / max_gauss) < factor ; 
         i++, gauss->width -= 2);
    gaussderiv->width = MAX_KERNEL_WIDTH;
    for (i = -hw ; fabs(gaussderiv->data[i+hw] / max_gaussderiv) < factor ; 
         i++, gaussderiv->width -= 2);
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

void	cam_keypoints_tracking_convolve_horiz(CamFloatImage *imgin, CamConvolutionKernel *kernel, CamFloatImage *imgout)
{
  float *ptrrow = imgin->data;
  register float *ptrout = imgout->data,
    *ppp;
  register float sum;
  register int radius = kernel->width / 2;
  register int ncols = imgin->ncols, nrows = imgin->nrows;
  register int i, j, k;
  
  for (j = 0 ; j < nrows ; j++)
    {
      
      for (i = 0 ; i < radius ; i++)
	*ptrout++ = 0.0;
      
      for ( ; i < ncols - radius ; i++)
	{
	  ppp = ptrrow + i - radius;
	  sum = 0.0;
	  //for (k = 0 ; k < kernel->width ; k++)
	  for (k = kernel->width-1 ; k >= 0 ; k--)
	    sum += *ppp++ * kernel->data[k];
	  *ptrout++ = sum;
	}

      for ( ; i < ncols ; i++)
	*ptrout++ = 0.0;
      
      ptrrow += ncols;
    }
}

void	cam_keypoints_tracking_convolve_vert(CamFloatImage *imgin, CamConvolutionKernel *kernel, CamFloatImage *imgout)
{
  float *ptrcol = imgin->data;
  register float *ptrout = imgout->data,
    *ppp;
  register float sum;
  register int radius = kernel->width / 2;
  register int ncols = imgin->ncols, nrows = imgin->nrows;
  register int i, j, k;

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
	  //for (k = 0 ; k < kernel->width-1 ; k++)
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

void	cam_keypoints_tracking_convolve_separate(CamFloatImage *imgin, CamConvolutionKernel *horiz_kernel, CamConvolutionKernel *vert_kernel, CamFloatImage *imgout)
{
  CamFloatImage tmpimg;

  CamAllocateFloatImage(&tmpimg, imgin->ncols, imgin->nrows);
  cam_keypoints_tracking_convolve_horiz(imgin, horiz_kernel, &tmpimg);
  cam_keypoints_tracking_convolve_vert(&tmpimg, vert_kernel, imgout);
  CamDisallocateFloatImage(&tmpimg);
}

void			cam_keypoints_tracking_copy_image_to_float_image(CamFloatImage *dst, CamImage *src)
{
  register unsigned int	i;

  CamAllocateFloatImage(dst, src->width, src->height);
  for (i = 0 ; i < src->width * src->height ; ++i)
    dst->data[i] = (float)src->imageData[i];
}

void		cam_keypoints_tracking_compute_gradients(CamFloatImage *img, CamFloatImage *gradx, CamFloatImage *grady, CamConvolutionKernel *gauss_kernel, CamConvolutionKernel *gaussderiv_kernel)
{
  cam_keypoints_tracking_convolve_separate(img, gaussderiv_kernel, gauss_kernel, gradx);
  cam_keypoints_tracking_convolve_separate(img, gauss_kernel, gaussderiv_kernel, grady);
}

float cam_keypoints_tracking_interpolate(float x, float y, CamFloatImage *img)
{
  int xt;
  int yt;
  float ax;
  float ay;
  float *ptr;

  xt = (int) x;
  yt = (int) y;
  ax = x - xt;
  ay = y - yt;
  ptr = img->data + (img->ncols*yt) + xt;
  return ( (1-ax) * (1-ay) * *ptr +
           ax   * (1-ay) * *(ptr+1) +
           (1-ax) *   ay   * *(ptr+(img->ncols)) +
           ax   *   ay   * *(ptr+(img->ncols)+1) );
}

void			cam_keypoints_tracking_compute_intensity_difference(CamFloatImage *img1, CamFloatImage *img2, float x1, float y1, float x2, float y2, int width, int height, CamFloatImage *imgdiff)
{
  register int		hw;
  register int		hh;
  float			g1;
  float			g2;
  register int		i;
  register int		j;
  register float	*ptrimgdiff;
  
  hw = width / 2;
  hh = height / 2;
  ptrimgdiff = imgdiff->data;
  for (j = -hh ; j <= hh ; j++)
    for (i = -hw ; i <= hw ; i++)
      {
	g1 = cam_keypoints_tracking_interpolate(x1+i, y1+j, img1);
	g2 = cam_keypoints_tracking_interpolate(x2+i, y2+j, img2);
	*ptrimgdiff++ = g1 - g2;
      }
}

void		cam_keypoints_tracking_compute_gradient_sum(CamFloatImage *gradx1, CamFloatImage *grady1, CamFloatImage *gradx2, CamFloatImage *grady2, float x1, float y1, float x2, float y2, int width, int height, CamFloatImage * gradx, CamFloatImage *grady)
{
  register int	hw;
  register int	hh;
  float		g1;
  float		g2;
  register int	i;
  register int	j;
  float		*ptrgradx;
  float		*ptrgrady;

  hw = width/2;
  hh = height/2;
  ptrgradx = gradx->data;
  ptrgrady = grady->data;
  for (j = -hh ; j <= hh ; j++)
    {
      for (i = -hw ; i <= hw ; i++)
	{
	  g1 = cam_keypoints_tracking_interpolate(x1+i, y1+j, gradx1);
	  g2 = cam_keypoints_tracking_interpolate(x2+i, y2+j, gradx2);
	  *ptrgradx++ = g1 + g2;
	  g1 = cam_keypoints_tracking_interpolate(x1+i, y1+j, grady1);
	  g2 = cam_keypoints_tracking_interpolate(x2+i, y2+j, grady2);
	  *ptrgrady++ = g1 + g2;
	}
    }
}

void			cam_keypoints_tracking_compute_2by2_gradient_matrix(CamFloatImage *gradx, CamFloatImage *grady, int width, int height,
							 float *gxx, float *gxy, float *gyy) 

{
  register float	gx;
  register float	gy;
  register		int i;
  float			*ptrgradx;
  float			*ptrgrady;

  ptrgradx = gradx->data;
  ptrgrady = grady->data;
  *gxx = 0.0;
  *gxy = 0.0;
  *gyy = 0.0;
  for (i = 0 ; i < width * height ; i++)
    {
      gx = *ptrgradx++;
      gy = *ptrgrady++;
      *gxx += gx*gx;
      *gxy += gx*gy;
      *gyy += gy*gy;
    }
}

void			cam_keypoints_tracking_compute_2by1_error_vector(CamFloatImage *imgdiff, CamFloatImage *gradx, CamFloatImage *grady, int width, int height, float step_factor, float *ex, float *ey)
{
  register float	diff;
  register int		i;
  float			*ptrgradx;
  float			*ptrgrady;
  float			*ptrimgdiff;

  ptrgradx = gradx->data;
  ptrgrady = grady->data;
  ptrimgdiff = imgdiff->data;
  /* Compute values */
  *ex = 0;  *ey = 0;  
  for (i = 0 ; i < width * height ; i++)  {
    diff = *ptrimgdiff++;
    *ex += diff * (*ptrgradx++);
    *ey += diff * (*ptrgrady++);
  }
  *ex *= step_factor;
  *ey *= step_factor;
}

void	cam_keypoints_tracking_solve_mouvement_equation(float gxx, float gxy, float gyy, float ex, float ey,
					      float small, float *dx, float *dy)
{
  float det = gxx*gyy - gxy*gxy;

  
  if (det < small)  return;

  *dx = (gyy*ex - gxy*ey)/det;
  *dy = (gxx*ey - gxy*ex)/det;
}

void		cam_keypoints_tracking_compute_local_image_displacement(float x1, float y1, float *x2, float *y2, CamFloatImage *img1, CamFloatImage *gradx1, CamFloatImage *grady1, CamFloatImage *img2, CamFloatImage *gradx2, CamFloatImage *grady2, int width, int height, float step_factor, float small)
{
  CamFloatImage	imgdiff;
  CamFloatImage	gradx;
  CamFloatImage	grady;
  float		gxx, gxy, gyy, ex, ey, dx, dy;
  int		iteration = 0;
  int		hw = width/2;
  int		hh = height/2;
  int		nc = img1->ncols;
  int		nr = img1->nrows;
  float		one_plus_eps = 1.001f;
#ifdef CAM_TRACKING_SUBTIMINGS
  int			deltaTimers1;
  struct timeval	tv1;
  struct timeval	tv2;
  int			deltaTimers2;
  struct timeval	tv3;
  struct timeval	tv4;
#endif

  CamAllocateFloatImage(&imgdiff, width, height);
  CamAllocateFloatImage(&gradx, width, height);
  CamAllocateFloatImage(&grady, width, height);
    
  do  {

    if (  x1-hw < 0.0f || nc-( x1+hw) < one_plus_eps ||
	  *x2-hw < 0.0f || nc-(*x2+hw) < one_plus_eps ||
          y1-hh < 0.0f || nr-( y1+hh) < one_plus_eps ||
	  *y2-hh < 0.0f || nr-(*y2+hh) < one_plus_eps) {
      break;
    }

#ifdef CAM_TRACKING_SUBTIMINGS
  gettimeofday(&tv1, NULL);      
#endif
  cam_keypoints_tracking_compute_intensity_difference(img1, img2, x1, y1, *x2, *y2, width, height, &imgdiff);
#ifdef CAM_TRACKING_SUBTIMINGS
  gettimeofday(&tv2, NULL);
  deltaTimers1 = tv2.tv_usec - tv1.tv_usec;
  if (deltaTimers1 < 0)
    deltaTimers1 = 1000000 + deltaTimers2;
  printf("Intensity computation : %ius\n", deltaTimers1);
#endif

    cam_keypoints_tracking_compute_gradient_sum(gradx1, grady1, gradx2, grady2, x1, y1, *x2, *y2, width, height, &gradx, &grady);

    cam_keypoints_tracking_compute_2by2_gradient_matrix(&gradx, &grady, width, height, &gxx, &gxy, &gyy);
    cam_keypoints_tracking_compute_2by1_error_vector(&imgdiff, &gradx, &grady, width, height, step_factor, &ex, &ey);
    
    dx = 0.0f;
    dy = 0.0f;
    cam_keypoints_tracking_solve_mouvement_equation(gxx, gxy, gyy, ex, ey, small, &dx, &dy);

    *x2 += dx;
    *y2 += dy;
    iteration++;

  }  while ((fabs(dx)>=0.1f || fabs(dy)>=0.1f) && iteration < 10);

#ifdef CAM_TRACKING_DEBUG_3
  printf("Iterations : %i\n", iteration);
#endif
  
  CamDisallocateFloatImage(&imgdiff);
  CamDisallocateFloatImage(&gradx);
  CamDisallocateFloatImage(&grady);
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
  float			x1, y1, x2, y2;
  int			windowSize;
  CamConvolutionKernel	gaussKernel;
  CamConvolutionKernel	gaussDerivKernel;
  int			kernelSize;
  CamFloatImage		gradX1, gradY1, gradX2, gradY2, img1, img2;
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
  
  CamAllocateFloatImage(&gradX1, image->width, image->height);
  CamAllocateFloatImage(&gradY1, image->width, image->height);
  CamAllocateFloatImage(&gradX2, image->width, image->height);
  CamAllocateFloatImage(&gradY2, image->width, image->height);

  // in the specified research volume, finds the more coherent point corresponding to the current seed <=> highest detector value
  l = 0;
  windowSize = max(min(tc->rw.height, tc->rw.width) / 3, 5);
  if (!(windowSize % 2))
    windowSize++;
  /*
    kernelSize = max(windowSize / 2, 3);
  if (!(kernelSize % 2))
    kernelSize++;
  */
  kernelSize = 5;
  printf("windowSize : %i kernelSize : %i\n", windowSize, kernelSize);
#ifdef CAM_TRACKING_SUBTIMINGS
  gettimeofday(&tv3, NULL);      
#endif
  cam_keypoints_tracking_compute_kernels(1.0f, &gaussKernel, &gaussDerivKernel);
#ifdef CAM_TRACKING_SUBTIMINGS
  gettimeofday(&tv4, NULL);
  deltaTimers2 = tv4.tv_usec - tv3.tv_usec;
  if (deltaTimers2 < 0)
    deltaTimers2 = 1000000 + deltaTimers2;
  printf("Kernels computation : %ius\n", deltaTimers2);
#endif

  cam_keypoints_tracking_copy_image_to_float_image(&img1, tc->previousImage);
  cam_keypoints_tracking_copy_image_to_float_image(&img2, image);
  cam_keypoints_tracking_compute_gradients(&img1, &gradX1, &gradY1, &gaussKernel, &gaussDerivKernel);
  cam_keypoints_tracking_compute_gradients(&img2, &gradY1, &gradY2, &gaussKernel, &gaussDerivKernel);

  for (i = 0 ; i < tc->nbSeeds ; ++i)
    {
      x1 = (float)tc->previousFeatures->keypoint[seedsIndexes[i]]->x;
      y1 = (float)tc->previousFeatures->keypoint[seedsIndexes[i]]->y;
      x2 = x1;
      y2 = y1;
      cam_keypoints_tracking_compute_local_image_displacement(x1, y1, &x2, &y2, &img1, &gradX1, &gradY1, &img2, &gradX2, &gradY2, 7, 7, 1.0f, 0.001);
      maxOnEachScale = cam_keypoints_tracking_extract_overall_max_on_each_scale(tc, integralImage, detectorValue, seedsIndexes[i], (int)(x2 - x1), (int)(y2 - y1), 0);
#ifdef CAM_TRACKING_DEBUG_2
      printf("shiftx : %f shifty : %f\n", x2 - x1, y2 - y1);
#endif
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
  CamDisallocateFloatImage(&img1);
  CamDisallocateFloatImage(&img2);
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
	  printf("Missmatch of %ith: %i %i %i %i, %i %i %i %i\n", i, featuresMatches->pairs[i].p1->x, featuresMatches->pairs[i].p1->y, featuresMatches->pairs[i].p1->scale, featuresMatches->pairs[i].p1->value, featuresMatches->pairs[i].p2->x,featuresMatches->pairs[i].p2->y,featuresMatches->pairs[i].p2->scale, featuresMatches->pairs[i].p2->value);
#endif
	  featuresMatches->pairs[i].mark = 0;
	  featuresMatches->pairs[i].error = currentDistance; // the first better found, maybe others
	  ++featuresMatches->nbOutliers;
	}
    }
#ifdef	CAM_TRACKING_DEBUG_2
  printf("Total inlier nb : %i // Total outlier nb : %i\n", featuresMatches->nbMatches, featuresMatches->nbOutliers);
#endif
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

      for (i = 0 ; i < tc->previousFeatures->nbPoints - tc->nbSeeds ; ++i)
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
  tc->rv.width = 3;
  tc->rv.height = 3;
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
  //char			img1[] = "./resources/rover/translation1.bmp";
  //char			img2[] = "./resources/rover/translation2.bmp";
  char			img1[] = "./resources/rover/rotation1.bmp";
  char			img2[] = "./resources/rover/rotation2.bmp";
  //  char			img1[] = "./resources/klt/img0.bmp";
  //char			img2[] = "./resources/klt/img2.bmp";

  /* begin images building */
  modelImage.imageData = NULL;
  camLoadBMP(&modelImage, img1);
  camAllocateYUVImage(&firstImage, modelImage.width, modelImage.height);
  camRGB2YUV(&modelImage, &firstImage);
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
