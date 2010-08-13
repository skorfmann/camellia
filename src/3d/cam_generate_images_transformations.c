
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cam_matrix.h"
#include "cam_2d_points.h"
#include "cam_list.h"
#include "misc.h"
#include "cam_ppm_to_points.h"
#include "cam_points_to_ppm.h"
#include "cam_ppm_to_matrix.h"
#include "cam_matrix_to_ppm.h"
#include "cam_matrix_to_points.h"
#include "cam_points_to_matrix.h"
#include "cam_interpolate_missing_image_data.h"

void		cam_homography_to_file(char *dstFile, CamMatrix *H)
{
  FILE		*file;

  file = fopen(dstFile, "w");
  if (!file)
    {
      printf("cam_homography_to_file : unable to create file %s\n", dstFile);
      exit (-1);
    }
  fwrite(&H->ncols, sizeof(int), 1, file);
  fwrite(&H->nrows, sizeof(int), 1, file);
  fwrite(H->data, sizeof(POINTS_TYPE), H->ncols * H->nrows, file);
  fclose (file);
}

CamList		*cam_apply_transformation_list(CamList *points, CamMatrix *H)
{
  CamList	*res;
  CamMatrix	pt;
  CamMatrix	ptTransformed;
  CamList	*ptr;

  cam_allocate_matrix(&pt, 1, 3);
  cam_allocate_matrix(&ptTransformed, 1, 3);
  res = NULL;
  ptr = points;
  cam_matrix_set_value(&pt, 0, 2, 1.0f);
  while (ptr)
    {
      cam_matrix_set_value(&pt, 0, 0, ((CamColorized2dPoint *)(ptr->data))->point.x);
      cam_matrix_set_value(&pt, 0, 1, ((CamColorized2dPoint *)(ptr->data))->point.y);
      cam_matrix_multiply(&ptTransformed, H, &pt);
      res = cam_add_to_linked_list(res, (CamColorized2dPoint *)malloc(sizeof(CamColorized2dPoint)));
      ((CamColorized2dPoint *)(res->data))->point.x = cam_matrix_get_value(&ptTransformed, 0, 0);
      ((CamColorized2dPoint *)(res->data))->point.y = cam_matrix_get_value(&ptTransformed, 0, 1);
      ((CamColorized2dPoint *)(res->data))->color.r = ((CamColorized2dPoint *)(ptr->data))->color.r;
      ((CamColorized2dPoint *)(res->data))->color.g = ((CamColorized2dPoint *)(ptr->data))->color.g;
      ((CamColorized2dPoint *)(res->data))->color.b = ((CamColorized2dPoint *)(ptr->data))->color.b;
      ptr = ptr->next;
    }
  cam_disallocate_matrix(&pt);
  cam_disallocate_matrix(&ptTransformed);
  return (res);
}

CamImageMatrix		*cam_apply_transformation_imagematrix(CamImageMatrix *img, CamMatrix *H, unsigned char bgR, unsigned char bgG, unsigned char bgB)
{
  CamImageMatrix	*res;
  CamMatrix		pt;
  CamMatrix		ptTransformed;
  int			i;
  int			j;

  cam_allocate_matrix(&pt, 1, 3);
  cam_allocate_matrix(&ptTransformed, 1, 3);
  res = (CamImageMatrix *)malloc(sizeof(CamImageMatrix));
  cam_allocate_imagematrix(res, img->r.ncols, img->r.nrows);
  cam_matrix_set_value(&pt, 0, 2, 1.0f);
  for (j = 0 ; j < img->r.nrows ; ++j)
    {
      for (i = 0 ; i < img->r.ncols ; ++i)
	{
	  cam_matrix_set_value(&res->r, i, j, bgR);
	  cam_matrix_set_value(&res->g, i, j, bgG);
	  cam_matrix_set_value(&res->b, i, j, bgB);
	}
    }
  for (j = 0 ; j < img->r.nrows ; ++j)
    {
      for (i = 0 ; i < img->r.ncols ; ++i)
	{
	  cam_matrix_set_value(&pt, 0, 0, (POINTS_TYPE)i - img->r.ncols / 2);
	  cam_matrix_set_value(&pt, 0, 1, (POINTS_TYPE)j - img->r.nrows / 2);
	  cam_matrix_multiply(&ptTransformed, H, &pt);
	  cam_matrix_add_value(&ptTransformed, 0, 0, (POINTS_TYPE)img->r.ncols / 2);
	  cam_matrix_add_value(&ptTransformed, 0, 1, (POINTS_TYPE)img->r.nrows / 2);
	  if ((int)cam_matrix_get_value(&ptTransformed, 0, 0) >= 0 && (int)cam_matrix_get_value(&ptTransformed, 0, 0) < img->r.ncols &&
	      (int)cam_matrix_get_value(&ptTransformed, 0, 1) >= 0 && (int)cam_matrix_get_value(&ptTransformed, 0, 1) < img->r.nrows)
	    {
	      cam_matrix_set_value(&res->r, cam_matrix_get_value(&ptTransformed, 0, 0), cam_matrix_get_value(&ptTransformed, 0, 1),
				   cam_matrix_get_value(&img->r, i, j));
	      cam_matrix_set_value(&res->g, cam_matrix_get_value(&ptTransformed, 0, 0), cam_matrix_get_value(&ptTransformed, 0, 1),
				   cam_matrix_get_value(&img->g, i, j));
	      cam_matrix_set_value(&res->b, cam_matrix_get_value(&ptTransformed, 0, 0), cam_matrix_get_value(&ptTransformed, 0, 1),
				   cam_matrix_get_value(&img->b, i, j));
	    }
	}
    }
  cam_disallocate_matrix(&pt);
  cam_disallocate_matrix(&ptTransformed);
  return (res);
}

CamList		*cam_affine_transform_list(CamList *points, char *dir, char *file, POINTS_TYPE lambda1, POINTS_TYPE lambda2, POINTS_TYPE phi, POINTS_TYPE theta, POINTS_TYPE tx, POINTS_TYPE ty)
{
  CamMatrix	H;
  CamMatrix	D;
  CamMatrix	Rtheta;
  CamMatrix	Rphi;
  CamMatrix	Rminusphi;
  CamMatrix	tmp1;
  CamMatrix	tmp2;
  CamList	*res;
  char		*outputPath;

  outputPath = (char *)malloc((strlen(dir) + strlen(file) + strlen("tr") + 3) * sizeof(char) );
  sprintf(outputPath,"%s/%s.%s", dir, file, "tr");
  cam_allocate_matrix(&H, 3, 3);
  cam_allocate_matrix(&D, 2, 2);
  cam_allocate_matrix(&Rtheta, 2, 2);
  cam_allocate_matrix(&Rphi, 2, 2);
  cam_allocate_matrix(&Rminusphi, 2, 2);
  cam_allocate_matrix(&tmp1, 2, 2);
  cam_allocate_matrix(&tmp2, 2, 2);

  cam_matrix_set_value(&D, 0, 0, lambda1);
  cam_matrix_set_value(&D, 0, 1, 0.0f);
  cam_matrix_set_value(&D, 1, 0, 0.0f);
  cam_matrix_set_value(&D, 1, 1, lambda2);

  cam_matrix_set_value(&Rtheta, 0, 0, cos(theta));
  cam_matrix_set_value(&Rtheta, 0, 1, -sin(theta));
  cam_matrix_set_value(&Rtheta, 1, 0, sin(theta));
  cam_matrix_set_value(&Rtheta, 1, 1, cos(theta));

  cam_matrix_set_value(&Rphi, 0, 0, cos(phi));
  cam_matrix_set_value(&Rphi, 0, 1, -sin(phi));
  cam_matrix_set_value(&Rphi, 1, 0, sin(phi));
  cam_matrix_set_value(&Rphi, 1, 1, cos(phi));
  
  cam_matrix_set_value(&Rminusphi, 0, 0, cos(-phi));
  cam_matrix_set_value(&Rminusphi, 0, 1, -sin(-phi));
  cam_matrix_set_value(&Rminusphi, 1, 0, sin(-phi));
  cam_matrix_set_value(&Rminusphi, 1, 1, cos(-phi));

  cam_matrix_multiply(&tmp1, &Rtheta, &Rminusphi);
  cam_matrix_multiply(&tmp2, &tmp1, &D);
  cam_matrix_multiply(&tmp1, &tmp2, &Rphi);

  cam_matrix_set_value(&H, 0, 0, cam_matrix_get_value(&tmp1, 0, 0));
  cam_matrix_set_value(&H, 0, 1, cam_matrix_get_value(&tmp1, 0, 1));
  cam_matrix_set_value(&H, 0, 2, 0.0f);
  cam_matrix_set_value(&H, 1, 0, cam_matrix_get_value(&tmp1, 1, 0));
  cam_matrix_set_value(&H, 1, 1, cam_matrix_get_value(&tmp1, 1, 1));
  cam_matrix_set_value(&H, 1, 2, 0.0f);
  cam_matrix_set_value(&H, 2, 0, tx);
  cam_matrix_set_value(&H, 2, 1, ty);
  cam_matrix_set_value(&H, 2, 2, 1.0f);

  res = cam_apply_transformation_list(points, &H);
  cam_homography_to_file(outputPath, &H);

  cam_disallocate_matrix(&H);
  cam_disallocate_matrix(&D);
  cam_disallocate_matrix(&Rtheta);
  cam_disallocate_matrix(&Rphi);
  cam_disallocate_matrix(&Rminusphi);
  cam_disallocate_matrix(&tmp1);
  cam_disallocate_matrix(&tmp2);
  free(outputPath);
  return (res);
}

CamList		*cam_similarity_transform_list(CamList *points, char *dir, char *file, POINTS_TYPE s, POINTS_TYPE theta, POINTS_TYPE tx, POINTS_TYPE ty)
{
  CamMatrix	H;
  CamList	*res;
  char		*outputPath;

  outputPath = (char *)malloc((strlen(dir) + strlen(file) + strlen("tr") + 3) * sizeof(char) );
  sprintf(outputPath,"%s/%s.%s", dir, file, "tr");
  cam_allocate_matrix(&H, 3, 3);
  cam_matrix_set_value(&H, 0, 0, s * cos(theta));
  cam_matrix_set_value(&H, 0, 1, s * -sin(theta));
  cam_matrix_set_value(&H, 0, 2, 0.0f);
  cam_matrix_set_value(&H, 1, 0, s * sin(theta));
  cam_matrix_set_value(&H, 1, 1, s * cos(theta));
  cam_matrix_set_value(&H, 1, 2, 0.0f);
  cam_matrix_set_value(&H, 2, 0, tx);
  cam_matrix_set_value(&H, 2, 1, ty);
  cam_matrix_set_value(&H, 2, 2, 1.0f);
  res = cam_apply_transformation_list(points, &H);
  cam_homography_to_file(outputPath, &H);
  cam_disallocate_matrix(&H);
  free(outputPath);
  return (res);
}

CamImageMatrix		*cam_similarity_transform_imagematrix(CamImageMatrix *img, char *dir, char *file, POINTS_TYPE s, POINTS_TYPE theta, POINTS_TYPE tx, POINTS_TYPE ty, unsigned char bgR, unsigned char bgG, unsigned char bgB)
{
  CamMatrix		H;
  CamImageMatrix	*res;
  char			*outputPath;

  outputPath = (char *)malloc((strlen(dir) + strlen(file) + strlen("tr") + 3) * sizeof(char) );
  sprintf(outputPath,"%s/%s.%s", dir, file, "tr");
  cam_allocate_matrix(&H, 3, 3);
  cam_matrix_set_value(&H, 0, 0, s * cos(theta));
  cam_matrix_set_value(&H, 0, 1, s * -sin(theta));
  cam_matrix_set_value(&H, 0, 2, 0.0f);
  cam_matrix_set_value(&H, 1, 0, s * sin(theta));
  cam_matrix_set_value(&H, 1, 1, s * cos(theta));
  cam_matrix_set_value(&H, 1, 2, 0.0f);
  cam_matrix_set_value(&H, 2, 0, tx);
  cam_matrix_set_value(&H, 2, 1, ty);
  cam_matrix_set_value(&H, 2, 2, 1.0f);
  res = cam_apply_transformation_imagematrix(img, &H, bgR, bgG, bgB);
  cam_homography_to_file(outputPath, &H);
  cam_disallocate_matrix(&H);
  free(outputPath);
  return (res);
}

CamList		*cam_euclidian_transform_list(CamList *points, char *dir, char *file, POINTS_TYPE theta, POINTS_TYPE tx, POINTS_TYPE ty)
{
  CamMatrix	H;
  CamList	*res;
  char		*outputPath;

  outputPath = (char *)malloc((strlen(dir) + strlen(file) + strlen("tr") + 3) * sizeof(char) );
  sprintf(outputPath, "%s/%s.%s", dir, file, "tr");
  cam_allocate_matrix(&H, 3, 3);
  cam_matrix_set_value(&H, 0, 0, cos(theta));
  cam_matrix_set_value(&H, 0, 1, -sin(theta));
  cam_matrix_set_value(&H, 0, 2, 0.0f);
  cam_matrix_set_value(&H, 1, 0, sin(theta));
  cam_matrix_set_value(&H, 1, 1, cos(theta));
  cam_matrix_set_value(&H, 1, 2, 0.0f);
  cam_matrix_set_value(&H, 2, 0, tx);
  cam_matrix_set_value(&H, 2, 1, ty);
  cam_matrix_set_value(&H, 2, 2, 1.0f);
  res = cam_apply_transformation_list(points, &H);
  cam_homography_to_file(outputPath, &H);
  cam_disallocate_matrix(&H);
  free(outputPath);

  return (res);
}

int			main()
{
  char			outputDir[] = "data/tracking";
  CamImageMatrix	*image;
  CamImageMatrix	*imageTransformed2;
  CamImageMatrix	*imageInterpoled;
  CamList		*pts;
  CamMatrix		inverseHomography;
  POINTS_TYPE		angle1;
  POINTS_TYPE		tx1;
  POINTS_TYPE		ty1;
  POINTS_TYPE		ratio;
  unsigned char		bgR, bgG, bgB;
  CamImageMatrix	tst;
  
 
  int	i;
  int	j;

  angle1 = PI/70;
  ratio = 1.0f;
  tx1 = 0.0f;
  ty1 = 0.0f;

  bgR = 0;
  bgG = 0;
  bgB = 255;

  cam_allocate_imagematrix(&tst, 255, 255);
  for (j = 0 ; j < 255 ; ++j)
    {
      for (i = 0 ; i < 255 ; ++i)
	{
	  if (i > 128)
	    cam_matrix_set_value(&tst.r, i , j, (POINTS_TYPE)0);
	  else
	    cam_matrix_set_value(&tst.r, i , j, (POINTS_TYPE)255);
	  if (j > 128)
	    cam_matrix_set_value(&tst.g, i , j, (POINTS_TYPE)255);
	  else
	    cam_matrix_set_value(&tst.g, i , j, (POINTS_TYPE)0);
	  cam_matrix_set_value(&tst.b, i , j, 0.0f);
	}
    }
  cam_matrix_to_pathppm(outputDir, "testimg", &tst);

  image = cam_ppm_to_matrix("data/tracking/img0.ppm");
  cam_allocate_matrix(&inverseHomography, 3, 3);
  pts = cam_matrix_to_points(image);
  imageTransformed2 = cam_similarity_transform_imagematrix(image, outputDir, "homo1", ratio, angle1, tx1, ty1, bgR, bgG, bgB);
  cam_matrix_set_value(&inverseHomography, 0, 0, cos(-angle1) / ratio);
  cam_matrix_set_value(&inverseHomography, 0, 1, -sin(-angle1) / ratio);
  cam_matrix_set_value(&inverseHomography, 0, 2, 0.0f);
  cam_matrix_set_value(&inverseHomography, 1, 0, sin(-angle1) / ratio);
  cam_matrix_set_value(&inverseHomography, 1, 1, cos(-angle1) / ratio);
  cam_matrix_set_value(&inverseHomography, 1, 2, 0.0f);
  cam_matrix_set_value(&inverseHomography, 2, 0, -tx1);
  cam_matrix_set_value(&inverseHomography, 2, 1, -ty1);
  cam_matrix_set_value(&inverseHomography, 2, 2, 1.0f);
  imageInterpoled = cam_interpolate_missing_image_data_bilinear(imageTransformed2, image, &inverseHomography, bgR, bgG, bgB);
  cam_matrix_to_pathppm(outputDir, "transformed", imageTransformed2);
  cam_matrix_to_pathppm(outputDir, "img1", imageInterpoled);
  
  cam_disallocate_linked_list(pts);
  cam_disallocate_imagematrix(image);
  cam_disallocate_imagematrix(imageTransformed2);
  cam_disallocate_imagematrix(imageInterpoled);
  cam_disallocate_matrix(&inverseHomography);
  free(image);
  free(imageTransformed2);
  free(imageInterpoled);
  return (0);
}
