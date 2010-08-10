
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


CamList		*cam_apply_transformation(CamList *points, CamMatrix *H)
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

CamList		*cam_affine_transform(CamList *points, char *dir, char *file, POINTS_TYPE lambda1, POINTS_TYPE lambda2, POINTS_TYPE phi, POINTS_TYPE theta, POINTS_TYPE tx, POINTS_TYPE ty)
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

  res = cam_apply_transformation(points, &H);
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

CamList		*cam_similarity_transform(CamList *points, char *dir, char *file, POINTS_TYPE s, POINTS_TYPE theta, POINTS_TYPE tx, POINTS_TYPE ty)
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
  res = cam_apply_transformation(points, &H);
  cam_homography_to_file(outputPath, &H);
  cam_disallocate_matrix(&H);
  free(outputPath);
  return (res);
}

CamList		*cam_euclidian_transform(CamList *points, char *dir, char *file, POINTS_TYPE theta, POINTS_TYPE tx, POINTS_TYPE ty)
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
  res = cam_apply_transformation(points, &H);
  cam_homography_to_file(outputPath, &H);
  cam_disallocate_matrix(&H);
  free(outputPath);

  return (res);
}

int			main()
{
  char			outputDir[] = "data/tracking";
  POINTS_TYPE		mask_data[] = {1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f};
  CamMatrix		mask;
  CamImageMatrix	*image;
  CamImageMatrix	imageTransformed;
  CamImageMatrix	*imageTransformed2;
  CamImageMatrix	*imageInterpoled;
  CamList		*pts;
  CamList		*ptsTransformed;
  
  image = cam_ppm_to_matrix("data/tracking/img0.ppm");
  cam_allocate_imagematrix(&imageTransformed, image->r.ncols, image->r.nrows);
  cam_allocate_matrix(&mask, 3, 3);
  memcpy(mask.data,mask_data,9 * sizeof(POINTS_TYPE));

  cam_matrix_image_convolution(&imageTransformed, image, &mask, 9.0f);
  /*  cam_matrix_to_pathppm(outputDir, "image", image);
      cam_matrix_to_pathppm(outputDir, "imageTransformed", &imageTransformed);*/

  /*  pts = cam_matrix_to_points(&imageTransformed);*/
  pts = cam_matrix_to_points(image);
  ptsTransformed = cam_similarity_transform(pts, outputDir, "homo1", 1.0f, PI/4, 0.0f, 0.0f);
  imageTransformed2 = cam_points_to_matrix(ptsTransformed, image->r.ncols, image->r.nrows);
  /*cam_matrix_to_pathppm(outputDir, "img", imageTransformed2);*/
  imageInterpoled = cam_interpolate_missing_image_data_median(imageTransformed2, 1, 0, 0, 0);
  cam_matrix_to_pathppm(outputDir, "img1", imageInterpoled);
  
  cam_disallocate_linked_list(pts);
  cam_disallocate_linked_list(ptsTransformed);
  cam_disallocate_imagematrix(image);
  cam_disallocate_imagematrix(&imageTransformed);
  cam_disallocate_imagematrix(imageTransformed2);
  cam_disallocate_imagematrix(imageInterpoled);
  cam_disallocate_matrix(&mask);
  free(image);
  free(imageTransformed2);
  free(imageInterpoled);
  return (0);
}
