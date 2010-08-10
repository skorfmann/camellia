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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cam_matrix.h"
#include "cam_list.h"
#include "cam_evaluate_tracking.h"
#include "misc.h"

/* TODO : utiliser les CamKeypointsMatch ? */

CamMatrix	*cam_file_to_homography(char *srcFile)
{
  int		ncols;
  int		nrows;
  FILE		*file;
  CamMatrix	*res;

  file = fopen(srcFile, "r");
  if (!file)
    {
      printf("cam_file_to_homography : unable to open file %s\n", srcFile);
      exit (-1);
    }
  if (!fread(&ncols, sizeof(int), 1, file) || !fread(&nrows, sizeof(int), 1, file))
    {
      printf("cam_file_to_homography : incorrect homography header (%s)\n", srcFile);
      exit (-1);
    }
    
  res = (CamMatrix *)malloc(sizeof(CamMatrix));
  cam_allocate_matrix(res, ncols, nrows);
  if (fread(res->data, sizeof(POINTS_TYPE), ncols * nrows, file) != (size_t)(ncols * nrows))
    {
      printf("cam_file_to_homography : incorrect homography data (%s)\n", srcFile);
      exit (-1);
    }

  fclose(file);
  return (res);
}

CamList		*cam_load_points(char *srcFile)
{
  CamList	*res;
  int		nbMatches;
  FILE		*file;
  int		index;
  PointsMatch	tmp;  

  file = fopen(srcFile, "r");
  if (!file)
    {
      printf("cam_load_points : unable to open %s\n", srcFile);
      exit (-1);
    }
  if (!fread(&nbMatches, sizeof(int), 1, file))
    {
      printf("cam_load_points : incorrecct header (%s)\n", srcFile);
      exit (-1);
    }
  index = 0;
  res = NULL;
  while (index < nbMatches)
    {
      res = cam_add_to_linked_list(res, (PointsMatch *)malloc(sizeof(PointsMatch)));
      if (fread(res->data, sizeof(tmp), 1, file) != 1)
	{
	  printf("cam_load_points : %s is not a valid matches file\n", srcFile);
	  exit (-1);
	}
      ++index;
    }
  fclose (file);
  return (res);
}

POINTS_TYPE	cam_euclidian_distance(POINTS_TYPE dx, POINTS_TYPE dy)
{
  return (sqrt((dx * dx) + (dy * dy)));
}

int error_cmp(const void *a, const void *b)
{
  if (*(POINTS_TYPE*)(a) < *(POINTS_TYPE*)(b))
    return (-1);
  if (*(POINTS_TYPE*)(a) == *(POINTS_TYPE*)(b))
    return (0);
  return (1);
}

POINTS_TYPE		*cam_compute_tracking_errors(CamMatrix *H, CamList *points)
{
  CamMatrix		pt1;
  CamMatrix		pt2;
  CamList		*ptr;
  POINTS_TYPE		dx;
  POINTS_TYPE		dy;
  POINTS_TYPE		*err;
  int			index;

  cam_allocate_matrix(&pt1, 1, 3);
  cam_allocate_matrix(&pt2, 1, 3);
  cam_matrix_set_value(&pt1, 0, 2, 1.0f);
  ptr = points;
  index = 0;
  err = (POINTS_TYPE *)malloc(points->index * sizeof(POINTS_TYPE));
  while (ptr)
    {
      cam_matrix_set_value(&pt1, 0, 0, ((PointsMatch *)ptr->data)->pt1.x);
      cam_matrix_set_value(&pt1, 0, 1, ((PointsMatch *)ptr->data)->pt1.y);
      cam_matrix_multiply(&pt2, H, &pt1);
      dx = ABSF(cam_matrix_get_value(&pt2, 0, 0) - ((PointsMatch *)ptr->data)->pt2.x);
      dy = ABSF(cam_matrix_get_value(&pt2, 0, 1) - ((PointsMatch *)ptr->data)->pt2.y);
      err[index] = cam_euclidian_distance(dx, dy);
      ++index;
      ptr = ptr->next;
    }
  cam_disallocate_matrix(&pt1);
  cam_disallocate_matrix(&pt2);
  return (err);
}

void	cam_errors_to_file(char *dir, char *fileName, POINTS_TYPE *err, int nb)
{
  FILE	*file;
  char	*outputPath;
  int	index;

  outputPath = (char *)malloc((strlen(dir) + strlen(fileName) + strlen("err") + 3) * sizeof(char) );
  sprintf(outputPath,"%s/%s.%s", dir, fileName, "err");
  file = fopen(outputPath, "w+");
  index = 0;
  while (index < nb)
    {
      fprintf(file, "%f\n", err[index]);
      ++index;
    }
  free(outputPath);
  fclose(file);
}

int		main()
{
  CamMatrix	*H;
  CamList	*points;
  POINTS_TYPE	*err;

  H = cam_file_to_homography("data/tracking/homo1.tr");
  points = cam_load_points("data/tracking/matches0.matches");
  err = cam_compute_tracking_errors(H, points);
  qsort(err, points->index, sizeof(POINTS_TYPE), error_cmp);
  cam_errors_to_file("data/tracking","errors",err, points->index);

  cam_disallocate_matrix(H);
  free(H);
  cam_disallocate_linked_list(points);
  free (err);

  return (0);
}
