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
      printf("cam_load_points : incorrect header (%s)\n", srcFile);
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

POINTS_TYPE		*cam_compute_tracking_endpoint_errors(int width, int height, CamMatrix *H, CamList *points)
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
      cam_matrix_set_value(&pt1, 0, 0, ((PointsMatch *)ptr->data)->pt1.x - (POINTS_TYPE)(width / 2));
      cam_matrix_set_value(&pt1, 0, 1, ((PointsMatch *)ptr->data)->pt1.y - (POINTS_TYPE)(height / 2));
      cam_matrix_multiply(&pt2, H, &pt1);
      dx = ABSF(cam_matrix_get_value(&pt2, 0, 0) - (((PointsMatch *)ptr->data)->pt2.x - (POINTS_TYPE)(width / 2)));
      dy = ABSF(cam_matrix_get_value(&pt2, 0, 1) - (((PointsMatch *)ptr->data)->pt2.y - (POINTS_TYPE)(height / 2)));
      err[index] = cam_euclidian_distance(dx, dy);
      ++index;
      ptr = ptr->next;
    }
  cam_disallocate_matrix(&pt1);
  cam_disallocate_matrix(&pt2);
  return (err);
}

POINTS_TYPE		*cam_compute_tracking_angular_errors(int width, int height, CamMatrix *H, CamList *points)
{
  CamMatrix		pt1;
  CamMatrix		pt2;
  CamList		*ptr;
  POINTS_TYPE		u, v;
  POINTS_TYPE		uGT, vGT;
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
      cam_matrix_set_value(&pt1, 0, 0, ((PointsMatch *)ptr->data)->pt1.x - (POINTS_TYPE)(width / 2));
      cam_matrix_set_value(&pt1, 0, 1, ((PointsMatch *)ptr->data)->pt1.y - (POINTS_TYPE)(height / 2));
      cam_matrix_multiply(&pt2, H, &pt1);
      u = (((PointsMatch *)ptr->data)->pt2.x - (POINTS_TYPE)(width / 2)) - (((PointsMatch *)ptr->data)->pt1.x - (POINTS_TYPE)(width / 2));
      v = (((PointsMatch *)ptr->data)->pt2.y - (POINTS_TYPE)(height / 2)) - (((PointsMatch *)ptr->data)->pt1.y - (POINTS_TYPE)(height / 2));
      uGT = cam_matrix_get_value(&pt2, 0, 0) - (((PointsMatch *)ptr->data)->pt1.x - (POINTS_TYPE)(width / 2));
      vGT = cam_matrix_get_value(&pt2, 0, 1) - (((PointsMatch *)ptr->data)->pt1.y - (POINTS_TYPE)(height / 2));
      err[index] = acos( (1.0f + u * uGT + v * vGT) /  (sqrt(1.0f + u * u + v * v) * sqrt(1.0f + uGT * uGT + vGT * vGT) + 0.00001) );
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

void		write_scales_histogram(char *path, int index, int nb)
{
  FILE		*file;
  int		i;
  char		file_name[50];
  int		*tab;
  int		max;
  POINTS_TYPE	*histo;

  tab = malloc(nb * sizeof(int));
  sprintf(file_name, "%s/scales%i.scales", path, index);
  file = fopen(file_name, "r");
  if (!file)
    {
      printf("write_scales_histogram : unable to open %s\n", file_name);
      exit (-1);
    }
  i = 0;
  max = 0;
  while (i < nb)
    {
      if (fread(tab + i, sizeof(*tab), 1, file) != 1)
	{
	  printf("write_scales_histogram : %s is not a valid matches file\n", file_name);
	  exit (-1);
	}
      tab[i] = tab[i] / 4;
      if (tab[i] > max)
	max = tab[i];
      ++i;
    }
  histo = malloc(max* sizeof(*histo));
  i = 0;
  while (i < max)
    {
      histo[i] = 0.0f;
      ++i;
    }
  i = 0;
  while (i < nb)
    {
      histo[tab[i] - 1] += 1.0f;
      ++i;
    }
  fclose (file);
  sprintf(file_name, "%s/histo_scales%i.histo", path, index);
  file = fopen(file_name, "w+");
  i = 0;
  while (i < max)
    {
      fprintf(file, "%f\n", histo[i]);
      ++i;
    }
  fclose(file);
  free(histo);
  free(tab);
}

void	print_evaluate_tracking_usage()
{
  printf("./evaluate_tracking number_of_evaluations width height\n");
  printf("Nb : the .tr and .matches files have to be in data/tracking and have their figure starting from 0\n");
  printf("eg : homography0.tr matches0.matches homography1.tr matches1.matches\n");
}

int		main(int ac, char **av)
{
  int		i;
  CamMatrix	*H;
  CamList	*points;
  POINTS_TYPE	*err;
  int		nbFiles;
  int		index;
  int		width;
  int		height;
  char		path[] = "data/tracking";
  char		file[50];
  POINTS_TYPE	treshold = 10;
  int		foo;
  POINTS_TYPE	*angular_averages;
  POINTS_TYPE	*endpoint_averages;
  POINTS_TYPE	*false_positive;

  if (ac != 4)
    {
      print_evaluate_tracking_usage();
      exit (-1);
    }
  nbFiles = atoi(av[1]);
  width = atoi(av[2]);
  height = atoi(av[3]);
  angular_averages = (POINTS_TYPE *)malloc(nbFiles * sizeof(POINTS_TYPE));
  endpoint_averages = (POINTS_TYPE *)malloc(nbFiles * sizeof(POINTS_TYPE));
  false_positive = (POINTS_TYPE *)malloc(nbFiles * sizeof(POINTS_TYPE));
  index = 0;
  while (index < nbFiles)
    {
      sprintf(file, "%s/homography%i.tr", path, index);
      H = cam_file_to_homography(file);
      sprintf(file, "%s/matches%i.matches", path, index);
      points = cam_load_points(file);

      err = cam_compute_tracking_endpoint_errors(width, height, H, points);
      qsort(err, points->index, sizeof(POINTS_TYPE), error_cmp);
      foo = 0;
      endpoint_averages[index] = 0.0f;
      while (err[foo] < treshold && foo < points->index)
	{
	  endpoint_averages[index] += err[foo];
	  ++foo;
	}
      endpoint_averages[index] /= foo;
      sprintf(file, "endpoint_errors%i", index);
      cam_errors_to_file(path, file, err, points->index);
      free (err);

      err = cam_compute_tracking_angular_errors(width, height, H, points);
      qsort(err, points->index, sizeof(POINTS_TYPE), error_cmp);
      angular_averages[index] = 0.0f;
      i = 0;
      while (i < foo)
	{
	  angular_averages[index] += err[i];
	  ++i;
	}
      angular_averages[index] /= foo;
      sprintf(file, "angular_errors%i", index);
      cam_errors_to_file(path, file, err, points->index);
      free (err);

      false_positive[index] = (POINTS_TYPE)((points->index - foo) * 100 / points->index);

      write_scales_histogram(path, index, points->index);

      cam_disallocate_matrix(H);
      free(H);
      cam_disallocate_linked_list(points);
      ++index;
    }

  
  cam_errors_to_file(path, "angular_averages.avg", angular_averages, nbFiles);
  cam_errors_to_file(path, "endpoint_averages.avg", endpoint_averages, nbFiles);
  cam_errors_to_file(path, "false_positive.miss", false_positive, nbFiles);
      
  free(angular_averages);
  free(endpoint_averages);
  free(false_positive);
  return (0);
}
