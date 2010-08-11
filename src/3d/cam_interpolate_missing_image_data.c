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

#include <stdlib.h>
#include <stdio.h>
#include "misc.h"
#include "cam_2d_points.h"
#include "cam_interpolate_missing_image_data.h"

POINTS_TYPE	cam_round(POINTS_TYPE a)
{
  POINTS_TYPE	b;

  b = (POINTS_TYPE)((int)a);
  if (ABSF(a - b) >= 0.5f)
    return (b + 1.0f);
  return (b);
}

int	compare_grayscale(const void *a, const void *b)
{
  if (((RGBandGRAY *)a)->grey < ((RGBandGRAY *)b)->grey)
    return (-1);
  if (((RGBandGRAY *)a)->grey == ((RGBandGRAY *)b)->grey)
    return (0);
  return (-1);
}

/* takes median from neightbours */
CamImageMatrix		*cam_interpolate_missing_image_data_median(CamImageMatrix *img, int neighborhood, unsigned char bgR, unsigned char bgG, unsigned char bgB)
{
  int			i;
  int			j;
  CamImageMatrix	*res;
  int			k;
  int			l;
  int			nb;
  RGBandGRAY		*pts;

  pts = (RGBandGRAY *)malloc((2*neighborhood+1)*(2*neighborhood+1)*sizeof(RGBandGRAY));
  res = (CamImageMatrix *)malloc(sizeof(CamImageMatrix));
  cam_allocate_imagematrix(res, img->r.ncols, img->r.nrows);
  for (j = neighborhood ; j < img->r.nrows - neighborhood ; ++j)
    {
      for (i = neighborhood ; i < img->r.ncols - neighborhood ; ++i)
	{
	  if ((unsigned char)cam_matrix_get_value(&img->r, i ,j) == bgR &&
	      (unsigned char)cam_matrix_get_value(&img->g, i ,j) == bgG &&
	      (unsigned char)cam_matrix_get_value(&img->b, i ,j) == bgB)
	    {
	      nb = 0;
	      for (k = -neighborhood ; k <= neighborhood ; ++k)
		{
		  for (l = -neighborhood ; l <= neighborhood ; ++l)
		    {
		      if (k || l)
			{
			  if ((unsigned char)cam_matrix_get_value(&img->r, i + k,j + l) != bgR &&
			      (unsigned char)cam_matrix_get_value(&img->g, i + k,j + l) != bgG &&
			      (unsigned char)cam_matrix_get_value(&img->b, i + k,j + l) != bgB)
			    {
			      pts[nb].red = cam_matrix_get_value(&img->r, i + k, j + l);
			      pts[nb].green = cam_matrix_get_value(&img->g, i + k, j + l);
			      pts[nb].blue = cam_matrix_get_value(&img->b, i + k, j + l);
			      pts[nb].grey = pts[nb].red * 0.3f + pts[nb].green * 0.59f + pts[nb].blue * 0.11f;
			      ++nb;
			    }
			}
		    }
		}
	      if (!nb)
		{
		  cam_matrix_set_value(&res->r, i ,j, cam_matrix_get_value(&img->r, i, j));
		  cam_matrix_set_value(&res->g, i ,j, cam_matrix_get_value(&img->g, i, j));
		  cam_matrix_set_value(&res->b, i ,j, cam_matrix_get_value(&img->b, i, j));
		}
	      else
		{
		  qsort(pts, nb, sizeof(POINTS_TYPE), compare_grayscale);
		  cam_matrix_set_value(&res->r, i ,j, pts[nb/2].red);
		  cam_matrix_set_value(&res->g, i ,j, pts[nb/2].green);
		  cam_matrix_set_value(&res->b, i ,j, pts[nb/2].blue);
		}
	    }
	  else
	    {
	      cam_matrix_set_value(&res->r, i ,j, cam_matrix_get_value(&img->r, i, j));
	      cam_matrix_set_value(&res->g, i ,j, cam_matrix_get_value(&img->g, i, j));
	      cam_matrix_set_value(&res->b, i ,j, cam_matrix_get_value(&img->b, i, j));
	    }
	}
    }
  free (pts);
  return (res);
}

CamImageMatrix		*cam_interpolate_missing_image_data_back_projection(CamImageMatrix *img, CamImageMatrix *origin, CamMatrix *inverseHomography, unsigned char bgR, unsigned char bgG, unsigned char bgB)
{
  CamImageMatrix	*res;
  int			i;
  int			j;
  CamMatrix		pt1;
  CamMatrix		pt2;
  int			roundX;
  int			roundY;

  res = (CamImageMatrix *)malloc(sizeof(CamImageMatrix));
  cam_allocate_imagematrix(res, img->r.ncols, img->r.nrows);
  cam_allocate_matrix(&pt1, 1, 3);
  cam_allocate_matrix(&pt2, 1, 3);
  cam_matrix_set_value(&pt1, 0 , 2 , 1.0f);
  for (j = 0 ; j < img->r.nrows ; ++j)
    {
      for (i = 0 ; i < img->r.ncols ; ++i)
	{
	  if ((unsigned char)cam_matrix_get_value(&img->r, i ,j) == bgR &&
	      (unsigned char)cam_matrix_get_value(&img->g, i ,j) == bgG &&
	      (unsigned char)cam_matrix_get_value(&img->b, i ,j) == bgB)
	    {
	      cam_matrix_set_value(&pt1, 0, 0, (POINTS_TYPE)(i - img->r.ncols / 2));
	      cam_matrix_set_value(&pt1, 0, 1, (POINTS_TYPE)(j - img->r.nrows / 2));
	      cam_matrix_multiply(&pt2, inverseHomography, &pt1);
	      roundX = (int)cam_round(cam_matrix_get_value(&pt2, 0, 0)) + img->r.ncols / 2;
	      roundY = (int)cam_round(cam_matrix_get_value(&pt2, 0, 1)) + img->r.nrows / 2;
	      if (roundX >= 0 && roundX < img->r.ncols && roundY >= 0 && roundY < img->r.nrows)
		{
		  cam_matrix_set_value(&res->r, i, j, cam_matrix_get_value(&origin->r, roundX, roundY));
		  cam_matrix_set_value(&res->g, i, j, cam_matrix_get_value(&origin->g, roundX, roundY));
		  cam_matrix_set_value(&res->b, i, j, cam_matrix_get_value(&origin->b, roundX, roundY));
		}
	      else
	      {
		cam_matrix_set_value(&res->r, i ,j, (POINTS_TYPE)bgR);
		cam_matrix_set_value(&res->g, i ,j, (POINTS_TYPE)bgG);
		cam_matrix_set_value(&res->b, i ,j, (POINTS_TYPE)bgB);
	      }
	    }
	  else
	    {
	      cam_matrix_set_value(&res->r, i ,j, cam_matrix_get_value(&img->r, i, j));
	      cam_matrix_set_value(&res->g, i ,j, cam_matrix_get_value(&img->g, i, j));
	      cam_matrix_set_value(&res->b, i ,j, cam_matrix_get_value(&img->b, i, j));
	    }
	}
    }
  cam_disallocate_matrix(&pt1);
  cam_disallocate_matrix(&pt2);
  return (res);
}

/* following : http://en.wikipedia.org/wiki/Bilinear_interpolation */
/* a--b */
/* |  | */
/* c--d */
CamColorized2dPoint	*cam_bilinear_interpolation(Cam2dPoint *dst, CamColorized2dPoint *a, CamColorized2dPoint *b, CamColorized2dPoint *c, CamColorized2dPoint *d)
{
  CamColorized2dPoint	*res;
  POINTS_TYPE		x1;
  POINTS_TYPE		x2;
  POINTS_TYPE		y1;
  POINTS_TYPE		y2;
  POINTS_TYPE		x;
  POINTS_TYPE		y;
  POINTS_TYPE		denominator;

  x = dst->x;
  y = dst->y;
  x1 = a->point.x;
  x2 = b->point.x;
  y1 = d->point.y;
  y2 = b->point.y;
  denominator = (x2 - x1) * (y2 - y1);
  res = (CamColorized2dPoint *)malloc(sizeof(CamColorized2dPoint));
  res->point.x = dst->x;
  res->point.y = dst->y;
  res->color.r = c->color.r * (x2 - x) * (y2 - y) / denominator +
    d->color.r * (x - x1) * (y2 - y) / denominator+
    a->color.r * (x2 - x) * (y - y1) / denominator+
    b->color.r * (x - x1) * (y - y1) / denominator;
  res->color.g = c->color.g * (x2 - x) * (y2 - y) / denominator +
    d->color.g * (x - x1) * (y2 - y) / denominator+
    a->color.g * (x2 - x) * (y - y1) / denominator+
    b->color.g * (x - x1) * (y - y1) / denominator;
  res->color.b = c->color.b * (x2 - x) * (y2 - y) / denominator +
    d->color.b * (x - x1) * (y2 - y) / denominator+
    a->color.b * (x2 - x) * (y - y1) / denominator+
    b->color.b * (x - x1) * (y - y1) / denominator;

  return (res);
}

CamImageMatrix		*cam_interpolate_missing_image_data_bilinear(CamImageMatrix *img, unsigned char bgR, unsigned char bgG, unsigned char bgB)
{
  CamColorized2dPoint	a, b, c, d, *interpoled;
  CamImageMatrix	*res;
  Cam2dPoint		pt;
  int			i;
  int			j;
  
  res = (CamImageMatrix *)malloc(sizeof(CamImageMatrix));
  cam_allocate_imagematrix(res, img->r.ncols, img->r.nrows);
  for (j = 1 ; j < img->r.nrows - 1 ; ++j)
    {
      for (i = 1 ; i < img->r.ncols - 1 ; ++i)
	{
	  if ((unsigned char)cam_matrix_get_value(&img->r, i ,j) == bgR &&
	      (unsigned char)cam_matrix_get_value(&img->g, i ,j) == bgG &&
	      (unsigned char)cam_matrix_get_value(&img->b, i ,j) == bgB)
	    {
	      if ((unsigned char)cam_matrix_get_value(&img->r, i - 1, j - 1) != bgR &&
		  (unsigned char)cam_matrix_get_value(&img->g, i - 1, j - 1) != bgG &&
		  (unsigned char)cam_matrix_get_value(&img->b, i - 1, j - 1) != bgB &&

		  (unsigned char)cam_matrix_get_value(&img->r, i - 1,j+1) != bgR &&
		  (unsigned char)cam_matrix_get_value(&img->g, i - 1,j+1) != bgG &&
		  (unsigned char)cam_matrix_get_value(&img->b, i - 1,j+1) != bgB &&

		  (unsigned char)cam_matrix_get_value(&img->r, i + 1,j-1) != bgR &&
		  (unsigned char)cam_matrix_get_value(&img->g, i + 1,j-1) != bgG &&
		  (unsigned char)cam_matrix_get_value(&img->b, i +1,j-1) != bgB &&

		  (unsigned char)cam_matrix_get_value(&img->r, i +1,j+1) != bgR &&
		  (unsigned char)cam_matrix_get_value(&img->g, i +1,j+1) != bgG &&
		  (unsigned char)cam_matrix_get_value(&img->b, i +1,j+1) != bgB)
		{
		  /*printf("%i %i\n", i, j);*/
		  pt.x = i;
		  pt.y = j;
		  c.point.x = i - 1;
		  c.point.y = j - 1;
		  c.color.r = (unsigned char)cam_matrix_get_value(&img->r, (int)c.point.x , (int)c.point.y);
		  c.color.g = (unsigned char)cam_matrix_get_value(&img->g, (int)c.point.x , (int)c.point.y);
		  c.color.b = (unsigned char)cam_matrix_get_value(&img->b, (int)c.point.x , (int)c.point.y);
		  /*printf("c : %f %f // %i %i %i\n", c.point.x, c.point.y, c.color.r, c.color.g, c.color.b);*/
		  b.point.x = i + 1;
		  b.point.y = j + 1;
		  b.color.r = (unsigned char)cam_matrix_get_value(&img->r, (int)b.point.x , (int)c.point.y);
		  b.color.g = (unsigned char)cam_matrix_get_value(&img->g, (int)b.point.x , (int)c.point.y);
		  b.color.b = (unsigned char)cam_matrix_get_value(&img->b, (int)b.point.x , (int)c.point.y);
		  /*printf("b : %f %f // %i %i %i\n", b.point.x, b.point.y, b.color.r, b.color.g, b.color.b);*/
		  a.point.x = i - 1;
		  a.point.y = j + 1;
		  a.color.r = (unsigned char)cam_matrix_get_value(&img->r, (int)a.point.x , (int)a.point.y);
		  a.color.g = (unsigned char)cam_matrix_get_value(&img->g, (int)a.point.x , (int)c.point.y);
		  a.color.b = (unsigned char)cam_matrix_get_value(&img->b, (int)a.point.x , (int)a.point.y);
		  /*printf("a : %f %f // %i %i %i\n", a.point.x, a.point.y, a.color.r, a.color.g, a.color.b);*/
		  d.point.x = i + 1;
		  d.point.y = j - 1;
		  d.color.r = (unsigned char)cam_matrix_get_value(&img->r, (int)d.point.x , (int)d.point.y);
		  d.color.g = (unsigned char)cam_matrix_get_value(&img->g, (int)d.point.x , (int)d.point.y);
		  d.color.b = (unsigned char)cam_matrix_get_value(&img->b, (int)d.point.x , (int)d.point.y);
		  /*printf("d : %f %f // %i %i %i\n", d.point.x, d.point.y, d.color.r, d.color.g, d.color.b);*/
		  interpoled = cam_bilinear_interpolation(&pt, &a, &b, &c, &d);
		  /*printf("res : %i %i %i\n", interpoled->color.r, interpoled->color.g, interpoled->color.b);*/
		  cam_matrix_set_value(&res->r, i ,j, (POINTS_TYPE)interpoled->color.r);
		  cam_matrix_set_value(&res->g, i ,j, (POINTS_TYPE)interpoled->color.g);
		  cam_matrix_set_value(&res->b, i ,j, (POINTS_TYPE)interpoled->color.b);
		  free (interpoled);
		}
	      else
		{
		  printf ("here\n");
		  cam_matrix_set_value(&res->r, i ,j, cam_matrix_get_value(&img->r, i, j));
		  cam_matrix_set_value(&res->g, i ,j, cam_matrix_get_value(&img->g, i, j));
		  cam_matrix_set_value(&res->b, i ,j, cam_matrix_get_value(&img->b, i, j));
		}
	    } 
	  else
	    {
	      cam_matrix_set_value(&res->r, i ,j, cam_matrix_get_value(&img->r, i, j));
	      cam_matrix_set_value(&res->g, i ,j, cam_matrix_get_value(&img->g, i, j));
	      cam_matrix_set_value(&res->b, i ,j, cam_matrix_get_value(&img->b, i, j));
	    }
	}
    }
  return (res);
}
