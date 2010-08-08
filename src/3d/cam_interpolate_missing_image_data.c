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
#include "cam_interpolate_missing_image_data.h"

int	compare_grayscale(const void *a, const void *b)
{
  if (((RGBandGRAY *)a)->grey < ((RGBandGRAY *)b)->grey)
    return (-1);
  if (((RGBandGRAY *)a)->grey == ((RGBandGRAY *)b)->grey)
    return (0);
  return (-1);
}

/* takes median from neightbours */
CamImageMatrix		*cam_interpolate_missing_image_data(CamImageMatrix *img, unsigned char bgR, unsigned char bgG, unsigned char bgB)
{
  int			i;
  int			j;
  CamImageMatrix	*res;
  int			k;
  int			l;
  int			nb;
  RGBandGRAY		pts[8];

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
	      nb = 0;
	      for (k = -1 ; k <= 1 ; ++k)
		{
		  for (l = -1 ; l <= 1 ; ++l)
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
	      qsort(pts, nb, sizeof(POINTS_TYPE), compare_grayscale);
	      cam_matrix_set_value(&res->r, i ,j, pts[nb/2].red);
	      cam_matrix_set_value(&res->g, i ,j, pts[nb/2].green);
	      cam_matrix_set_value(&res->b, i ,j, pts[nb/2].blue);
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
