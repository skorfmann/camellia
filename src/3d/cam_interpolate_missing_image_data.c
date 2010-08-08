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

/* takes average from neightbours */
void		cam_interpolate_missing_image_data(CamImageMatrix *img, unsigned char bgR, unsigned char bgG, unsigned char bgB)
{
  int		i;
  int		j;
  int		k;
  int		l;
  POINTS_TYPE	avgR;
  POINTS_TYPE	avgG;
  POINTS_TYPE	avgB;

  for (j = 1 ; j < img->r.nrows - 1 ; ++j)
    {
      for (i = 1 ; i < img->r.ncols - 1 ; ++i)
	{
	  if ((unsigned char)cam_matrix_get_value(&img->r, i ,j) == bgR &&
	      (unsigned char)cam_matrix_get_value(&img->g, i ,j) == bgG &&
	      (unsigned char)cam_matrix_get_value(&img->b, i ,j) == bgB)
	    {
	      avgR = 0.0f;
	      avgG = 0.0f;
	      avgB = 0.0f;
	      for (k = -1 ; k <= 1 ; ++k)
		{
		  for (l = -1 ; l <= 1 ; ++l)
		    {
		      if (k || l)
			{
			  if ((unsigned char)cam_matrix_get_value(&img->r, i + k,j + l) != bgR &&
			      (unsigned char)cam_matrix_get_value(&img->g, i + k,j + l) != bgG &&
			      (unsigned char)cam_matrix_get_value(&img->b, i + k,j + l) != bgB)
			  avgR += cam_matrix_get_value(&img->r, i + k, j + l);
			  avgG += cam_matrix_get_value(&img->r, i + k, j + l);
			  avgB += cam_matrix_get_value(&img->r, i + k, j + l);
			}
		    }
		}
	      avgR /= 8.0f;
	      avgG /= 8.0f;
	      avgB /= 8.0f;
	      cam_matrix_set_value(&img->r, i ,j, avgR);
	      cam_matrix_set_value(&img->g, i ,j, avgG);
	      cam_matrix_set_value(&img->b, i ,j, avgB);
	    }
	}
    }
    
}
