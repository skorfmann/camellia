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
#include <string.h>
#include "cam_pgm_to_matrix.h"

void	cam_allocate_imagematrix(CamImageMatrix *m, int ncols, int nrows)
{
  cam_allocate_matrix(&m->r, ncols, nrows);
  cam_allocate_matrix(&m->g, ncols, nrows);
  cam_allocate_matrix(&m->b, ncols, nrows);
}

void	cam_disallocate_imagematrix(CamImageMatrix *m)
{
  cam_disallocate_matrix(&m->r);
  cam_disallocate_matrix(&m->g);
  cam_disallocate_matrix(&m->b);
}

CamImageMatrix		*cam_pgm_to_matrix(char *path)
{
  int			height;
  int			width;
  CamImageMatrix	*res;
  FILE			*file;
  char			header[3];
  int			c;
  int			i;
  int			j;

  res = (CamImageMatrix *)malloc(sizeof(CamImageMatrix));
  file = fopen(path, "r");
  if (!file)
    {
      printf("cam_pgm_to_points : unable to open file %s\n", path);
      exit (-1);
    }
  fread(header, sizeof(char), 3, file);
  if (strncmp(header, "P6\n", 3))
    {
      printf("cam_pgm_to_points : %s not a PGM file\n", path);
      exit (-1);
    }

  width = 0;
  c = fgetc(file);
  while (c != (int)(' '))
    {
      width = 10 * width + (c - (int)('0'));
      c = fgetc(file);
    }

  height = 0;
  c = fgetc(file);
  while (c != (int)('\n'))
    {
      height = 10 * height + (c - (int)('0'));
      c = fgetc(file);
    }

  if (fgetc(file) != (int)('2') || fgetc(file) != (int)('5') || fgetc(file) != (int)('5') || fgetc(file) != (int)('\n'))
    {
      printf("cam_pgm_to_points : %s not a valid PGM file\n", path);
      exit (-1);
    }

  cam_allocate_imagematrix(res, width, height);
  for (j = 0 ; j < height ; ++j)
    {
      for (i = 0 ; i < width ; ++i)
	{
	  if ((c = fgetc(file)) != EOF)
	    cam_matrix_set_value(&res->r, i, j, (POINTS_TYPE)c);
	  else
	    {
	      printf("cam_pgm_to_points : malformed pgm file (red) %s (%i, %i)\n", path, i, j);
	      exit (-1);
	    }
	  if ((c = fgetc(file)) != EOF)
	    cam_matrix_set_value(&res->g, i, j, (POINTS_TYPE)c);
	  else
	    {
	      printf("cam_pgm_to_points : malformed pgm file (green) %s (%i, %i)\n", path, i, j);
	      exit (-1);
	    }
	  if ((c = fgetc(file)) != EOF)
	    cam_matrix_set_value(&res->b, i, j, (POINTS_TYPE)c);
	  else
	    {
	      printf("cam_pgm_to_points : malformed pgm file (blue) %s (%i, %i)\n", path, i, j);
	      exit (-1);
	    }
	}
    }
  fclose(file);
  return (res);
}
