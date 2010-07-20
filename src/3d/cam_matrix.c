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
#include "cam_matrix.h"

void	cam_allocate_matrix(CamMatrix *m, int ncols, int nrows)
{
  m->ncols = ncols;
  m->nrows = nrows;
  m->data = (POINTS_TYPE *)calloc(nrows * ncols, sizeof(POINTS_TYPE));
}

void	cam_disallocate_matrix(CamMatrix *m)
{
  free(m->data);
  m->data = NULL;
  m->nrows = 0;
  m->ncols = 0;
}

void	cam_matrix_set_value(CamMatrix *m, int x, int y, POINTS_TYPE value)
{
  m->data[y * m->ncols + x] = value;
  return ;
}

POINTS_TYPE	cam_matrix_get_value(CamMatrix *m, int x, int y)
{
  return (m->data[y * m->ncols + x]);
}

void	cam_matrix_add_value(CamMatrix *m, int x, int y, POINTS_TYPE value)
{
  m->data[y * m->ncols + x] += value;
  return ;
}

void		cam_print_matrix(CamMatrix *mat, char *name)
{
  register int	i;
  register int	j;

  if (name)
    printf("%s\n", name);
  for (j = 0 ; j < mat->nrows ; ++j)
    {
      for (i = 0 ; i < mat->ncols ; ++i)
	{
	  printf("%f\t", cam_matrix_get_value(mat, i, j));
	}
      printf("\n");
    }
}

void	cam_matrix_multiply(CamMatrix *res, CamMatrix *m1, CamMatrix *m2)
{
  register int	i;
  register int	j;
  register int	k;

  if (!res->data)
    cam_allocate_matrix(res, m1->nrows, m2->ncols);
  if (m1->ncols != m2->nrows)
    {
      printf("CamMatrixMultiply : m1->ncols != m2->nrows");
      exit(-1);
    }
  for (j = 0 ; j < res->nrows ; ++j)
    {
      for (i = 0 ; i < res->ncols ; ++i)
	{
	  for (k = 0 ; k < m1->ncols ; ++k)
	    {
	      cam_matrix_add_value(res, i, j, cam_matrix_get_value(m1, k, j) * cam_matrix_get_value(m2, i, k));
	    }
	}
    }
}

void	cam_matrix_copy(CamMatrix *dst, CamMatrix *src)
{
  register int	i;
  register int	j;

  if (dst->ncols != src->ncols || dst->nrows != src->nrows)
    {
      printf("CamMatrixCopy : m1->ncols != m2->cols || m1->nrows != m2->nrows");
      exit (-1);
    }
  for (j = 0 ; j < src->nrows ; ++j)
    {
      for (i = 0 ; i < src->ncols ; ++i)
	{
	  cam_matrix_set_value(dst, i, j, cam_matrix_get_value(src, i, j));
	}
    }
}
