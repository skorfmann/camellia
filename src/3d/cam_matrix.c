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

void		cam_matrix_convolution(CamMatrix *dst, CamMatrix *src, CamMatrix *mask, POINTS_TYPE factor)
{
  int		i;
  int		j;
  int		k;
  int		l;
  POINTS_TYPE	tmp;

  if (!(mask->ncols % 2) || !(mask->nrows % 2))
    {
      printf("cam_convolution_matrix : incorrect convolution matrix\n");
      exit (-1);
    }
  if (src->ncols < mask->ncols || src->nrows < mask->ncols)
    {
      printf("cam_convolution_matrix : matrix smaller than convolution mask\n");
      exit (-1);
    }
  if  ((src->nrows != dst->nrows) || (src->ncols != dst->ncols))
    {
      printf("cam_convolution_matrix : source matrix and destination matrix have different sizes\n");
      exit (-1);
    }
  for (j = mask->nrows / 2 ; j < src->nrows - mask->nrows / 2; ++j)
    {  
      for (i = mask->ncols / 2 ; i < src->ncols - mask->ncols / 2 ; ++i)
	{
	  tmp = 0.0f;
	  for (l = 0 ; l < mask->nrows ; ++l)
	    {
	      for (k = 0 ; k < mask->ncols ; ++k)
		{
		  tmp += cam_matrix_get_value(src, i - mask->ncols / 2 + k, j - mask->nrows / 2 +l) * cam_matrix_get_value(mask, k, l);
		}
	    }
	  tmp /= factor;
	  cam_matrix_set_value(dst, i, j, tmp);
	}
    }
}

void		cam_matrix_image_convolution(CamImageMatrix *dst, CamImageMatrix *src, CamMatrix *mask, POINTS_TYPE factor)
{
  cam_matrix_convolution(&dst->r, &src->r, mask, factor);
  cam_matrix_convolution(&dst->g, &src->g, mask, factor);
  cam_matrix_convolution(&dst->b, &src->b, mask, factor);
}

void	cam_allocate_matrix(CamMatrix *m, int ncols, int nrows)
{
  m->ncols = ncols;
  m->nrows = nrows;
  if (!(ncols * nrows))
    {
      printf("cam_allocate_matrix : dimension error (%i,%i)\n", ncols, nrows);
      exit (-1);
    }
  m->data = (POINTS_TYPE *)calloc(nrows * ncols, sizeof(POINTS_TYPE));
  if (!m->data)
    {
      printf("cam_allocate_matrix : alloc error (%i,%i)\n", ncols, nrows);
      exit (-1);
    }
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

void	cam_matrix_add(CamMatrix *res, CamMatrix *m1, CamMatrix *m2)
{
  register int	i;
  register int	j;
  
  if (m1->ncols != m2->ncols || m1->nrows != m2->nrows)
    {
      printf("CamMatrixAdd : m1->ncols != m2->cols || m1->nrows != m2->nrows");
      exit(-1);
    }
  for (j = 0 ; j < m1->nrows ; ++j)
    {
      for (i = 0 ; i < m1->ncols ; ++i)
	{
	  cam_matrix_set_value(res, i, j, cam_matrix_get_value(m1, i, j) + cam_matrix_get_value(m2, i, j));
	}
    }
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
      printf("CamMatrixMultiply : m1->ncols != m2->nrows\n");
      exit(-1);
    }
  if (m1->nrows != res->nrows || m2->ncols != res->ncols)
    {
      printf("CamMatrixMultiply : m1->nrows != res->nrows || m2->ncols != res->ncols\n");
      exit(-1);
    }
  for (j = 0 ; j < res->nrows ; ++j)
    {
      for (i = 0 ; i < res->ncols ; ++i)
	{
	  cam_matrix_set_value(res, i, j, 0.0f);
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

void		cam_matrix_transpose(CamMatrix *dst, CamMatrix *src)
{
  register int	i;
  register int	j;

  if (dst->ncols != src->nrows || dst->nrows != src->ncols)
    {
      printf("cam_matrix_transpose : dst->ncols != src->nrows || dst->nrows != src->ncols\n");
      exit (-1);
    }
  for (j = 0 ; j < src->nrows ; ++j)
    for (i = 0 ; i < src->ncols ; ++i)
      {
	cam_matrix_set_value(dst, j, i, cam_matrix_get_value(src, i , j));
      }
}
