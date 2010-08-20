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
#include "cam_p_from_f.h"

/* TODO : test me */

/**********************/
/* solves	      */
/* [f1 f2 f3]         */
/* [f4 f5 f6] * e = 0 */
/* [f7 f8 f9]         */
/**********************/

CamMatrix	*cam_compute_epipole(CamMatrix *f)
{
  CamMatrix	*e;
  POINTS_TYPE	x;
  POINTS_TYPE	y;
  POINTS_TYPE	f1, f1c;
  POINTS_TYPE	f2, f2c;
  POINTS_TYPE	f3, f3c;
  POINTS_TYPE	f4, f4c;
  POINTS_TYPE	f5, f5c;
  POINTS_TYPE	f6, f6c;
  POINTS_TYPE	f7c;
  POINTS_TYPE	f8c;
  POINTS_TYPE	f9c;
 
  e = (CamMatrix *)malloc(sizeof(CamMatrix));
  cam_allocate_matrix(e, 1, 3);
  f1c = cam_matrix_get_value(f, 0, 0);
  f2c = cam_matrix_get_value(f, 1, 0);
  f3c = cam_matrix_get_value(f, 2, 0);
  f4c = cam_matrix_get_value(f, 0, 1);
  f5c = cam_matrix_get_value(f, 1, 1);
  f6c = cam_matrix_get_value(f, 2, 1);
  f7c = cam_matrix_get_value(f, 0, 2);
  f8c = cam_matrix_get_value(f, 1, 2);
  f9c = cam_matrix_get_value(f, 2, 2);

  /* 1&2 */
  if (ABSF(f1c * f5c - f2c * f4c) > SUP0)
    {
      f1 = f1c;
      f2 = f2c;
      f3 = f3c;
      f4 = f4c;
      f5 = f5c;
      f6 = f6c;
    }
  /* 1&3 */
  else if (ABSF(f1c * f8c - f2c * f7c) > SUP0)
    {
      f1 = f1c;
      f2 = f2c;
      f3 = f3c;
      f4 = f7c;
      f5 = f8c;
      f6 = f9c;
    }
  /* 2&3 */
  else if (ABSF(f4c * f8c - f5c * f7c) > SUP0)
    {
      f1 = f4c;
      f2 = f5c;
      f3 = f6c;
      f4 = f7c;
      f5 = f8c;
      f6 = f9c;
    }
  else
    {
      printf("cam_compute_epipole : unable to determine y\n");
      exit (-1);
    }

  y = (f3 * f4 - f1 * f6) / (f1 * f5 - f2 * f4);

  if (ABSF(f1c) > SUP0)
    {
      f1 = f1c;
      f2 = f2c;
      f3 = f3c;
    }
  else if (ABSF(f4c) > SUP0)
    {
      f1 = f4c;
      f2 = f5c;
      f3 = f6c;
    }
  else if (ABSF(f7c) > SUP0)
    {
      f1 = f7c;
      f2 = f8c;
      f3 = f9c;
    }
  else
    {
      printf("cam_compute_epipole : unable to determine x\n");
      exit (-1);
    }
  x = (-f2 * y -f3) / f1;

  cam_matrix_set_value(e, 0, 0, x);
  cam_matrix_set_value(e, 0, 1, y);
  cam_matrix_set_value(e, 0, 2, 1.0f);

  return (e);
}

CamMatrix	*cam_compute_e(CamMatrix *f)
{
  CamMatrix	*res;

  res = cam_compute_epipole(f);
  return (res);
}

CamMatrix	*cam_compute_eprime(CamMatrix *f)
{
  CamMatrix	fTranspose;
  CamMatrix	*res;

  cam_allocate_matrix(&fTranspose, 3, 3);
  cam_matrix_transpose(&fTranspose, f);
  res = cam_compute_epipole(&fTranspose);
  cam_disallocate_matrix(&fTranspose);
  return (res);
}

/* TODO : test me */

void		cam_compute_p_from_f(CamMatrix *f, CamProjectionsPair *p)
{
  CamMatrix	*ePrime;
  CamMatrix	ePrimeSkew;
  CamMatrix	tmp;
  int		i;
  int		j;

  ePrime = cam_compute_eprime(f);
  cam_allocate_matrix(&ePrimeSkew, 3, 3);
  cam_allocate_matrix(&tmp, 3, 3);

  cam_matrix_set_value(&ePrimeSkew, 0, 0, 0.0f);
  cam_matrix_set_value(&ePrimeSkew, 1, 0, -cam_matrix_get_value(ePrime, 0, 2));
  cam_matrix_set_value(&ePrimeSkew, 2, 0, cam_matrix_get_value(ePrime, 0, 1));
  cam_matrix_set_value(&ePrimeSkew, 0, 1, cam_matrix_get_value(ePrime, 0, 2));
  cam_matrix_set_value(&ePrimeSkew, 1, 1, 0.0f);
  cam_matrix_set_value(&ePrimeSkew, 2, 1, -cam_matrix_get_value(ePrime, 0, 0));
  cam_matrix_set_value(&ePrimeSkew, 0, 2, -cam_matrix_get_value(ePrime, 0, 1));
  cam_matrix_set_value(&ePrimeSkew, 1, 2, cam_matrix_get_value(ePrime, 0, 0));
  cam_matrix_set_value(&ePrimeSkew, 2, 2, 0.0f);
  
  cam_matrix_set_value(&p->p1, 0, 0, 1.0f);
  cam_matrix_set_value(&p->p1, 1, 0, 0.0f);
  cam_matrix_set_value(&p->p1, 2, 0, 0.0f);
  cam_matrix_set_value(&p->p1, 0, 1, 0.0f);
  cam_matrix_set_value(&p->p1, 1, 1, 1.0f);
  cam_matrix_set_value(&p->p1, 2, 1, 0.0f);
  cam_matrix_set_value(&p->p1, 0, 2, 0.0f);
  cam_matrix_set_value(&p->p1, 1, 2, 0.0f);
  cam_matrix_set_value(&p->p1, 2, 2, 1.0f);
  cam_matrix_set_value(&p->p1, 3, 0, 0.0f);
  cam_matrix_set_value(&p->p1, 3, 1, 0.0f);
  cam_matrix_set_value(&p->p1, 3, 2, 0.0f);

  cam_matrix_multiply(&tmp, &ePrimeSkew, f);

  for (j = 0 ; j < 3 ; ++j)
    {
      for (i = 0 ; i < 3 ; ++i)
        {
          cam_matrix_set_value(&p->p2, i, j, cam_matrix_get_value(&tmp, i, j));
        }
      cam_matrix_set_value(&p->p2, i, j, -cam_matrix_get_value(ePrime, 0, j));
    }
  cam_disallocate_matrix(&ePrimeSkew);
  cam_disallocate_matrix(ePrime);
  cam_disallocate_matrix(&tmp);
  free(ePrime);
}

void	cam_disallocate_projections_pair(CamProjectionsPair *p)
{
  cam_disallocate_matrix(&p->p1);
  cam_disallocate_matrix(&p->p2);
}
