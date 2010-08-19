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
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_poly.h>
#include "cam_reprojection_error_minimization.h"

/* implementation algorithm p.318 in multiple view geometry in computer vision */

/* begin step 1 */
CamMatrix	*cam_transformation_matrix(CamMatrix *pt)
{
  CamMatrix	*res;

  res = (CamMatrix *)malloc(sizeof(CamMatrix));
  cam_allocate_matrix(res, 3, 3);
  cam_matrix_set_value(res, 0, 0, 1.0f);
  cam_matrix_set_value(res, 0, 1, 0.0f);
  cam_matrix_set_value(res, 0, 2, 0.0f);
  cam_matrix_set_value(res, 1, 0, 0.0f);
  cam_matrix_set_value(res, 1, 1, 1.0f);
  cam_matrix_set_value(res, 1, 2, 0.0f);
  cam_matrix_set_value(res, 2, 0, -cam_matrix_get_value(pt, 0, 0));
  cam_matrix_set_value(res, 2, 1, -cam_matrix_get_value(pt, 0, 1));
  cam_matrix_set_value(res, 2, 2, 1.0f);
  return (res);
}
/* end step 1 */

/* begin step 2 */
void		cam_transform_f_with_transformation(CamMatrix *F, CamMatrix *T, CamMatrix *Tprime)
{
  CamMatrix	tmp;
  CamMatrix	Tinverse;
  CamMatrix	Tprimeinverse;
  CamMatrix	Tprimeinversetranspose;

  cam_allocate_matrix(&tmp, 3, 3);
  cam_allocate_matrix(&Tinverse, 3, 3);
  cam_allocate_matrix(&Tprimeinverse, 3, 3);
  cam_allocate_matrix(&Tprimeinversetranspose, 3, 3);
  
  cam_matrix_inverse_3x3(&Tinverse, T);
  cam_matrix_inverse_3x3(&Tprimeinverse, Tprime);
  cam_matrix_transpose(&Tprimeinversetranspose, &Tprimeinverse);
  
  cam_matrix_multiply(&tmp, &Tprimeinversetranspose, F);
  cam_matrix_multiply(F, &tmp, &Tinverse);
  
  cam_disallocate_matrix(&tmp);
  cam_disallocate_matrix(&Tinverse);
  cam_disallocate_matrix(&Tprimeinverse);
  cam_disallocate_matrix(&Tprimeinversetranspose);
}
/* end step 2 */

/* begin step 3 */
void		cam_normalize_epipole(CamMatrix *e)
{
  POINTS_TYPE	s;

  s = sqrt(1.0f / (cam_matrix_get_value(e, 0, 0) * cam_matrix_get_value(e, 0, 0) +
		   cam_matrix_get_value(e, 0, 1) * cam_matrix_get_value(e, 0, 1)));
  cam_matrix_set_value(e, 0, 0, s * cam_matrix_get_value(e, 0, 0));
  cam_matrix_set_value(e, 0, 1, s * cam_matrix_get_value(e, 0, 1));
}

void	cam_normalize_epipoles(CamEpipoles *e)
{
  cam_normalize_epipole(&e->eFirst);
  cam_normalize_epipole(&e->eSecond);
}

CamEpipoles	*cam_compute_normalized_epipoles(CamMatrix *f)
{
  CamEpipoles	*res;

  res = cam_compute_epipoles(f);
  cam_normalize_epipoles(res);
  return (res);
}
/* end step 3 */

/* begin step 4 */
CamMatrix	*cam_form_rotation_matrix(CamMatrix *epipole)
{
  CamMatrix	*res;

  res = malloc(sizeof(CamMatrix));
  cam_allocate_matrix(res, 3 , 3);
  cam_matrix_set_value(res, 0, 0, cam_matrix_get_value(epipole, 0, 0));
  cam_matrix_set_value(res, 1, 0, cam_matrix_get_value(epipole, 0, 1));
  cam_matrix_set_value(res, 2, 0, 0.0f);
  cam_matrix_set_value(res, 0, 1, -cam_matrix_get_value(epipole, 0, 1));
  cam_matrix_set_value(res, 1, 1, cam_matrix_get_value(epipole, 0, 0));
  cam_matrix_set_value(res, 2, 1, 0.0f);
  cam_matrix_set_value(res, 0, 2, 0.0f);
  cam_matrix_set_value(res, 1, 2, 0.0f);
  cam_matrix_set_value(res, 2, 2, 1.0f);
  return (res);
}
/* end step 4 */

/* begin step 5 */
void		cam_transform_f_with_rotation(CamMatrix *F, CamMatrix *R, CamMatrix *Rprime)
{
  CamMatrix	tmp;
  CamMatrix	Rtranspose;

  cam_allocate_matrix(&tmp, 3, 3);
  cam_allocate_matrix(&Rtranspose, 3, 3);
  cam_matrix_transpose(&Rtranspose, R);
  cam_matrix_multiply(&tmp, Rprime, F);
  cam_matrix_multiply(F, &tmp, &Rtranspose);
  cam_disallocate_matrix(&tmp);
  cam_disallocate_matrix(&Rtranspose);
}
/* end step 5 */

/* begin step 6 */
void		cam_replace_data_in_f(CamMatrix *F, CamMatrix *e, CamMatrix *eprime)
{
  POINTS_TYPE	f;
  POINTS_TYPE	fprime;
  POINTS_TYPE	a;
  POINTS_TYPE	b;
  POINTS_TYPE	c;
  POINTS_TYPE	d;
  
  f = cam_matrix_get_value(e, 0, 2);
  fprime = cam_matrix_get_value(eprime, 0, 2);
  a = cam_matrix_get_value(F, 1, 1);
  b = cam_matrix_get_value(F, 2, 1);
  c = cam_matrix_get_value(F, 1, 2);
  d = cam_matrix_get_value(F, 2, 2);
  cam_matrix_set_value(F, 0, 0, f * fprime * d);
  cam_matrix_set_value(F, 0, 1, -f * b);
  cam_matrix_set_value(F, 0, 2, -f * d);
  cam_matrix_set_value(F, 1, 0, -fprime * c);
  cam_matrix_set_value(F, 2, 0, -fprime * d);
}
/* end step 6 */

/* begin step 7 */
void				cam_solve_polynom(CamMatrix *F)
{
  double			coeffs[7];
  double			z[12];
  POINTS_TYPE			a;
  POINTS_TYPE			b;
  POINTS_TYPE			c;
  POINTS_TYPE			d;
  POINTS_TYPE			f;
  POINTS_TYPE			fprime;
  gsl_poly_complex_workspace	*w;

  w = gsl_poly_complex_workspace_alloc(7);

  a = cam_matrix_get_value(F, 1, 1);
  b = cam_matrix_get_value(F, 2, 1);
  c = cam_matrix_get_value(F, 1, 2);
  d = cam_matrix_get_value(F, 2, 2);
  f = -cam_matrix_get_value(F, 0, 1) / b;
  fprime = -cam_matrix_get_value(F, 1, 0) / c;

  /* begin integrity check */
  if (-cam_matrix_get_value(F, 2, 0) / d != fprime || -cam_matrix_get_value(F, 0, 2) / d != f || cam_matrix_get_value(F, 0, 0) / d != f * fprime)
    {
      printf("cam_solve_polynom : failed the integrity check\n");
      exit (-1);
    }
  /* end integrity check */
  
  coeffs[6] = -pow(f, 4) * pow(a, 2) * c * d + pow(f, 4) * a * b * pow(c, 2);
  coeffs[5] = -pow(f, 4) * pow(a, 2) * pow(d, 2) + pow(f, 4)* pow(b, 2) * pow(c, 2) + pow(a, 4) +
    2 * pow(fprime, 2) * pow(a, 2) * pow(c, 2) + pow(fprime, 4) * pow(c, 4);;
  coeffs[4] = - pow(f, 4) * pow(a, 2) * pow(d, 2) + pow(f, 4) * pow(b , 2) * c * d - 2 * pow(f, 2) * pow(a, 2) * c * d +
    2 * pow(f, 2) * a * b * pow(c, 2) + 4 * pow(fprime, 2) * a *b * pow(c, 2) + 4 * pow(a, 3) * b + 4 * pow(fprime, 2) * pow(a, 2) * c * d +
    4 * pow(fprime, 4) * pow(c, 3) * d;
  coeffs[3] = -2 * pow(f, 2) * pow(a, 2) * pow(d, 2) + 2 * pow(f, 2) * pow(b, 2) * pow(c, 2) + 8 * pow(fprime, 2) * a * b *c * d +
    2 * pow(fprime, 2) * pow(b ,2) * pow(c, 2) + 6 * pow(a, 2) * pow(b, 2) + 2 * pow(fprime, 2) * pow(a, 2) * pow(d, 2) +
    6 * pow(fprime, 4) * pow(c, 2) * pow(d, 2);;
  coeffs[2] = -2 * pow(f, 2) * a * b * pow(d, 2) + 2 * pow(f, 2) * pow(b, 2) * c * d + 4 * a * pow(b, 3) + 4* pow(fprime, 2) * a * b * pow(d, 2) +
    4 * pow(fprime, 2) * pow(b , 2) * c * d - pow(a, 2) * c * d + pow(a, 2) * b * d + a * b * pow(c, 2) + 4 * pow(fprime, 4) * c * pow(d, 3);
  coeffs[1] = pow(b, 4) - pow(a, 2) * pow(d, 2) + pow(b, 2) * pow(c, 2) + 2 * pow(fprime, 2) * pow(b, 2) * pow(d, 2) + pow(fprime, 4) * pow(d, 4);
  coeffs[0] = -a * b * pow(d, 2) + pow(b, 2) * c * d;

  gsl_poly_complex_solve(coeffs, 7, w, z);
  gsl_poly_complex_workspace_free(w);
}
