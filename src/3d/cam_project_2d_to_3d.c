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
#include "cam_project_2d_to_3d.h"

#define	ABSF(x) (x >= 0.0f ? (x) : -(x))

/* solve : [ x y 1 ] = Rt * [ X Y Z 1 ] */
void		cam_compute_vector_to_3d_point(CamMatrix *v, CamMatrix *t, CamMatrix *Rt, Cam2dPoint *pt)
{
  POINTS_TYPE	x;
  POINTS_TYPE	y;
  POINTS_TYPE	X;
  POINTS_TYPE	Y;
  POINTS_TYPE	Z;
  POINTS_TYPE	r1;
  POINTS_TYPE	r2;
  POINTS_TYPE	r3;
  POINTS_TYPE	r4;
  POINTS_TYPE	r5;
  POINTS_TYPE	r6;
  POINTS_TYPE	r7;
  POINTS_TYPE	r8;
  POINTS_TYPE	r9;
  POINTS_TYPE	t1;
  POINTS_TYPE	t2;
  POINTS_TYPE	t3;
  POINTS_TYPE	a;
  POINTS_TYPE	b;
  CamMatrix	test;
  
  x = pt->x;
  y = pt->y;
  r1 = cam_matrix_get_value(Rt, 0, 0);
  r2 = cam_matrix_get_value(Rt, 1, 0);
  r3 = cam_matrix_get_value(Rt, 2, 0);
  r4 = cam_matrix_get_value(Rt, 0, 1);
  r5 = cam_matrix_get_value(Rt, 1, 1);
  r6 = cam_matrix_get_value(Rt, 2, 1);
  r7 = cam_matrix_get_value(Rt, 0, 2);
  r8 = cam_matrix_get_value(Rt, 1, 2);
  r9 = cam_matrix_get_value(Rt, 2, 2);
  t1 = cam_matrix_get_value(Rt, 3, 0);
  t2 = cam_matrix_get_value(Rt, 3, 1);
  t3 = cam_matrix_get_value(Rt, 3, 2);

  if (ABSF(r5 * r1 - r2 * r4) >= 0.001f)
    {
      if (ABSF(r1) >= 0.001f)
	{
	  a = r7 / r1 * (x - r2 * ( (r1 * (y - t2) - r4 * (x - t1)) / (r5 * r1 - r2 * r4)) - t1 ) +
	    r8 * ( (r1 * (y - t2) - r4 * (x - t1)) / (r5 * r1 - r2 * r4) ) + t3;
	  b = r9 + (r7 * r2 * r6 - r8 * r6 * r1 + r8 * r4 *r3) / (r1 * r5 - r2 * r4) - r3 * r7 / r1 - (r7 * r2 * r4 * r3) / (r1 * (r1 * r5 - r2 * r4));
	}
      else if (ABSF(r4) >= 0.001f)
	{
	  a = r7 / r4 * (y - r5 * ( (r4 * (x - t1) - r1 * (y - t2)) / (r2 * r4 - r5 * r1)) - t2 ) +
	    r8 * ( (r4 * (x - t1) - r1 * (y - t2)) / (r2 * r4 - r5 * r1) ) + t3;
	  b = r9 + (r7 * r5 * r3 - r8 * r3 * r4 + r8 * r1 * r6) / (r4 * r2 - r5 * r1) - r6 * r7 / r4 - (r7 * r5 * r1 * r6) / (r4 * (r4 * r2 - r5 * r1));
	}
      if (ABSF(b) >= 0.001f)
	Z = (1.0f - a) / b;
      else
	printf("Erf\n");
    }
  else if (ABSF(r8 * r1 - r2 * r7) >= 0.001f)
    {
      if (ABSF(r1) >= 0.001f)
	{
	  a = r4 / r1 * (x - r2 * ( (r1 * (1 - t3) - r7 * (x - t1)) / (r8 * r1 - r2 * r7)) - t1 ) +
	    r5 * ( (r1 * (1 - t3) - r7 * (x - t1)) / (r8 * r1 - r2 * r7) ) + t2;
	  b = r6 + (r4 * r2 * r9 - r5 * r9 * r1 + r5 * r7 * r3) / (r1 * r8 - r2 * r7) - r3 * r4 / r1 - (r4 * r2 * r7 * r3) / (r1 * (r1 * r8 - r2 * r7));
	  Z = (y - a) / b;
	}
      else if (ABSF(r7) >= 0.001f)
	{
	  a = r7 / r1 * (x - r2 * ( (r1 * (y - t2) - r4 * (x - t1)) / (r5 * r1 - r2 * r4)) - t1 ) +
	    r8 * ( (r1 * (y - t2) - r4 * (x - t1)) / (r5 * r1 - r2 * r4) ) + t3;
	  b = r9 + (r7 * r2 * r6 - r8 * r6 * r1 + r8 * r4 *r3) / (r1 * r5 - r2 * r4) - r3 * r7 / r1 - (r7 * r2 * r4 * r3) / (r1 * (r1 * r5 - r2 * r4));
	  /*  printf("here4\n");*/
	  a = r4 / r7 * (1 - r8 * ( (r7 * (x - t1) - r1 * (1 - t3)) / (r2 * r7 - r8 * r1)) - t3 ) +
	    r5 * ( (r7 * (x - t1) - r1 * (1 - t3)) / (r2 * r7 - r8 * r1) ) + t2;
	  b = r6 + (r4 * r8 * r3 - r5 * r3 * r7 + r5 * r1 *r9) / (r7 * r2 - r8 * r1) - r9 * r4 / r7 - (r4 * r8 * r1 * r9) / (r7 * (r7 * r2 - r8 * r1));
	}
      if (ABSF(b) >= 0.001f)
	Z = (y - a) / b;
      else
	printf("Erf\n");
    }
  else if (ABSF(r4 * r8 - r7 * r5) >= 0.001f)
    {
      if (ABSF(r7) >= 0.001f)
	{
	  a = r7 / r1 * (x - r2 * ( (r1 * (y - t2) - r4 * (x - t1)) / (r5 * r1 - r2 * r4)) - t1 ) +
	    r8 * ( (r1 * (y - t2) - r4 * (x - t1)) / (r5 * r1 - r2 * r4) ) + t3;
	  b = r9 + (r7 * r2 * r6 - r8 * r6 * r1 + r8 * r4 *r3) / (r1 * r5 - r2 * r4) - r3 * r7 / r1 - (r7 * r2 * r4 * r3) / (r1 * (r1 * r5 - r2 * r4));
	  /*	  printf("here5\n");*/
	  a = r1 / r7 * (1 - r8 * ( (r7 * (y - t2) - r4 * (1 - t3)) / (r5 * r7 - r8 * r4)) - t3 ) +
	    r2 * ( (r7 * (y - t2) - r4 * (1 - t3)) / (r5 * r7 - r8 * r4) ) + t1;
	  b = r3 + (r1 * r8 * r6 - r2 * r6 * r7 + r2 * r4 *r9) / (r7 * r5 - r8 * r4) - r9 * r1 / r7 - (r1 * r8 * r4 * r9) / (r7 * (r7 * r5 - r8 * r4));
	}
      else if (ABSF(r4) >= 0.001f)
	{
	  a = r1 / r4 * (y - r5 * ( (r4 * (1 - t3) - r7 * (y - t2)) / (r8 * r4 - r5 * r7)) - t2 ) +
	    r2 * ( (r4 * (1 - t3) - r7 * (y - t2)) / (r8 * r4 - r5 * r7) ) + t1;
	  b = r3 + (r1 * r5 * r9 - r2 * r9 * r4 + r2 * r7 * r6) / (r4 * r8 - r5 * r7) - r6 * r1 / r4 - (r1 * r5 * r7 * r6) / (r4 * (r4 * r8 - r5 * r7));
	}
      if (ABSF(b) >= 0.001f)
	Z = (x - a) / b;
      else
	printf("Erf\n");
    }
  else
    {
      printf("cam_compute_vector_to_3d_point : unable to determine Z\n");
      exit(-1);
    }


  if (ABSF(r5 * r1 - r2 * r4) >= 0.001f)
    Y = (r1 * (y - r6 * Z - t2) - r4 * (x - r3 * Z - t1) ) / (r5 * r1 - r2 * r4);
  else if (ABSF(r8 * r4 - r5 * r7) >= 0.001f)
    Y = (r4 * (1 - r9 * Z - t3) - r7 * (y - r6 * Z - t2) ) / (r8 * r4 - r5 * r7);
  else if (ABSF(r8 * r1 - r2 * r7) >= 0.001f)
    Y = (r1 * (1 - r9 * Z - t3) - r7 * (x - r3 * Z - t1) ) / (r8 * r1 - r2 * r7);
  else
    {
      printf("cam_compute_vector_to_3d_point : unable to determine Y\n");
      exit(-1);
    }

  if (ABSF(r1) >= 0.001f)
    X = (x - r2 * Y - r3 * Z - t1) / r1;
  else if (ABSF(r4) >= 0.001f)
    X = (y - r5 * Y - r6 * Z - t2) / r4;
  else if (ABSF(r7) >= 0.001f)
    X = (1 - r8 * Y - r9 * Z - t2) / r7;
  else
    {
      printf("cam_compute_vector_to_3d_point : unable to determine X\n");
      exit(-1);
    }
 
  /* check */
  /*
  cam_allocate_matrix(&test, 1, 3);
  cam_matrix_set_value(v, 0, 0, X);
  cam_matrix_set_value(v, 0, 1, Y);
  cam_matrix_set_value(v, 0, 2, Z);
  cam_matrix_set_value(v, 0, 3, 1);
  cam_matrix_multiply(&test, Rt, v);
  if (ABSF(cam_matrix_get_value(&test, 0, 2)) >= 1.001 ||
      (ABSF(cam_matrix_get_value(&test, 0, 0)) >= 4.001f && ABSF(cam_matrix_get_value(&test, 0, 0)) <= 3.999f) ||
      (ABSF(cam_matrix_get_value(&test, 0, 1)) >= 12.001f && ABSF(cam_matrix_get_value(&test, 0, 1)) <= 11.999f))
      cam_print_matrix(&test, "test");
  */
  /* end check */

  cam_matrix_set_value(v, 0, 0, cam_matrix_get_value(t, 0, 0) - X);
  cam_matrix_set_value(v, 0, 1, cam_matrix_get_value(t, 0, 1) - Y);
  cam_matrix_set_value(v, 0, 2, cam_matrix_get_value(t, 0, 2) - Z);
  cam_matrix_set_value(v, 0, 3, 1);
}

/* TODO */
Cam3dPoint	*cam_vectors_intersection(CamMatrix *v1, CamMatrix *t1, CamMatrix *v2, CamMatrix *t2)
{
  Cam3dPoint	*res;
  POINTS_TYPE	X;
  POINTS_TYPE	Y;
  POINTS_TYPE	Z;
  POINTS_TYPE	v1x;
  POINTS_TYPE	v1y;
  POINTS_TYPE	v1z;
  POINTS_TYPE	v2x;
  POINTS_TYPE	v2y;
  POINTS_TYPE	v2z;
  POINTS_TYPE	t1x;
  POINTS_TYPE	t1y;
  POINTS_TYPE	t1z;
  POINTS_TYPE	t2x;
  POINTS_TYPE	t2y;
  POINTS_TYPE	t2z;
  POINTS_TYPE	beta;
  
  v1x = cam_matrix_get_value(v1, 0, 0);
  v1y = cam_matrix_get_value(v1, 0, 1);
  v1z = cam_matrix_get_value(v1, 0, 2);
  v2x = cam_matrix_get_value(v2, 0, 0);
  v2y = cam_matrix_get_value(v2, 0, 1);
  v2z = cam_matrix_get_value(v2, 0, 2);
  t1x = cam_matrix_get_value(t1, 0, 0);
  t1y = cam_matrix_get_value(t1, 0, 1);
  t1z = cam_matrix_get_value(t1, 0, 2);
  t2x = cam_matrix_get_value(t2, 0, 0);
  t2y = cam_matrix_get_value(t2, 0, 1);
  t2z = cam_matrix_get_value(t2, 0, 2);

  res = (Cam3dPoint *)malloc(sizeof(Cam3dPoint));

  beta = (t1y - t2y + v1y * (t2x - t1x) / v1x) / (v2y + v2x * v1y / v1x);
  
  X = t2x + beta * v2x;
  Y = t2y + beta * v2y;
  Z = t2y + beta * v2y;

  res->x = X;
  res->y = Y;
  res->z = Z;
  res->dist = 1;

  return (res);
}

Cam3dPoint	*cam_triangulate_one_3d_point(CamProjectionsPair *projectionPair, CamMatrix *t1, CamMatrix *t2, CamMatrix *K, Cam2dPoint *a, Cam2dPoint *b)
{
  CamMatrix	v1;
  CamMatrix	v2;
  Cam3dPoint	*res;

  cam_allocate_matrix(&v1, 1 ,4);
  cam_allocate_matrix(&v2, 1 ,4);

  /* TODO : compute K-1 . P pour recuperer Rt */
  cam_compute_vector_to_3d_point(&v1, t1, &projectionPair->p1, a);
  cam_compute_vector_to_3d_point(&v2, t2, &projectionPair->p2, b);
  res = cam_vectors_intersection(&v1, t1, &v2, t2);
  cam_disallocate_matrix(&v1);
  cam_disallocate_matrix(&v2);
  return (res);
}
