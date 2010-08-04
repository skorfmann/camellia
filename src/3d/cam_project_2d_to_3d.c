
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
#include "cam_project_2d_to_3d.h"

/*#define	DEBUG*/

/*******************************************/
/* solve :                                 */
/*             [r1 r2 r3 t1]               */
/* [ x y 1 ] = [r4 r5 r6 t2] * [ X Y Z 1 ] */
/*             [r7 r8 r9 t3]               */
/*******************************************/

void		cam_compute_vector_to_3d_point(CamMatrix *v, CamMatrix *t, CamMatrix *Rt, Cam2dPoint *pt)
{
  POINTS_TYPE	x;
  POINTS_TYPE	y;
  POINTS_TYPE	X;
  POINTS_TYPE	Y;
  POINTS_TYPE	Z;
  POINTS_TYPE	r1, r1c;
  POINTS_TYPE	r2, r2c;
  POINTS_TYPE	r3, r3c;
  POINTS_TYPE	r4, r4c;
  POINTS_TYPE	r5, r5c;
  POINTS_TYPE	r6, r6c;
  POINTS_TYPE	r7, r7c;
  POINTS_TYPE	r8, r8c;
  POINTS_TYPE	r9, r9c;
  POINTS_TYPE	t1, t1c;
  POINTS_TYPE	t2, t2c;
  POINTS_TYPE	t3, t3c;
  POINTS_TYPE	a;
  POINTS_TYPE	b;
  POINTS_TYPE	tmp1;
  POINTS_TYPE	tmp2;
#ifdef DEBUG
  CamMatrix	test;
  FILE		*debug_file;
#endif
  
  x = pt->x;
  y = pt->y;
  r1c = cam_matrix_get_value(Rt, 0, 0);
  r2c = cam_matrix_get_value(Rt, 1, 0);
  r3c = cam_matrix_get_value(Rt, 2, 0);
  r4c = cam_matrix_get_value(Rt, 0, 1);
  r5c = cam_matrix_get_value(Rt, 1, 1);
  r6c = cam_matrix_get_value(Rt, 2, 1);
  r7c = cam_matrix_get_value(Rt, 0, 2);
  r8c = cam_matrix_get_value(Rt, 1, 2);
  r9c = cam_matrix_get_value(Rt, 2, 2);
  t1c = cam_matrix_get_value(Rt, 3, 0);
  t2c = cam_matrix_get_value(Rt, 3, 1);
  t3c = cam_matrix_get_value(Rt, 3, 2);

  if (ABSF(r5c * r1c - r2c * r4c) >= SUP0)
    {
      /* 1-1 2-2 3-3 */
      if (ABSF(r1c) >= SUP0)
	{
#ifdef DEBUG
	  printf("here1\n");
#endif
	  t1 = t1c;
	  t2 = t2c;
	  t3 = t3c;
	  r1 = r1c;
	  r2 = r2c;
	  r3 = r3c;
	  r4 = r4c;
	  r5 = r5c;
	  r6 = r6c;
	  r7 = r7c;
	  r8 = r8c;
	  r9 = r9c;
	}
      /* 1-2 2-1 3-3 */
      else if (ABSF(r4c) >= SUP0)
	{
#ifdef DEBUG
	  printf("here2\n");
#endif
	  t1 = t2c;
	  t2 = t1c;
	  t3 = t3c;
	  r1 = r4c;
	  r2 = r5c;
	  r3 = r6c;
	  r4 = r1c;
	  r5 = r2c;
	  r6 = r3c;
	  r7 = r7c;
	  r8 = r8c;
	  r9 = r9c;
	}
      tmp1 = 1.0f;
    }
  else if (ABSF(r8c * r1c - r2c * r7c) >= SUP0)
    {
      /* 1-1 2-3 3-2 */
      if (ABSF(r1c) >= SUP0)
	{
#ifdef DEBUG
	  printf("here3\n");
#endif
	  t1 = t1c;
	  t2 = t3c;
	  t3 = t2c;
	  r1 = r1c;
	  r2 = r2c;
	  r3 = r3c;
	  r4 = r7c;
	  r5 = r8c;
	  r6 = r9c;
	  r7 = r4c;
	  r8 = r5c;
	  r9 = r6c;
	}
      /* 1-3 2-1 3-2 */
      else if (ABSF(r7c) >= SUP0)
	{
#ifdef DEBUG
	  printf("here4\n");
#endif
	  t1 = t3c;
	  t2 = t1c;
	  t3 = t2c;
	  r1 = r7c;
	  r2 = r8c;
	  r3 = r9c;
	  r4 = r1c;
	  r5 = r2c;
	  r6 = r3c;
	  r7 = r4c;
	  r8 = r5c;
	  r9 = r6c;
	}
      tmp1 = y;
    }
  else if (ABSF(r4c * r8c - r7c * r5c) >= SUP0)
    {
      /* 1-3 2-2 3-1 */
      if (ABSF(r7c) >= SUP0)
	{
#ifdef DEBUG
	  printf("here5\n");
#endif
	  t1 = t3c;
	  t2 = t2c;
	  t3 = t1c;
	  r1 = r7c;
	  r2 = r8c;
	  r3 = r9c;
	  r4 = r4c;
	  r5 = r5c;
	  r6 = r6c;
	  r7 = r1c;
	  r8 = r2c;
	  r9 = r3c;
	}
      /* 1-2 2-3 3-1 */
      else if (ABSF(r4c) >= SUP0)
	{
#ifdef DEBUG
	  printf("here6\n");
#endif
	  t1 = t2c;
	  t2 = t3c;
	  t3 = t1c;
	  r1 = r4c;
	  r2 = r5c;
	  r3 = r6c;
	  r4 = r7c;
	  r5 = r8c;
	  r6 = r9c;
	  r7 = r1c;
	  r8 = r2c;
	  r9 = r3c;
	}
      tmp1 = x;
    }
  else
    {
      printf("cam_compute_vector_to_3d_point : unable to determine Z\n");
      exit(-1);
    }
  
  a = r7 / r1 * (x - r2 * ( (r1 * (y - t2) - r4 * (x - t1)) / (r5 * r1 - r2 * r4)) - t1 ) +
    r8 * ( (r1 * (y - t2) - r4 * (x - t1)) / (r5 * r1 - r2 * r4) ) + t3;
  b = r9 + (r7 * r2 * r6 - r8 * r6 * r1 + r8 * r4 *r3) / (r1 * r5 - r2 * r4) - r3 * r7 / r1 - (r7 * r2 * r4 * r3) / (r1 * (r1 * r5 - r2 * r4));
  Z = (tmp1 - a) / b;
  
  /* 1&2 */
  if (ABSF(r5c * r1c - r2c * r4c) >= SUP0)
    {
      tmp1 = x;
      tmp2 = y;
      t1 = t1c;
      t2 = t2c;
      r1 = r1c;
      r2 = r2c;
      r3 = r3c;
      r4 = r4c;
      r5 = r5c;
      r6 = r6c;
    }
  /* 2&3 */
  else if (ABSF(r8c * r4c - r5c * r7c) >= SUP0)
    {
      tmp1 = y;
      tmp2 = 1.0f;
      t1 = t2c;
      t2 = t3c;
      r1 = r4c;
      r2 = r5c;
      r3 = r6c;
      r4 = r7c;
      r5 = r8c;
      r6 = r9c;
    }
  /* 1&3 */
  else if (ABSF(r8c * r1c - r2c * r7c) >= SUP0)
    {
      tmp1 = x;
      tmp2 = 1.0f;
      t1 = t1c;
      t2 = t3c;
      r1 = r1c;
      r2 = r2c;
      r3 = r3c;
      r4 = r7c;
      r5 = r8c;
      r6 = r9c;
    }
  else
    {
      printf("cam_compute_vector_to_3d_point : unable to determine Y\n");
      exit(-1);
    }

  Y = (r1 * (tmp2 - r6 * Z - t2) - r4 * (tmp1 - r3 * Z - t1) ) / (r5 * r1 - r2 * r4);
  
  /* 1 */
  if (ABSF(r1c) >= SUP0)
    {
      tmp1 = x;
      r1 = r1c;
      r2 = r2c;
      r3 = r3c;
      t1 = t1c;
      X = (x - r2c * Y - r3c * Z - t1c) / r1c;
    }
  /* 2 */
  else if (ABSF(r4c) >= SUP0)
    {
      tmp1 = y;
      r1 = r4c;
      r2 = r5c;
      r3 = r6c;
      t1 = t2c;
      X = (y - r5c * Y - r6c * Z - t2c) / r4c;
    }
  /* 3 */
  else if (ABSF(r7c) >= SUP0)
    {
      tmp1 = 1.0f;
      r1 = r7c;
      r2 = r8c;
      r3 = r9c;
      t1 = t3c;
      X = (1 - r8c * Y - r9c * Z - t3c) / r7c;
    }
  else
    {
      printf("cam_compute_vector_to_3d_point : unable to determine X\n");
      exit(-1);
    }
 
  X = (tmp1 - r2 * Y - r3 * Z - t1) / r1;

  /* check if the 3d point intersect the plan at the point position */
#ifdef DEBUG
  cam_allocate_matrix(&test, 1, 3);
  cam_matrix_set_value(v, 0, 0, X);
  cam_matrix_set_value(v, 0, 1, Y);
  cam_matrix_set_value(v, 0, 2, Z);
  cam_matrix_set_value(v, 0, 3, 1);
  cam_matrix_multiply(&test, Rt, v);

  cam_print_matrix(Rt, "Rt");
  cam_print_matrix(&test, "test");
  printf("Should be on : %f %f\n", pt->x, pt->y);

  cam_disallocate_matrix(&test);
  printf("camera center %f %f %f\n", cam_matrix_get_value(t, 0, 0),
	 cam_matrix_get_value(t, 0, 1),
	 cam_matrix_get_value(t, 0, 2));
  printf("X : %f Y : %f Z : %f\n", X, Y, Z);

  debug_file = fopen("debug_vectors", "a");
  fwrite(t->data, sizeof(POINTS_TYPE), 3, debug_file);
  fwrite(&X, sizeof(POINTS_TYPE), 1, debug_file);
  fwrite(&Y, sizeof(POINTS_TYPE), 1, debug_file);
  fwrite(&Z, sizeof(POINTS_TYPE), 1, debug_file);
  fclose(debug_file);

 #endif
  /* end check */

  cam_matrix_set_value(v, 0, 0, cam_matrix_get_value(t, 0, 0) - X);
  cam_matrix_set_value(v, 0, 1, cam_matrix_get_value(t, 0, 1) - Y);
  cam_matrix_set_value(v, 0, 2, cam_matrix_get_value(t, 0, 2) - Z);
  cam_matrix_set_value(v, 0, 3, 1);
}

/********************************************************************************/
/* solving :									*/
/* [t1x t1y t1z] + alpha * [v1x v1y v1z] = [t2x t2y t2z] + beta * [v2x v2y v2z] */
/********************************************************************************/

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
  POINTS_TYPE	alpha;
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

  if (ABSF(v2y * v1x - v1y * v2x) >= SUP0)
    {
      beta = (v1x * (t1y - t2y) + v1y * (t2x - t1x)) / (v2y * v1x - v1y * v2x);
    }
  else if (ABSF(v2z * v1x - v1z * v2x) >= SUP0)
    {
      beta = (v1x * (t1z - t2z) + v1z * (t2x - t1x)) / (v2z * v1x - v1z * v2x);
    }
  else if (ABSF(v2y * v1z - v1y * v2z) >= SUP0)
    {
      beta = (v1z * (t1y - t2y) + v1y * (t2z - t1z)) / (v2y * v1z - v1y * v2z);
    }
  else
    {
      printf("Unable to compute Beta\n");
      return (NULL);
    }

  if (ABSF(v1x) >= SUP0)
    alpha = (t2x + beta * v2x - t1x) / v1x;
  else if (ABSF(v1y) >= SUP0)
    alpha = (t2y + beta * v2y -t1y) / v1y;
  else if (ABSF(v1z) >= SUP0)
    alpha = (t2z + beta * v2z -t1z) / v1z;
  else
    {
      printf("Unable to compute Alpha\n");
      return (NULL);
    }

#ifdef DEBUG
  printf("alpha : %f beta : %f\n", alpha, beta);
  printf("t1x : %f t1y : %f t1z : %f // t2x : %f t2y : %f t2z : %f\n", t1x, t1y, t1z, t2x, t2y, t2z);
  printf("v1x : %f v1y : %f v1z : %f // v2x : %f v2y : %f v2z : %f\n", v1x, v1y, v1z, v2x, v2y, v2z);
#endif
  if (ABSF((t1z + alpha * v1z) - (t2z + beta * v2z)) >= SUP0)
    {
      printf("Error on vectors intersection : %f\n",   (t1z + alpha * v1z) - (t2z + beta * v2z));
      exit (0);
      return (NULL);
    }
  X = t2x + beta * v2x;
  Y = t2y + beta * v2y;
  Z = t2z + beta * v2z;

  res = (Cam3dPoint *)malloc(sizeof(Cam3dPoint));

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
