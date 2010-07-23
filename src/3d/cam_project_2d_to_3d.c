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
#include "cam_project_2d_to_3d.h"

/* solve : [ x y 1 ] = Rt * [ X Y Z 1 ] */
void		cam_compute_vector_to_3d_point(CamMatrix *v, CamMatrix *Rt, Cam2dPoint *pt)
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
  
  a = r7 / r1 * (x - r2 * ( (r1 * (y - t2) - r4 * (x - t1)) / (r5 * r1 - r2 * r4)) - t1 ) +
    r8 * ( (r1 * (y - t2) - r4 * (x - t1)) / (r5 * r1 - r2 * r4) ) + t3;
  b = r9 + (r2 * r7 *(r6 - r3 * r4 / r1)) / (r1 * (r5 * r1 - r2 * r4)) - r3 * r7 / r1 - r8 * ((r6 - r3 * r4 / r1) / (r5 * r1 - r2 * r4));
  Z = (1.0f - a) / b;
  Y = (r1 * (y - Z * (r6 - r3 * r4 / r1 ) - t2) - r4 * (x - t1) ) / (r5 * r1 - r2 * r4);
  X = (x - r2 * Y - r3 * Z - t1) / r1;
  cam_matrix_set_value(v, 0, 0, X);
  cam_matrix_set_value(v, 0, 1, Y);
  cam_matrix_set_value(v, 0, 2, Z);
}

Cam3dPoint	*cam_triangulate(CamProjectionsPair *projectionPair, CamMatrix *K, Cam2dPoint *a, Cam2dPoint *b)
{
  CamMatrix	v1;
  CamMatrix	v2;

  cam_allocate_matrix(&v1, 1, 4);
  cam_allocate_matrix(&v2, 1, 4);
  /* TODO : compute K-1 . P pour recuperer Rt */
  cam_compute_vector_to_3d_point(&v1, &projectionPair->p1, a);
  cam_compute_vector_to_3d_point(&v2, &projectionPair->p2, b);
  cam_disallocate_matrix(&v1);
  cam_disallocate_matrix(&v2);
  return (NULL);
}
