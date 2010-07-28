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

#include	<string.h>
#include	<stdlib.h>
#include	"cam_2d_points.h"
#include	"cam_matrix.h"
#include	"cam_list.h"
#include	"cam_write_points_to_pgm.h"
#include	"cam_3d_points_loaders.h"
#include	"cam_project_3d_to_2d.h"
#include	"cam_project_2d_to_3d.h"
#include	"cam_projection_matrix.h"
#include	"cam_p_from_f.h"
#include	"misc.h"

void		main_project_and_write_points()
{
  CamMatrix	K;
  CamMatrix	*R;
  CamMatrix	t;
  CamMatrix	P;
  CamList	*points;
  POINTS_TYPE	Kdata[9] = {1,0,0,0,1,0,0,0,1};
  POINTS_TYPE	Tdata[3] = {1.0f,0.0f,0.0f};
  CamList	*res;

  cam_allocate_matrix(&K, 3, 3);
  cam_allocate_matrix(&P, 4, 3);
  R = compute_rotation_matrix(0.0f, 0.0f, 0.0f);
  cam_allocate_matrix(&t, 1, 3);
  memcpy(K.data, Kdata, 9 * sizeof(POINTS_TYPE));
  memcpy(t.data, Tdata, 3 * sizeof(POINTS_TYPE));
  points = loadPoints1("/home/splin/manny");
  cam_compute_projection_matrix(&P, &K, R, &t);
  res = cam_project_3d_to_2d(points, &P);
  cam_center_2d_points(res, 800, 600);
  cam_write_points_to_pgm("pts.pgm", res, 800, 600,
			  255, 255, 255,
			  0, 0, 0);
  cam_disallocate_linked_list(res);
  cam_disallocate_linked_list(points);
  cam_disallocate_matrix(&t);
  cam_disallocate_matrix(R);

  cam_disallocate_matrix(&K);
  cam_disallocate_matrix(&P);
  free(R);
}

void			main_triangulate_2d_points()
{
  CamProjectionsPair	projectionPair;
  Cam2dPoint		pt;
  CamMatrix		K;
  POINTS_TYPE		Kdata[9] = {1,0,0,0,1,0,0,0,1};
  POINTS_TYPE		Tdata[3] = {0.0f,0.0f,0.0f};  
  CamMatrix		*R;
  CamMatrix		t;
  POINTS_TYPE	rval[6] = {-PI/4, PI/2, 0, PI/2, PI/4};

  int a,b,c;

  cam_allocate_matrix(&projectionPair.p1, 4, 3);
  cam_allocate_matrix(&projectionPair.p2, 4, 3);
  cam_allocate_matrix(&K, 3, 3);
  cam_allocate_matrix(&t, 1, 3);
  memcpy(K.data, Kdata, 9 * sizeof(POINTS_TYPE));
  memcpy(t.data, Tdata, 3 * sizeof(POINTS_TYPE));
  pt.x = 0.0f;
  pt.y = 0.0f;
  for (a = 0; a < 6 ; ++a)
    for (b = 0; b < 6 ; ++b)
      for (c = 0; c < 6 ; ++c)
    {
      R = compute_rotation_matrix(rval[a], rval[b], rval[c]);
      cam_compute_projection_matrix(&projectionPair.p1, &K, R, &t);
      cam_compute_projection_matrix(&projectionPair.p2, &K, R, &t);
      cam_triangulate_one_3d_point(&projectionPair, &t, &t, &K, &pt, &pt);
    }
  cam_disallocate_projections_pair(&projectionPair);
}

int		main()
{
  /*main_project_and_write_points();*/
  main_triangulate_2d_points();
 
  return (0);
}