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
#include	<stdio.h>
#include	"cam_2d_points.h"
#include	"cam_matrix.h"
#include	"cam_list.h"
#include	"cam_points_to_ppm.h"
#include	"cam_3d_points_loaders.h"
#include	"cam_project_3d_to_2d.h"
#include	"cam_project_2d_to_3d.h"
#include	"cam_projection_matrix.h"
#include	"cam_p_from_f.h"
#include	"misc.h"

/*#define		DEBUG*/

void		main_project_and_write_points()
{
  CamMatrix	K;
  CamMatrix	*R;
  CamMatrix	t;
  CamMatrix	P;
  CamList	*points;
  POINTS_TYPE	Kdata[9] = {1,0,0,0,1,0,0,0,1};
  POINTS_TYPE	Tdata[3] = {0.0f,0.0f,0.0f};
  CamList	*res;

  cam_allocate_matrix(&K, 3, 3);
  cam_allocate_matrix(&P, 4, 3);
  R = compute_rotation_matrix(0.0f, 0.0f, 0.0f);
  cam_allocate_matrix(&t, 1, 3);
  memcpy(K.data, Kdata, 9 * sizeof(POINTS_TYPE));
  memcpy(t.data, Tdata, 3 * sizeof(POINTS_TYPE));

  /*points = loadPoints1("/home/splin/manny");*/
  points = NULL;
  points = cam_add_to_linked_list(points, (Cam3dPoint *)malloc(sizeof(Cam2dPoint)));
  ((Cam3dPoint *)(points->data))->x = 2.0f;
  ((Cam3dPoint *)(points->data))->y = 1.0f;
  ((Cam3dPoint *)(points->data))->z = 0.0f;
  
  cam_compute_projection_matrix(&P, &K, R, &t);
  res = cam_project_3d_to_2d(points, &P);
  
  printf("Result point : %f %f\n", ((Cam2dPoint *)(res->data))->x, ((Cam2dPoint *)(res->data))->y);

  cam_center_2d_points(res, 800, 600);
  cam_points_to_ppm("pts.ppm", res, 800, 600,
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

void			main_triangulate_2d_points_tests()
{
  CamProjectionsPair	projectionPair;
  Cam2dPoint		pt;
  CamMatrix		K;
  POINTS_TYPE		Kdata[9] = {1,0,0,0,1,0,0,0,1};
  POINTS_TYPE		Tdata[3] = {8.0f,0.1f,0.5f};  
  CamMatrix		*R;
  CamMatrix		t;

  float a,b,c;

  cam_allocate_matrix(&projectionPair.p1, 4, 3);
  cam_allocate_matrix(&projectionPair.p2, 4, 3);
  cam_allocate_matrix(&K, 3, 3);
  cam_allocate_matrix(&t, 1, 3);
  memcpy(K.data, Kdata, 9 * sizeof(POINTS_TYPE));
  memcpy(t.data, Tdata, 3 * sizeof(POINTS_TYPE));
  pt.x = -4.0f;
  pt.y = 12.0f;

  for (a = -PI; a < PI ; a+=PI/12)
    for (b = -PI; b < PI ; b+=PI/12)
      for (c = -PI; c < PI ; c+=PI/12)
	{
	  R = compute_rotation_matrix(a, b, c);
	  cam_compute_projection_matrix(&projectionPair.p1, &K, R, &t);
	  cam_compute_projection_matrix(&projectionPair.p2, &K, R, &t);
	  cam_triangulate_one_3d_point(&projectionPair, &t, &t, &K, &pt, &pt);
	}
    
  cam_disallocate_projections_pair(&projectionPair);
}

void			main_triangulate_2d_points()
{
  CamProjectionsPair	projectionPair;
  CamMatrix		K;
  POINTS_TYPE		Kdata[9] = {1.0f,0.0f,0.0f,0.0f,1.0f,0.0f,0.0f,0.0f,1.0f};
  POINTS_TYPE		Tdata1[3] = {0.0f,0.0f,0.0f};  
  POINTS_TYPE		Tdata2[3] = {1.0f,0.0f,0.0f};
  CamMatrix		*R;
  CamMatrix		t1;
  CamMatrix		t2;
  CamList		*pts1;
  CamList		*pts2;
  CamList		*ppts1;
  CamList		*ppts2;
  CamList		*points;
  Cam3dPoint		*pt3d;
  FILE			*file;
  char			*filename = "manny";
  char			*path1;
  char			*path2;
#ifdef DEBUG
  FILE			*debug_file;
  CamList		*tailPoints;
#endif

  path1 = (char*)malloc((strlen("data/3d/3dmodels/") + strlen(filename) + 1) * sizeof (char));
  path2 = (char*)malloc((strlen("data/3d/results/") + strlen(filename) + 1) * sizeof (char));
  sprintf(path1, "data/3d/3dmodels/%s", filename);
  sprintf(path2, "data/3d/results/%s", filename);
  points = loadPoints1(path1);
  file = fopen(path2, "w");
  free(path1);
  free(path2);

  cam_allocate_matrix(&projectionPair.p1, 4, 3);
  cam_allocate_matrix(&projectionPair.p2, 4, 3);
  cam_allocate_matrix(&K, 3, 3);
  cam_allocate_matrix(&t1, 1, 3);
  cam_allocate_matrix(&t2, 1, 3);
  memcpy(K.data, Kdata, 9 * sizeof(POINTS_TYPE));

  memcpy(t1.data, Tdata1, 3 * sizeof(POINTS_TYPE));
  R = compute_rotation_matrix(0.0f, 0.0f, 0.0f);
  cam_compute_projection_matrix(&projectionPair.p1, &K, R, &t1);
  pts1 = cam_project_3d_to_2d(points, &projectionPair.p1);
  cam_disallocate_matrix(R);
  free(R);

  memcpy(t2.data, Tdata2, 3 * sizeof(POINTS_TYPE));
  R = compute_rotation_matrix(0.0f, 0.0f, 0.0f);
  cam_compute_projection_matrix(&projectionPair.p2, &K, R, &t2);
  pts2 = cam_project_3d_to_2d(points, &projectionPair.p2);
  cam_disallocate_matrix(R);
  free(R);

#ifdef DEBUG
  tailPoints = points;
  while (tailPoints->next)
    tailPoints = tailPoints->next;
  printf("3D point : %f %f %f\n", ((Cam3dPoint *)tailPoints->data)->x, ((Cam3dPoint *)tailPoints->data)->y, ((Cam3dPoint *)tailPoints->data)->z);
  debug_file = fopen("debug_vectors", "w");
  fwrite(t1.data, sizeof(POINTS_TYPE), 3, debug_file);
  fwrite(tailPoints->data, sizeof(POINTS_TYPE), 3, debug_file);
  fwrite(t2.data, sizeof(POINTS_TYPE), 3, debug_file);
  fwrite(tailPoints->data, sizeof(POINTS_TYPE), 3, debug_file);
  fclose(debug_file);
#endif

  ppts1 = pts1;
  ppts2 = pts2;
  while (ppts1 && ppts2)
    {
      pt3d = cam_triangulate_one_3d_point(&projectionPair, &t1, &t2, &K, ppts1->data, ppts2->data);

      if (pt3d)
	{
	  fwrite(pt3d, sizeof(POINTS_TYPE), 3, file);
	  free(pt3d);
	}
      ppts1 = ppts1->next;
      ppts2 = ppts2->next;
    }

  path1 = (char*)malloc((strlen("data/3d/results/") + strlen(filename) + strlen("_viewpoint_X.ppm")  + 1) * sizeof (char));
  sprintf(path1, "data/3d/results/%s_viewpoint_1.ppm", filename);
  cam_center_2d_points(pts1, 800, 600);
  cam_points_to_ppm(path1, pts1, 800, 600,
			  255, 255, 255,
			  0, 0, 0);
  sprintf(path1, "data/3d/results/%s_viewpoint_2.ppm", filename);
  cam_center_2d_points(pts2, 800, 600);
  cam_points_to_ppm(path1, pts2, 800, 600,
			  255, 255, 255,
			  0, 0, 0);
  free(path1);

    
  cam_disallocate_projections_pair(&projectionPair);
  cam_disallocate_linked_list(pts1);
  cam_disallocate_linked_list(pts2);
  cam_disallocate_linked_list(points);
  cam_disallocate_matrix(&t1);
  cam_disallocate_matrix(&t2);
  cam_disallocate_matrix(&K);

  fclose(file);
}

int		main()
{
  /*main_project_and_write_points();*/
  /*main_triangulate_2d_points_tests();*/
  main_triangulate_2d_points();
  
  return (0);
}
