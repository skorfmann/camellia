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
#include	<cv.h>
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

typedef struct
{
  Cam2dPoint	pt2d;
  int		frameIndex;
} mask2d;

typedef struct
{
  Cam3dPoint	pt3d;
  CamList	*mask;
} sbaPoint;

/***********************************/
/* Begin of quaternions conversion */
/***********************************/

void		cam_write_quaternion_camera_projection(FILE *fd, CamMatrix *P)
{
  POINTS_TYPE	q[4];
  POINTS_TYPE	tmp[4];
  int		i;
  int		maxPos;
  POINTS_TYPE	mag;

  tmp[0] = 1.0f + cam_matrix_get_value(P, 0, 0) + cam_matrix_get_value(P, 1, 1) + cam_matrix_get_value(P, 2, 2);
  tmp[1] = 1.0f + cam_matrix_get_value(P, 0, 0) - cam_matrix_get_value(P, 1, 1) - cam_matrix_get_value(P, 2, 2);
  tmp[2] = 1.0f - cam_matrix_get_value(P, 0, 0) + cam_matrix_get_value(P, 1, 1) - cam_matrix_get_value(P, 2, 2);
  tmp[3] = 1.0f - cam_matrix_get_value(P, 0, 0) - cam_matrix_get_value(P, 1, 1) + cam_matrix_get_value(P, 2, 2);

  for (i = 0, mag = -1.0f ; i < 4 ; ++i)
    {
      if (mag < tmp[i])
	{
	  mag = tmp[i];
	  maxPos = i;
	}
    }

  switch (maxPos)
    {
    case 0:
      q[0] = sqrt(tmp[0]) * 0.5f;
      q[1] = (cam_matrix_get_value(P, 1, 2) - cam_matrix_get_value(P, 2, 1)) / (4.0f * q[0]);
      q[2] = (cam_matrix_get_value(P, 2, 0) - cam_matrix_get_value(P, 0, 2)) / (4.0f * q[0]);
      q[3] = (cam_matrix_get_value(P, 0, 1) - cam_matrix_get_value(P, 1, 0)) / (4.0f * q[0]);
      break;
    case 1:
      q[1] = sqrt(tmp[1]) * 0.5f;
      q[0] = (cam_matrix_get_value(P, 1, 2) - cam_matrix_get_value(P, 2, 1)) / (4.0f * q[1]);
      q[2] = (cam_matrix_get_value(P, 0, 1) + cam_matrix_get_value(P, 1, 0)) / (4.0f * q[1]);
      q[3] = (cam_matrix_get_value(P, 2, 0) + cam_matrix_get_value(P, 0, 2)) / (4.0f * q[1]);
      break;
    case 2:
      q[2] = sqrt(tmp[2]) * 0.5f;
      q[0] = (cam_matrix_get_value(P, 2, 0) - cam_matrix_get_value(P, 0, 2)) / (4.0f * q[2]);
      q[1] = (cam_matrix_get_value(P, 0, 1) + cam_matrix_get_value(P, 1, 0)) / (4.0f * q[2]);
      q[3] = (cam_matrix_get_value(P, 1, 2) + cam_matrix_get_value(P, 2, 1)) / (4.0f * q[2]);
      break;
    case 3:
      q[3]= sqrt(tmp[3]) * 0.5f;
      q[0]= (cam_matrix_get_value(P, 0, 1) - cam_matrix_get_value(P, 1, 0)) / (4.0f * q[3]);
      q[1]= (cam_matrix_get_value(P, 2, 0) + cam_matrix_get_value(P, 0, 2)) / (4.0f * q[3]);
      q[2]= (cam_matrix_get_value(P, 1, 2) + cam_matrix_get_value(P, 2, 1)) / (4.0f * q[3]);
      break;
    default:
      printf("Error in write_quaternion_camera\n");
      exit (-1);
    }

  mag = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
  mag= 1.0f / sqrt(mag);
  q[0] *= mag;
  q[1] *= mag;
  q[2] *= mag;
  q[3] *= mag;

  fprintf(fd, "%f %f %f %f %f %f %f\n", q[0], q[1], q[2], q[3],
	  cam_matrix_get_value(P, 3, 0),
	  cam_matrix_get_value(P, 3, 1),
	  cam_matrix_get_value(P, 3, 2));
  return ;
}

void		cam_write_quaternion_camera_params(FILE *fd, POINTS_TYPE *R, POINTS_TYPE *t)
{
  CamMatrix	*Rot;
  CamMatrix	P;
  int		i;
  int		j;
  
  cam_allocate_matrix(&P, 4, 3);
  Rot = compute_rotation_matrix(R[0], R[1], R[2]);
  for (j = 0; j < 3 ; ++j)
    for (i = 0; i < 3 ; ++i)
      {
	cam_matrix_set_value(&P, i, j, cam_matrix_get_value(Rot, i, j));
      }
  for (i = 0 ; i < 3 ; ++i)
    cam_matrix_set_value(&P, 3, i, t[i]);
  cam_write_quaternion_camera_projection(fd, &P);
  cam_disallocate_matrix(Rot);
  free(Rot);
  cam_disallocate_matrix(&P);
}

/*********************************/
/* End of quaternions conversion */
/*********************************/

void		write_point(FILE *fd, Cam3dPoint *pt3d, float *p1, float *p2)
{
  fprintf(fd, "%f %f %f 2 0 %f %f 1 %f %f\n", pt3d->x, pt3d->y, pt3d->z, p1[0], p1[1], p2[0], p2[1]);
}

void		cam_print_sba()
{

}

void		cam_print_sba_3d_points_text(CamList *points, char *file)
{
  FILE		*output;
  CamList	*tmp;

  output = fopen(file, "w");
  tmp = points;
  while (tmp)
    {
      fprintf(output, "%f %f %f\n", ((sbaPoint *)tmp->data)->pt3d.x, ((sbaPoint *)tmp->data)->pt3d.y, ((sbaPoint *)tmp->data)->pt3d.z);
      tmp = tmp->next;
    }
}

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
	  cam_triangulate_one_perfect_3d_point(&projectionPair, &t, &t, &K, &pt, &pt);
	}
    
  cam_disallocate_projections_pair(&projectionPair);
}

void			main_triangulate_2d_points_noisy_2_views()
{
  CamProjectionsPair	projectionPair;
  CamMatrix		K;
  POINTS_TYPE		Kdata[9] = {1.0f,0.0f,0.0f,0.0f,1.0f,0.0f,0.0f,0.0f,1.0f};
  POINTS_TYPE		Tdata1[3] = {0.0f,0.0f,0.0f};  
  POINTS_TYPE		Rdata1[3] = {0.0f,0.0f,0.0f};
  POINTS_TYPE		Tdata2[3] = {1.0f, 0.0f, 0.0f};
  POINTS_TYPE		Rdata2[3] = {PI/6, 0.0f,0.0f};
  CamMatrix		*R;
  CamMatrix		t1;
  CamMatrix		t2;
  CamMatrix		F;
  CamList		*pts1;
  CamList		*pts2;
  CamList		*ppts1;
  CamList		*ppts2;
  CamList		*points;
  Cam3dPoint		*pt3d;
  FILE			*file;
  FILE			*cameras;
  FILE			*pts;
  char			*filename = "manny";
  char			*path1;
  char			*path2;
  CvMat			*proj1;
  CvMat			*proj2;
  CvMat			*points4d;
  CvMat			*points1;
  CvMat			*points2;
  CvMat			*status;
  CvMat			*fundamental_matrix;
  CvMat			*new_points1;
  CvMat			*new_points2;
  int			fm_count;
  int			i;
  int			j;
  int			nbPoints;
  int			curNbPoints;
  Cam3dPoint		tmp3d;

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
  cam_allocate_matrix(&F, 3, 3);
  memcpy(K.data, Kdata, 9 * sizeof(POINTS_TYPE));

  memcpy(t1.data, Tdata1, 3 * sizeof(POINTS_TYPE));
  R = compute_rotation_matrix(Rdata1[0], Rdata1[1], Rdata1[2]);
  cam_compute_projection_matrix(&projectionPair.p1, &K, R, &t1);
  pts1 = cam_project_3d_to_2d(points, &projectionPair.p1);
  cam_disallocate_matrix(R);
  free(R);

  memcpy(t2.data, Tdata2, 3 * sizeof(POINTS_TYPE));
  R = compute_rotation_matrix(Rdata2[0], Rdata2[1], Rdata2[2]);
  cam_compute_projection_matrix(&projectionPair.p2, &K, R, &t2);
  pts2 = cam_project_3d_to_2d(points, &projectionPair.p2);
  cam_add_noise_to_2d_points(pts2, 500.0f, 1);
  cam_disallocate_matrix(R);
  free(R);

  ppts1 = pts1;
  ppts2 = pts2;
  points1 = cvCreateMat(1,pts1->index,CV_32FC2);
  points2 = cvCreateMat(1,pts1->index,CV_32FC2);
  status = cvCreateMat(1,pts1->index,CV_8UC1);
  new_points1 = cvCreateMat(2,pts1->index,CV_32FC1);
  new_points2 = cvCreateMat(2,pts1->index,CV_32FC1);

  for( i = 0; i < pts1->index; i++ )
    {
      points1->data.fl[i*2] = ((Cam2dPoint *)(ppts1->data))->x;
      points1->data.fl[i*2+1] = ((Cam2dPoint *)(ppts1->data))->y;
      points2->data.fl[i*2] = ((Cam2dPoint *)(ppts2->data))->x;
      points2->data.fl[i*2+1] = ((Cam2dPoint *)(ppts2->data))->y;
      ppts1 = ppts1->next;
      ppts2 = ppts2->next;
    }
  fundamental_matrix = cvCreateMat(3,3,CV_32FC1);
  fm_count = cvFindFundamentalMat( points1,points2,fundamental_matrix,
				   CV_FM_RANSAC,1.0,0.99,status );

  if (!fm_count)
    {
      printf("main_triangulate_2d_points_noisy : no fundamental matrix found\n");
      exit (-1);
    }

  proj1 = cvCreateMat(3,4,CV_32FC1);
  proj2 = cvCreateMat(3,4,CV_32FC1);

  for (j = 0 ; j < 3 ; ++j)
    for (i = 0 ; i < 3 ; ++i)
      cam_matrix_set_value(&F, i, j, cvmGet(fundamental_matrix, j, i));
  
  cam_compute_p_from_f(&F, &projectionPair);

  cameras = fopen("cameras.txt", "w");
  cam_write_quaternion_camera_projection(cameras, &projectionPair.p1);
  cam_write_quaternion_camera_projection(cameras, &projectionPair.p2);
  fclose(cameras);

  for (j = 0 ; j < 3 ; ++j)
    for (i = 0 ; i < 4 ; ++i)
      {
	cvmSet(proj1, j, i, cam_matrix_get_value(&projectionPair.p1, i, j));
	cvmSet(proj2, j, i, cam_matrix_get_value(&projectionPair.p2, i, j));
      }

  cvReleaseMat(&points1);
  cvReleaseMat(&points2);
  curNbPoints = 0;
  ppts1 = pts1;
  ppts2 = pts2;
  pts = fopen("points.txt", "w");
  while (curNbPoints < pts1->index)
    {
      nbPoints = MIN((pts1->index - curNbPoints), 10000);
      points1 = cvCreateMat(1,nbPoints,CV_32FC2);
      points2 = cvCreateMat(1,nbPoints,CV_32FC2);
      new_points1 = cvCreateMat(2,nbPoints,CV_32FC1);
      new_points2 = cvCreateMat(2,nbPoints,CV_32FC1);

      for( i = 0 ; i < nbPoints ; i++ )
	{
	  points1->data.fl[i*2] = ((Cam2dPoint *)(ppts1->data))->x;
	  points1->data.fl[i*2+1] = ((Cam2dPoint *)(ppts1->data))->y;
	  points2->data.fl[i*2] = ((Cam2dPoint *)(ppts2->data))->x;
	  points2->data.fl[i*2+1] = ((Cam2dPoint *)(ppts2->data))->y;
	  ppts1 = ppts1->next;
	  ppts2 = ppts2->next;
	}

      printf("before %i\n", curNbPoints);
      cvCorrectMatches(fundamental_matrix, points1, points2, NULL, NULL);
      points4d = cvCreateMat(4,nbPoints,CV_32FC1);
      printf("after %i\n", curNbPoints);      

      for (i = 0 ; i < nbPoints ; ++i)
	{
	  cvmSet(new_points1,0,i,points1->data.fl[i*2]);
	  cvmSet(new_points1,1,i,points1->data.fl[i*2+1]);
	  cvmSet(new_points2,0,i,points2->data.fl[i*2]);
	  cvmSet(new_points2,1,i,points2->data.fl[i*2+1]);
	}
      
      cvTriangulatePoints(proj1, proj2, new_points1, new_points2, points4d);

      for (i = 0 ; i < nbPoints ; ++i)
	{
	  tmp3d.x = (POINTS_TYPE)cvmGet(points4d,0,i);
	  tmp3d.y = (POINTS_TYPE)cvmGet(points4d,1,i);
	  tmp3d.z = (POINTS_TYPE)cvmGet(points4d,2,i);
	  write_point(pts, &tmp3d, points1->data.fl + i*2, points2->data.fl + i*2);
	  fwrite(&tmp3d, sizeof(POINTS_TYPE), 3, file);
	}

      curNbPoints += nbPoints;
      cvReleaseMat(&points1);
      cvReleaseMat(&points2);
      cvReleaseMat(&new_points1);
      cvReleaseMat(&new_points2);
      cvReleaseMat(&points4d);
    }
  fclose (pts);

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

  cvReleaseMat(&proj1);
  cvReleaseMat(&proj2);
  cvReleaseMat(&status);
  cvReleaseMat(&fundamental_matrix);
  cvReleaseMat(&new_points1);
  cvReleaseMat(&new_points2);
  cam_disallocate_projections_pair(&projectionPair);
  cam_disallocate_linked_list(pts1);
  cam_disallocate_linked_list(pts2);
  cam_disallocate_linked_list(points);
  cam_disallocate_matrix(&t1);
  cam_disallocate_matrix(&t2);
  cam_disallocate_matrix(&K);
  cam_disallocate_matrix(&F);
  fclose(file);
}

/***********************************************************/
/* Begin of reconstruction with everything known (perfect) */
/***********************************************************/

CamList			*triangulate_2d_points_perfect_2_views(POINTS_TYPE *Kdata,
							       POINTS_TYPE *Rdata1, POINTS_TYPE *Tdata1,
							       POINTS_TYPE *Rdata2, POINTS_TYPE *Tdata2,
							       char *filename)
{
  CamMatrix		K;
  CamMatrix		*R;
  CamMatrix		t1;
  CamMatrix		t2;
  CamList		*pts1;
  CamList		*pts2;
  CamList		*ppts1;
  CamList		*ppts2;
  CamList		*points;
  CamList		*res;
  Cam3dPoint		*pt3d;
  FILE			*file;
  char			*path1;
  char			*path2;
  CamProjectionsPair	projectionPair;

  path1 = (char*)malloc((strlen("data/3d/3dmodels/") + strlen(filename) + 1) * sizeof (char));
  path2 = (char*)malloc((strlen("data/3d/results/") + strlen(filename) + 1) * sizeof (char));
  sprintf(path1, "data/3d/3dmodels/%s", filename);
  sprintf(path2, "data/3d/results/%s", filename);
  points = loadPoints1(path1);
  file = fopen(path2, "w");
  free(path1);
  free(path2);

  cam_allocate_matrix(&K, 3, 3);
  memcpy(K.data, Kdata, 9 * sizeof(POINTS_TYPE));
  cam_allocate_matrix(&projectionPair.p1, 4, 3);
  cam_allocate_matrix(&projectionPair.p2, 4, 3);
  cam_allocate_matrix(&t1, 1, 3);
  cam_allocate_matrix(&t2, 1, 3);  

  memcpy(t1.data, Tdata1, 3 * sizeof(POINTS_TYPE));
  R = compute_rotation_matrix(Rdata1[0], Rdata1[1], Rdata1[2]);
  cam_compute_projection_matrix(&projectionPair.p1, &K, R, &t1);
  pts1 = cam_project_3d_to_2d(points, &projectionPair.p1);
  cam_disallocate_matrix(R);
  free(R);
  memcpy(t2.data, Tdata2, 3 * sizeof(POINTS_TYPE));
  R = compute_rotation_matrix(Rdata2[0], Rdata2[0], Rdata2[0]);
  cam_compute_projection_matrix(&projectionPair.p2, &K, R, &t2);
  pts2 = cam_project_3d_to_2d(points, &projectionPair.p2);
  cam_disallocate_matrix(R);
  free(R);

  ppts1 = pts1;
  ppts2 = pts2;
  res = NULL;
  while (ppts1 && ppts2)
    {
      pt3d = cam_triangulate_one_perfect_3d_point(&projectionPair, &t1, &t2, &K, ppts1->data, ppts2->data);
      if (pt3d)
	{
	  res = cam_add_to_linked_list(res, malloc(sizeof(sbaPoint)));
	  memcpy(&((sbaPoint *)res->data)->pt3d, pt3d, sizeof(Cam3dPoint));
	  ((sbaPoint *)res->data)->mask = NULL;
	  ((sbaPoint *)res->data)->mask = cam_add_to_linked_list(((sbaPoint *)res->data)->mask, malloc(sizeof(mask2d)));
	  ((mask2d *)((sbaPoint *)res->data)->mask->data)->frameIndex = 0;
	  memcpy(&((mask2d *)((sbaPoint *)res->data)->mask->data)->pt2d, ppts1->data, sizeof(Cam2dPoint));
	  ((sbaPoint *)res->data)->mask = cam_add_to_linked_list(((sbaPoint *)res->data)->mask, malloc(sizeof(mask2d)));
	  ((mask2d *)((sbaPoint *)res->data)->mask->data)->frameIndex = 1;
	  memcpy(&((mask2d *)((sbaPoint *)res->data)->mask->data)->pt2d, ppts2->data, sizeof(Cam2dPoint));
	  free (pt3d);
	}
      ppts1 = ppts1->next;
      ppts2 = ppts2->next;
    }

  path1 = (char*)malloc((strlen("data/3d/results/") + strlen(filename) + strlen("_viewpoint_X.ppm")  + 1) * sizeof (char));
  sprintf(path1, "data/3d/results/%s_viewpoint_%i.ppm", filename, 1);
  cam_center_2d_points(pts1, 800, 600);
  cam_points_to_ppm(path1, pts1, 800, 600,
			  255, 255, 255,
			  0, 0, 0);
  sprintf(path1, "data/3d/results/%s_viewpoint_%i.ppm", filename, 2);
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

  return (res);
}

void			triangulate_2d_points_perfect_n_views(CamList	*Projections)
{
  CamMatrix		*p1;
  CamMatrix		*p2;
  CamList		*tmpList;
  int			nb;

  tmpList = Projections;
  p1 = (CamMatrix *)tmpList->data;
  tmpList = tmpList->next;
  p2 = (CamMatrix *)tmpList->data;
  tmpList = tmpList->next;
  
  while (p1 && p2)
    {
      
      p1 = p2;
      p2 = (CamMatrix *)tmpList->data;
      tmpList = tmpList->next;
    }
}

void		main_triangulate_2d_points_perfect_n_views()
{
  
}

void			main_triangulate_2d_points_perfect_2_views()
{
  POINTS_TYPE	Kdata[9] = {1.0f,0.0f,0.0f,0.0f,1.0f,0.0f,0.0f,0.0f,1.0f};
  POINTS_TYPE	Tdata1[3] = {0.0f,0.0f,0.0f};  
  POINTS_TYPE	Tdata2[3] = {1.0f,0.0f,0.0f};
  POINTS_TYPE	Rdata1[3] = {0.0f,0.0f,0.0f};
  POINTS_TYPE	Rdata2[3] = {0.0f,0.0f,0.0f};
  FILE		*cameras;	
  CamList	*points;
    
  cameras = fopen("cameras.txt", "w");
  cam_write_quaternion_camera_params(cameras, Rdata1, Tdata1);
  cam_write_quaternion_camera_params(cameras, Rdata2, Tdata2);
  fclose(cameras);
  points = triangulate_2d_points_perfect_2_views(Kdata, Rdata1, Tdata1, Rdata2, Tdata2, "manny");

  //cam_print_sba_3d_points_text(points,"points.txt");
  cam_print_sba(points,"points.txt");

  cam_disallocate_linked_list(points);
}

/*********************************************************/
/* End of reconstruction with everything known (perfect) */
/*********************************************************/

int		main()
{
  /*main_project_and_write_points();*/
  /*main_triangulate_2d_points_tests();*/
  main_triangulate_2d_points_perfect_2_views();
  //main_triangulate_2d_points_perfect_n_views();
  //main_triangulate_2d_points_noisy_2_views();
  /*main_bundle_2d_points_noisy();*/

  return (0);
}
