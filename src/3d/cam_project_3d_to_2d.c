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

#include	<stdio.h>
#include	<math.h>
#include	<string.h>
#include	<stdlib.h>
#include	"cam_list.h"
#include	"cam_matrix.h"
#include	"cam_vector.h"
#include	"cam_2d_points.h"
#include	"cam_3d_points.h"
#include	"cam_write_points_to_pgm.h"
#include	"misc.h"

#define PRINT_MATRIX
#define	PRINT_VECTOR
#define PI		3.1415926535897932384626433832795

/* absolute translation and rotation in the 3d space */
CamList		*cam_project_3d_to_2d(CamList *points, CamMatrix *K, CamMatrix *R, CamVector *t)
{
  CamList	*res;
  CamList	*pts;
  CamMatrix	pt3d;
  CamMatrix	pt2d;
  CamMatrix	P;
  CamMatrix	Rt;
  register int	i;
  register int	j;

  cam_allocate_matrix(&Rt, 4, 3);
  for (j = 0 ; j < 3 ; ++j)
    {
      for (i = 0 ; i < 3 ; ++i)
	{
	  cam_matrix_set_value(&Rt, i, j, cam_matrix_get_value(R, i, j));
	}
      cam_matrix_set_value(&Rt, i, j, cam_vector_get_value(t, j));
    }
#ifdef PRINT_MATRIX
  cam_print_matrix(K, "Calibration");
  cam_print_matrix(R, "Rotation");
#endif
#ifdef PRINT_VECTOR
  cam_print_vector(t, "Translation");
#endif
#ifdef PRINT_MATRIX
  cam_print_matrix(&Rt, "Rotation + Translation");
#endif
  cam_allocate_matrix(&P, 4, 3);
  res = NULL;
  cam_matrix_multiply(&P, K, &Rt);
#ifdef PRINT_MATRIX
  cam_print_matrix(&P, "Projection Matrix");
#endif

  cam_allocate_matrix(&pt3d, 1, 4);
  cam_allocate_matrix(&pt2d, 1, 3);
  cam_matrix_set_value(&pt3d, 0, 3, 1.0f);
  pts = points;
  while (pts)
    {
      cam_matrix_set_value(&pt3d, 0, 0, ((Cam3dPoint *)(pts->data))->x);
      cam_matrix_set_value(&pt3d, 0, 1, ((Cam3dPoint *)(pts->data))->y);
      cam_matrix_set_value(&pt3d, 0, 2, ((Cam3dPoint *)(pts->data))->z);
      cam_matrix_multiply(&pt2d, &P, &pt3d);
      res = cam_add_to_linked_list(res, (Cam2dPoint *)malloc(sizeof(Cam2dPoint)));
      ((Cam2dPoint *)(res->data))->x = cam_matrix_get_value(&pt2d, 0, 0);
      ((Cam2dPoint *)(res->data))->y = cam_matrix_get_value(&pt2d, 0, 1);
      pts = pts->next;
    }

  cam_disallocate_matrix(&pt3d);
  cam_disallocate_matrix(&Rt);
  cam_disallocate_matrix(&P);
  return (res);
}

CamMatrix	*compute_rotation_matrix(POINTS_TYPE rx, POINTS_TYPE ry, POINTS_TYPE rz)
{
  CamMatrix	rotx;
  CamMatrix	roty;
  CamMatrix	rotz;
  CamMatrix	*res;

  cam_allocate_matrix(&rotx, 3, 3);
  cam_allocate_matrix(&roty, 3, 3);
  cam_allocate_matrix(&rotz, 3, 3);
  res = (CamMatrix *)malloc(sizeof(CamMatrix));
  cam_allocate_matrix(res, 3, 3);
  /* Rx */
  cam_matrix_set_value(&rotx, 0, 0, (POINTS_TYPE)1.0f);
  cam_matrix_set_value(&rotx, 1, 0, (POINTS_TYPE)0.0f);
  cam_matrix_set_value(&rotx, 2, 0, (POINTS_TYPE)0.0f);
  cam_matrix_set_value(&rotx, 0, 1, (POINTS_TYPE)0.0f);
  cam_matrix_set_value(&rotx, 1, 1, (POINTS_TYPE)cos(rx));
  cam_matrix_set_value(&rotx, 2, 1, (POINTS_TYPE)sin(rx));
  cam_matrix_set_value(&rotx, 0, 2, (POINTS_TYPE)0.0f);
  cam_matrix_set_value(&rotx, 1, 2, (POINTS_TYPE)-sin(rx));
  cam_matrix_set_value(&rotx, 2, 2, (POINTS_TYPE)cos(rx));
#ifdef PRINT_MATRIX
  cam_print_matrix(&rotx, "Rx");
#endif
  /* Ry */
  cam_matrix_set_value(&roty, 0, 0, (POINTS_TYPE)cos(ry));
  cam_matrix_set_value(&roty, 1, 0, (POINTS_TYPE)0.0f);
  cam_matrix_set_value(&roty, 2, 0, (POINTS_TYPE)-sin(ry));
  cam_matrix_set_value(&roty, 0, 1, (POINTS_TYPE)0.0f);
  cam_matrix_set_value(&roty, 1, 1, (POINTS_TYPE)1.0f);
  cam_matrix_set_value(&roty, 2, 1, (POINTS_TYPE)0.0f);
  cam_matrix_set_value(&roty, 0, 2, (POINTS_TYPE)sin(ry));
  cam_matrix_set_value(&roty, 1, 2, (POINTS_TYPE)0.0f);
  cam_matrix_set_value(&roty, 2, 2, (POINTS_TYPE)cos(ry));
#ifdef PRINT_MATRIX
  cam_print_matrix(&roty, "Ry");
#endif
  /* Rz */
  cam_matrix_set_value(&rotz, 0, 0, (POINTS_TYPE)cos(rz));
  cam_matrix_set_value(&rotz, 1, 0, (POINTS_TYPE)sin(rz));
  cam_matrix_set_value(&rotz, 2, 0, (POINTS_TYPE)0.0f);
  cam_matrix_set_value(&rotz, 0, 1, (POINTS_TYPE)-sin(rz));
  cam_matrix_set_value(&rotz, 1, 1, (POINTS_TYPE)cos(rz));
  cam_matrix_set_value(&rotz, 2, 1, (POINTS_TYPE)0.0f);
  cam_matrix_set_value(&rotz, 0, 2, (POINTS_TYPE)0.0f);
  cam_matrix_set_value(&rotz, 1, 2, (POINTS_TYPE)0.0f);
  cam_matrix_set_value(&rotz, 2, 2, (POINTS_TYPE)1.0f);
#ifdef PRINT_MATRIX
  cam_print_matrix(&rotz, "Rz");
#endif

  cam_matrix_copy(res, &rotx);
  cam_matrix_multiply(&rotx, res, &roty);
  cam_matrix_multiply(res, &rotx, &rotz);

  cam_disallocate_matrix(&rotx);
  cam_disallocate_matrix(&roty);
  cam_disallocate_matrix(&rotz);
  return (res);
}

POINTS_TYPE	extractNextValue(FILE *fd)
{
  POINTS_TYPE	res;
  POINTS_TYPE	sign;
  char		symbol;
  BOOL		integerPart;
  POINTS_TYPE	power;

  res = 0;
  symbol = (char)fgetc(fd);
  integerPart = TRUE;
  power = (POINTS_TYPE)10;
  if (symbol == '-')
    sign = (POINTS_TYPE)-1;
  else
    {
      sign = (POINTS_TYPE)1;
      res = (POINTS_TYPE)symbol - '0';
    }
  while (1)
    {
      symbol = (char)fgetc(fd);
      if (symbol == ' ')
	break;
      if (symbol == '.')
	{
	  integerPart = FALSE;
	}
      else
	{
	  if (integerPart == TRUE)
	    {
	      res *= (POINTS_TYPE)10;
	      res += (POINTS_TYPE)symbol - '0';
	    }
	  else
	    {
	      res += (POINTS_TYPE)(symbol - '0') / power;
	      power *= 10;
	    }
	}
    }
  res *= sign;
  return (res);
}

CamList		*loadPoints(char *file)
{
  FILE		*dataFd;
  int		spacingSymbol;
  CamList	*pointsList;

  dataFd = fopen(file, "r");
  pointsList = NULL;
  if (!dataFd)
    {
      printf("Unable to open the file %s, which is containing the data\n", file);
      perror("");
      exit(-1);
    }
  while (1)
    {
      pointsList = cam_add_to_linked_list(pointsList, (Cam3dPoint *)malloc(sizeof(Cam3dPoint)));
      ((Cam3dPoint *)pointsList->data)->x = extractNextValue(dataFd);
      ((Cam3dPoint *)pointsList->data)->y = extractNextValue(dataFd);
      ((Cam3dPoint *)pointsList->data)->z = extractNextValue(dataFd);
      fgetc(dataFd);
      fgetc(dataFd);
      spacingSymbol = fgetc(dataFd);
      if (spacingSymbol == EOF)
	break;
      else
	fseek(dataFd, -1, SEEK_CUR);
    }
  printf("Load done\n");
  fclose(dataFd);
  return (pointsList);
}

int		main()
{
  CamMatrix	K;
  CamMatrix	*R;
  CamVector	t;
  CamList	*points;
  POINTS_TYPE	Kdata[9] = {1,0,0,0,1,0,0,0,1};
  POINTS_TYPE	Tdata[3] = {0,0,0};
  CamList	*res;

  cam_allocate_matrix(&K, 3, 3);
  R = compute_rotation_matrix(0.0f, 0.0f, PI);
  cam_allocate_vector(&t, 3);
  memcpy(K.data, Kdata, 9 * sizeof(POINTS_TYPE));
  memcpy(t.data, Tdata, 3 * sizeof(POINTS_TYPE));

  points = loadPoints("/home/splin/manny");

  res = cam_project_3d_to_2d(points, &K, R, &t);
  cam_write_points_to_pgm("pts.pgm", res, 800, 600,
			  255, 0, 0,
			  0, 255, 0);
  cam_disallocate_linked_list(res);
  cam_disallocate_linked_list(points);
  cam_disallocate_vector(&t);
  cam_disallocate_matrix(R);
  cam_disallocate_matrix(&K);
  free(R);
  return (0);
}
