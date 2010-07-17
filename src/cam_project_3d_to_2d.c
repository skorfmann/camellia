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
#include	"camellia.h"

#define PRINT_MATRIX
#define	PRINT_VECTORS
#define	POINTS_TYPE	float
#define PI		3.1415926535897932384626433832795

typedef struct		s_camList
{
  void			*data;
  struct s_camList	*next;
  int			index;
}			CamList;

typedef struct
{
  POINTS_TYPE	x;
  POINTS_TYPE	y;
  POINTS_TYPE	z;
}		Cam3dPoint;

typedef struct
{
  POINTS_TYPE	*data;
  int		nrows;
  int		ncols;
}		CamMatrix;

typedef struct
{
  POINTS_TYPE	*data;
  int		nbElems;
}		CamVector;

inline void	CamAllocateMatrix(CamMatrix *m, int ncols, int nrows)
{
  m->ncols = ncols;
  m->nrows = nrows;
  m->data = (POINTS_TYPE *)calloc(nrows * ncols, sizeof(POINTS_TYPE));
}

inline void	CamDisallocateMatrix(CamMatrix *m)
{
  free(m->data);
  m->data = NULL;
  m->nrows = 0;
  m->ncols = 0;
}

inline void	CamAllocateVector(CamVector *m, int nbElems)
{
  m->nbElems = nbElems;
  m->data = (POINTS_TYPE *)malloc(nbElems * sizeof(POINTS_TYPE));
}

inline void	CamDisallocateVector(CamVector *m)
{
  free(m->data);
  m->data = NULL;
  m->nbElems = 0;
}

inline void	CamMatrixSetValue(CamMatrix *m, int x, int y, POINTS_TYPE value)
{
  m->data[y * m->ncols + x] = value;
  return ;
}

inline void	CamMatrixAddValue(CamMatrix *m, int x, int y, POINTS_TYPE value)
{
  m->data[y * m->ncols + x] += value;
  return ;
}


inline POINTS_TYPE	CamMatrixGetValue(CamMatrix *m, int x, int y)
{
  return (m->data[y * m->ncols + x]);
}

inline void	CamVectorSetValue(CamMatrix *m, int index, POINTS_TYPE value)
{
  m->data[index] = value;
  return ;
}

inline POINTS_TYPE	CamVectorGetValue(CamVector *m, int index)
{
  return (m->data[index]);
}

inline void	CamMatrixAdd(CamMatrix *res, CamMatrix *m1, CamMatrix *m2)
{
  register int	i;
  register int	j;
  
  if (m1->ncols != m2->ncols || m1->nrows != m2->nrows)
    camError("CamMatrixAdd", "m1->ncols != m2->cols || m1->nrows != m2->nrows");
  for (j = 0 ; j < m1->nrows ; ++j)
    {
      for (i = 0 ; i < m1->ncols ; ++i)
	{
	  CamMatrixSetValue(res, i, j, CamMatrixGetValue(m1, i, j) + CamMatrixGetValue(m2, i, j));
	}
    }
}

inline void	CamMatrixCopy(CamMatrix *dst, CamMatrix *src)
{
  register int	i;
  register int	j;

  if (dst->ncols != src->ncols || dst->nrows != src->nrows)
    camError("CamMatrixCopy", "m1->ncols != m2->cols || m1->nrows != m2->nrows");
  for (j = 0 ; j < src->nrows ; ++j)
    {
      for (i = 0 ; i < src->ncols ; ++i)
	{
	  CamMatrixSetValue(dst, i, j, CamMatrixGetValue(src, i, j));
	}
    }
}

inline void	CamMatrixMultiply(CamMatrix *res, CamMatrix *m1, CamMatrix *m2)
{
  register int	i;
  register int	j;
  register int	k;

  if (!res->data)
    CamAllocateMatrix(res, m1->nrows, m2->ncols);
  if (m1->ncols != m2->nrows)
    camError("CamMatrixMultiply", "m1->ncols != m2->nrows");
  for (j = 0 ; j < res->nrows ; ++j)
    {
      for (i = 0 ; i < res->ncols ; ++i)
	{
	  for (k = 0 ; k < m1->ncols ; ++k)
	    {
	      CamMatrixAddValue(res, i, j, CamMatrixGetValue(m1, k, j) * CamMatrixGetValue(m2, i, k));
	    }
	}
    }
}

#ifdef	PRINT_MATRIX
void		CamPrintMatrix(CamMatrix *mat, char *name)
{
  register int	i;
  register int	j;

  if (name)
    printf("%s\n", name);
  for (j = 0 ; j < mat->nrows ; ++j)
    {
      for (i = 0 ; i < mat->ncols ; ++i)
	{
	  printf("%f\t", CamMatrixGetValue(mat, i, j));
	}
      printf("\n");
    }
}
#endif

#ifdef	PRINT_VECTOR
void		CamPrintVector(CamVector *vect, char *name)
{
  register int	i;

  if (name)
    printf("%s\n", name);
  for (i = 0 ; i < vect->nbElems ; ++i)
    {
      printf("%f\t", CamVectorGetValue(vect, i));
    }
  printf("\n");
}
#endif

/* absolute translation and rotation in the 3d space */
CamList		*cam_project_3d_to_2d(CamList *points, CamMatrix *K, CamMatrix *R, CamVector *t)
{
  CamList	*res;
  CamMatrix	P;
  CamMatrix	Rt;
  register int	i;
  register int	j;

  CamAllocateMatrix(&Rt, 4, 3);
  for (j = 0 ; j < 3 ; ++j)
    {
      for (i = 0 ; i < 3 ; ++i)
	{
	  CamMatrixSetValue(&Rt, i, j, CamMatrixGetValue(R, i, j));
	}
      CamMatrixSetValue(&Rt, i, j, CamVectorGetValue(t, j));
    }
#ifdef PRINT_MATRIX
  CamPrintMatrix(K, "Calibration");
  CamPrintMatrix(R, "Rotation");
#endif
#ifdef PRINT_VECTOR
  CamPrintVector(t, "Translation");
#endif
#ifdef PRINT_MATRIX
  CamPrintMatrix(&Rt, "Rotation + Translation");
#endif
  CamAllocateMatrix(&P, 4, 3);
  res = (CamList *)malloc(points->index * sizeof(CamList));
  CamMatrixMultiply(&P, K, &Rt);
#ifdef PRINT_MATRIX
  CamPrintMatrix(&P, "Projection Matrix");
#endif
  CamDisallocateMatrix(&Rt);
  CamDisallocateMatrix(&P);
  return (res);
}

CamMatrix	*compute_rotation_matrix(POINTS_TYPE rx, POINTS_TYPE ry, POINTS_TYPE rz)
{
  CamMatrix	rotx;
  CamMatrix	roty;
  CamMatrix	rotz;
  CamMatrix	*res;

  CamAllocateMatrix(&rotx, 3, 3);
  CamAllocateMatrix(&roty, 3, 3);
  CamAllocateMatrix(&rotz, 3, 3);
  res = (CamMatrix *)malloc(sizeof(CamMatrix));
  CamAllocateMatrix(res, 3, 3);
  /* Rx */
  CamMatrixSetValue(&rotx, 0, 0, (POINTS_TYPE)1.0f);
  CamMatrixSetValue(&rotx, 1, 0, (POINTS_TYPE)0.0f);
  CamMatrixSetValue(&rotx, 2, 0, (POINTS_TYPE)0.0f);
  CamMatrixSetValue(&rotx, 0, 1, (POINTS_TYPE)0.0f);
  CamMatrixSetValue(&rotx, 1, 1, (POINTS_TYPE)cos(rx));
  CamMatrixSetValue(&rotx, 2, 1, (POINTS_TYPE)sin(rx));
  CamMatrixSetValue(&rotx, 0, 2, (POINTS_TYPE)0.0f);
  CamMatrixSetValue(&rotx, 1, 2, (POINTS_TYPE)-sin(rx));
  CamMatrixSetValue(&rotx, 2, 2, (POINTS_TYPE)cos(rx));
#ifdef PRINT_MATRIX
  CamPrintMatrix(&rotx, "Rx");
#endif
  /* Ry */
  CamMatrixSetValue(&roty, 0, 0, (POINTS_TYPE)cos(ry));
  CamMatrixSetValue(&roty, 1, 0, (POINTS_TYPE)0.0f);
  CamMatrixSetValue(&roty, 2, 0, (POINTS_TYPE)-sin(ry));
  CamMatrixSetValue(&roty, 0, 1, (POINTS_TYPE)0.0f);
  CamMatrixSetValue(&roty, 1, 1, (POINTS_TYPE)1.0f);
  CamMatrixSetValue(&roty, 2, 1, (POINTS_TYPE)0.0f);
  CamMatrixSetValue(&roty, 0, 2, (POINTS_TYPE)sin(ry));
  CamMatrixSetValue(&roty, 1, 2, (POINTS_TYPE)0.0f);
  CamMatrixSetValue(&roty, 2, 2, (POINTS_TYPE)cos(ry));
#ifdef PRINT_MATRIX
  CamPrintMatrix(&roty, "Ry");
#endif
  /* Rz */
  CamMatrixSetValue(&rotz, 0, 0, (POINTS_TYPE)cos(rz));
  CamMatrixSetValue(&rotz, 1, 0, (POINTS_TYPE)sin(rz));
  CamMatrixSetValue(&rotz, 2, 0, (POINTS_TYPE)0.0f);
  CamMatrixSetValue(&rotz, 0, 1, (POINTS_TYPE)-sin(rz));
  CamMatrixSetValue(&rotz, 1, 1, (POINTS_TYPE)cos(rz));
  CamMatrixSetValue(&rotz, 2, 1, (POINTS_TYPE)0.0f);
  CamMatrixSetValue(&rotz, 0, 2, (POINTS_TYPE)0.0f);
  CamMatrixSetValue(&rotz, 1, 2, (POINTS_TYPE)0.0f);
  CamMatrixSetValue(&rotz, 2, 2, (POINTS_TYPE)1.0f);
#ifdef PRINT_MATRIX
  CamPrintMatrix(&rotz, "Rz");
#endif

  CamMatrixCopy(res, &rotx);
  //  CamMatrixAdd(res, res, &roty);
  //CamMatrixAdd(res, res, &rotz);

  CamDisallocateMatrix(&rotx);
  CamDisallocateMatrix(&roty);
  CamDisallocateMatrix(&rotz);
  return (res);
}

void		test_cam_project_3d_to_2d()
{
  CamMatrix	K;
  CamMatrix	*R;
  CamVector	t;
  CamList	*points;
  POINTS_TYPE	Kdata[9] = {1,0,0,0,1,0,0,0,1};
  POINTS_TYPE	Tdata[3] = {0,0,0};
  CamList	*res;

  CamAllocateMatrix(&K, 3, 3);
  R = compute_rotation_matrix(PI / 4, PI / 4, 0.0f);
  CamAllocateVector(&t, 3);
  memcpy(K.data, Kdata, 9 * sizeof(POINTS_TYPE));
  memcpy(t.data, Tdata, 3 * sizeof(POINTS_TYPE));
  points = NULL;
  points = cam_keypoints_tracking2_add_to_linked_list(points, (void*)1);
  points = cam_keypoints_tracking2_add_to_linked_list(points, (void*)4);
  res = cam_project_3d_to_2d(points, &K, R, &t);
  cam_keypoints_tracking2_free_linked_list(points);
  free(res);
  CamDisallocateVector(&t);
  CamDisallocateMatrix(R);
  CamDisallocateMatrix(&K);
  free(R);
}
