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

#include <GL/glut.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdarg.h>

//#define CAM_3D_VIWER_DISPLAY_KEYS
//#define CAM_3D_VIWER_DISPLAY_MOUSE
//#define	CAM_3D_DEBUG

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
  POINTS_TYPE	dist;
}		Cam3dPoint;

typedef struct
{
  POINTS_TYPE	*data;
  int		nrows;
  int		ncols;
}		CamMatrix;

typedef enum
  {
    FALSE,
    TRUE
  }	BOOL;

#define FOV	60.0f
#define ZNEAR	0.0f
#define ZFAR	1000.0f
#define POSX	0.0f
#define	POSY	0.0f
#define	POSZ	0.0f
#define ROTX	0.0f
#define ROTY	0.0f
#define ROTZ	PI
#define	WIDTH	1024
#define	HEIGHT	768

#define ZOOM_FACTOR 1.03f
#define TRANSLATION_ACCELERATION 0.3f
#define ROTATI0N_ACCELERATION 0.3f

/* globals */
float		lastx;
float		lasty;
static CamList	*pointsList = NULL;
static Cam3dPoint *sortedPointsList;
static float	width;
static float	height;
static float	posX=POSX, posY=POSY, posZ=POSZ;
static float	rotX=ROTX, rotY=ROTY, rotZ= ROTZ;
static float	fov = FOV;
static float	zNear =	ZNEAR;
static float	zFar = ZFAR;
static BOOL	drawAxis = FALSE;
static BOOL	mouseLeftDown = FALSE;
static BOOL	mouseRightDown = FALSE;
static BOOL	printClosestPoint = FALSE;
static BOOL	printNormale = FALSE;
static CamMatrix currentRotation;
static CamMatrix currentNormal;
static int	lastX;
static int	lastY;
static int	pointIndex = 0;
static float	rightUpMouseXPos = 0.0f;
static float	rightUpMouseYPos = 0.0f;
static float	rightUpMouseZPos = 0.0f;
static int	info, mainWindow;

/* to be deleted*/
float oldX;
float oldY;
float oldZ;

/*************************/
/* Begin list operations */
/*************************/

inline CamList	*cam_3d_viewer_add_to_linked_list(CamList *l, void *data)
{
  CamList	*head;

  head = (CamList*)malloc(sizeof(CamList));
  head->data = data;
  head->next = l;
  if (l)
    head->index = l->index + 1;
  else
    head->index = 1;
  return (head);
}

void		cam_3d_viewer_free_linked_list(CamList *l)
{
  CamList	*ptr;

  ptr = l;
  while (ptr)
    {
      l = l->next;
      free (ptr);
      ptr = l;
    }
}

void		cam_3d_viewer_free_data_in_linked_list(CamList *l)
{
  CamList	*ptr;

  ptr = l;
  while (ptr)
    {
      free (ptr->data);
      ptr = ptr->next;
    }
}

void	cam_3d_viewer_disallocate_linked_list(CamList *l)
{
  cam_3d_viewer_free_data_in_linked_list(l);
  cam_3d_viewer_free_linked_list(l);
}

/***********************/
/* End list operations */
/***********************/

/***************************/
/* Begin matrix operations */
/***************************/
inline void	cam_3d_viewer_allocate_matrix(CamMatrix *m, int ncols, int nrows)
{
  m->ncols = ncols;
  m->nrows = nrows;
  m->data = (POINTS_TYPE *)calloc(nrows * ncols, sizeof(POINTS_TYPE));
}

inline void	cam_3d_viewer_disallocate_matrix(CamMatrix *m)
{
  free(m->data);
  m->data = NULL;
  m->nrows = 0;
  m->ncols = 0;
}

inline void	cam_3d_viewer_matrix_set_value(CamMatrix *m, int x, int y, POINTS_TYPE value)
{
  m->data[y * m->ncols + x] = value;
  return ;
}

inline POINTS_TYPE	cam_3d_viewer_matrix_get_value(CamMatrix *m, int x, int y)
{
  return (m->data[y * m->ncols + x]);
}

inline void	cam_3d_viewer_matrix_add_value(CamMatrix *m, int x, int y, POINTS_TYPE value)
{
  m->data[y * m->ncols + x] += value;
  return ;
}

void		cam_3d_viewer_print_matrix(CamMatrix *mat, char *name)
{
  register int	i;
  register int	j;

  if (name)
    printf("%s\n", name);
  for (j = 0 ; j < mat->nrows ; ++j)
    {
      for (i = 0 ; i < mat->ncols ; ++i)
	{
	  printf("%f\t", cam_3d_viewer_matrix_get_value(mat, i, j));
	}
      printf("\n");
    }
}

inline void	cam_3d_viewer_matrix_multiply(CamMatrix *res, CamMatrix *m1, CamMatrix *m2)
{
  register int	i;
  register int	j;
  register int	k;

  if (!res->data)
    cam_3d_viewer_allocate_matrix(res, m1->nrows, m2->ncols);
  if (m1->ncols != m2->nrows)
    camError("CamMatrixMultiply", "m1->ncols != m2->nrows");
  for (j = 0 ; j < res->nrows ; ++j)
    {
      for (i = 0 ; i < res->ncols ; ++i)
	{
	  for (k = 0 ; k < m1->ncols ; ++k)
	    {
	      cam_3d_viewer_matrix_add_value(res, i, j, cam_3d_viewer_matrix_get_value(m1, k, j) * cam_3d_viewer_matrix_get_value(m2, i, k));
	    }
	}
    }
}

inline void	cam_3d_viewer_matrix_copy(CamMatrix *dst, CamMatrix *src)
{
  register int	i;
  register int	j;

  if (dst->ncols != src->ncols || dst->nrows != src->nrows)
    camError("CamMatrixCopy", "m1->ncols != m2->cols || m1->nrows != m2->nrows");
  for (j = 0 ; j < src->nrows ; ++j)
    {
      for (i = 0 ; i < src->ncols ; ++i)
	{
	  cam_3d_viewer_matrix_set_value(dst, i, j, cam_3d_viewer_matrix_get_value(src, i, j));
	}
    }
}

/***************************/
/* End matrix operations */
/***************************/

void	camera()
{
  CamMatrix	rotx;
  CamMatrix	roty;
  CamMatrix	rotz;
  CamMatrix	tmp1;
  CamMatrix	tmp2;
  CamMatrix	res;
  CamMatrix	vect;

  cam_3d_viewer_allocate_matrix(&rotx, 3, 3);
  cam_3d_viewer_allocate_matrix(&roty, 3, 3);
  cam_3d_viewer_allocate_matrix(&rotz, 3, 3);
  cam_3d_viewer_allocate_matrix(&tmp1, 3, 3);
  cam_3d_viewer_allocate_matrix(&tmp2, 3, 3);
  cam_3d_viewer_allocate_matrix(&res, 1, 3);
  cam_3d_viewer_allocate_matrix(&vect, 1, 3);
  cam_3d_viewer_matrix_set_value(&vect, 0, 0, (POINTS_TYPE)0.0f);
  cam_3d_viewer_matrix_set_value(&vect, 0, 1, (POINTS_TYPE)0.0f);
  cam_3d_viewer_matrix_set_value(&vect, 0, 2, (POINTS_TYPE)1.0f);
  /* Rx */
  cam_3d_viewer_matrix_set_value(&rotx, 0, 0, (POINTS_TYPE)1.0f);
  cam_3d_viewer_matrix_set_value(&rotx, 1, 0, (POINTS_TYPE)0.0f);
  cam_3d_viewer_matrix_set_value(&rotx, 2, 0, (POINTS_TYPE)0.0f);
  cam_3d_viewer_matrix_set_value(&rotx, 0, 1, (POINTS_TYPE)0.0f);
  cam_3d_viewer_matrix_set_value(&rotx, 1, 1, (POINTS_TYPE)cos(rotX));
  cam_3d_viewer_matrix_set_value(&rotx, 2, 1, (POINTS_TYPE)sin(rotX));
  cam_3d_viewer_matrix_set_value(&rotx, 0, 2, (POINTS_TYPE)0.0f);
  cam_3d_viewer_matrix_set_value(&rotx, 1, 2, (POINTS_TYPE)-sin(rotX));
  cam_3d_viewer_matrix_set_value(&rotx, 2, 2, (POINTS_TYPE)cos(rotX));
  /* Ry */
  cam_3d_viewer_matrix_set_value(&roty, 0, 0, (POINTS_TYPE)cos(rotY));
  cam_3d_viewer_matrix_set_value(&roty, 1, 0, (POINTS_TYPE)0.0f);
  cam_3d_viewer_matrix_set_value(&roty, 2, 0, (POINTS_TYPE)-sin(rotY));
  cam_3d_viewer_matrix_set_value(&roty, 0, 1, (POINTS_TYPE)0.0f);
  cam_3d_viewer_matrix_set_value(&roty, 1, 1, (POINTS_TYPE)1.0f);
  cam_3d_viewer_matrix_set_value(&roty, 2, 1, (POINTS_TYPE)0.0f);
  cam_3d_viewer_matrix_set_value(&roty, 0, 2, (POINTS_TYPE)sin(rotY));
  cam_3d_viewer_matrix_set_value(&roty, 1, 2, (POINTS_TYPE)0.0f);
  cam_3d_viewer_matrix_set_value(&roty, 2, 2, (POINTS_TYPE)cos(rotY));
  /* Rz */
  cam_3d_viewer_matrix_set_value(&rotz, 0, 0, (POINTS_TYPE)cos(rotZ));
  cam_3d_viewer_matrix_set_value(&rotz, 1, 0, (POINTS_TYPE)sin(rotZ));
  cam_3d_viewer_matrix_set_value(&rotz, 2, 0, (POINTS_TYPE)0.0f);
  cam_3d_viewer_matrix_set_value(&rotz, 0, 1, (POINTS_TYPE)-sin(rotZ));
  cam_3d_viewer_matrix_set_value(&rotz, 1, 1, (POINTS_TYPE)cos(rotZ));
  cam_3d_viewer_matrix_set_value(&rotz, 2, 1, (POINTS_TYPE)0.0f);
  cam_3d_viewer_matrix_set_value(&rotz, 0, 2, (POINTS_TYPE)0.0f);
  cam_3d_viewer_matrix_set_value(&rotz, 1, 2, (POINTS_TYPE)0.0f);
  cam_3d_viewer_matrix_set_value(&rotz, 2, 2, (POINTS_TYPE)1.0f);

  cam_3d_viewer_matrix_multiply(&tmp1, &rotx, &roty);

  cam_3d_viewer_matrix_multiply(&tmp2, &tmp1, &rotz);
  
  cam_3d_viewer_matrix_copy(&currentRotation, &tmp2);
  cam_3d_viewer_matrix_multiply(&res, &currentRotation, &vect);

  cam_3d_viewer_matrix_copy(&currentNormal, &res);

  glLoadIdentity();
  gluLookAt(posX, posY, posZ,
	    posX + cam_3d_viewer_matrix_get_value(&currentNormal, 0, 0),
	    posY + cam_3d_viewer_matrix_get_value(&currentNormal, 0, 1),
	    posZ + cam_3d_viewer_matrix_get_value(&currentNormal, 0, 2),
	    0.0f, 1.0f, 0.0f);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glViewport(0, 0, width, height);
  gluPerspective(fov, width / height, zNear, zFar);
  glMatrixMode(GL_MODELVIEW);


  cam_3d_viewer_disallocate_matrix(&rotx);
  cam_3d_viewer_disallocate_matrix(&roty);
  cam_3d_viewer_disallocate_matrix(&rotz);
  cam_3d_viewer_disallocate_matrix(&tmp1);
  cam_3d_viewer_disallocate_matrix(&tmp2);
  cam_3d_viewer_disallocate_matrix(&res);
  cam_3d_viewer_disallocate_matrix(&vect);
}

void	resetView()
{
  fov = FOV;
  zNear = ZNEAR;
  zFar = ZFAR;
  posX = POSX;
  posY = POSY;
  posZ = POSZ;
  rotX = ROTX;
  rotY = ROTY;
  rotZ = ROTZ;
  pointIndex = 0;
#ifdef CAM_3D_DEBUG
  printf("Resetting view\n");
#endif
  camera();
}

/*************************/
/* Begin events handling */
/*************************/

void	processKeyboardKeys(unsigned char key, int x, int y)
{
#ifdef CAM_3D_VIWER_DISPLAY_KEYS
  printf("keyboard key : %i\n", (int)key);
#endif
  switch (key)
    {
    case 27:
      exit(0);
    case 'r':
      resetView();
      break;
    case 's':
      posZ -= 1 * TRANSLATION_ACCELERATION;
      break;
    case 'z':
      posZ += 1 * TRANSLATION_ACCELERATION;
      break;
    case 'q':
      posX += 1 * TRANSLATION_ACCELERATION;
      break;
    case 'd':
      posX -= 1 * TRANSLATION_ACCELERATION;
      break;
    case 'a':
      posY += 1 * TRANSLATION_ACCELERATION;
      break;
    case 'e':
      posY -= 1 * TRANSLATION_ACCELERATION;
      break;
    case 'i':
      rotZ += 0.1f;
      break;
    case 'k':
      rotZ -= 0.1f;
      break;
    case 'o':
      rotY += 0.1f;
      break;
    case 'l':
      rotY -= 0.1f;
      break; 
    case 'p':
      rotX += 0.1f;
      break;
    case 'm':
      rotX -= 0.1f;
      break;
    case 'h':
      drawAxis = (drawAxis + 1) % 2;
      break;
    case 'c':
      printClosestPoint = (printClosestPoint + 1) % 2;
      break;
    case 'b':
      pointIndex++;
      break;
    case 'n':
      pointIndex--;
      if (pointIndex < 0)
	pointIndex = 0;
      break;
    case 'v':
      printNormale = (printNormale + 1) % 2;
      break;
   default:
      break;
    }
}

int	camSort3dPoints(const void *p1, const void *p2)
{
  if (((Cam3dPoint*)p1)->dist < ((Cam3dPoint*)p2)->dist)
    return (-1);
  if (((Cam3dPoint*)p1)->dist == ((Cam3dPoint*)p2)->dist)
    return (0);
  return (1);
  //return (((Cam3dPoint*)p1)->dist - ((Cam3dPoint*)p2)->dist);
}

void	processMouseKeys(int button, int state, int x, int y)
{
  float	d1;
  float	dist1;
  float t1;
  float	NX;
  float	NY;
  float	NZ;
  int	i;

#ifdef CAM_3D_VIWER_DISPLAY_MOUSE
  printf("Mouse pos : %i %i\n", x, y);
#endif
  if (button == GLUT_LEFT_BUTTON)
    {
      if (state == GLUT_DOWN)
	{
#ifdef CAM_3D_VIWER_DISPLAY_MOUSE
	  printf("Mouse left down\n");
#endif
	  mouseLeftDown = TRUE;
	}
      if (state == GLUT_UP)
	{
#ifdef CAM_3D_VIWER_DISPLAY_MOUSE
	  printf("Mouse left up\n");
#endif
	  mouseLeftDown = FALSE;
	}
    }
  if (button == GLUT_RIGHT_BUTTON)
    {
      if (state = GLUT_DOWN)
	{
#ifdef CAM_3D_VIWER_DISPLAY_MOUSE
	  printf("Mouse right down\n");
	  mouseRightDown = TRUE;
#endif
	  
	}
      if (state = GLUT_UP)
	{
#ifdef CAM_3D_VIWER_DISPLAY_MOUSE
	  printf("Mouse right down\n");
#endif
	  mouseRightDown = FALSE;
	  rightUpMouseXPos = posX;
	  rightUpMouseYPos = posY;
	  rightUpMouseZPos = posZ;
	  /* to be deleted */
	  oldX = cam_3d_viewer_matrix_get_value(&currentNormal, 0, 0);
	  oldY = cam_3d_viewer_matrix_get_value(&currentNormal, 0, 1);
	  oldZ = cam_3d_viewer_matrix_get_value(&currentNormal, 0, 2);

	  if (printClosestPoint == TRUE)
	    {
	      NX = cam_3d_viewer_matrix_get_value(&currentNormal, 0, 0);
	      NY = cam_3d_viewer_matrix_get_value(&currentNormal, 0, 1);
	      NZ = cam_3d_viewer_matrix_get_value(&currentNormal, 0, 2);

	      i = 0;

	      while (i < pointsList->index)
		{
		  d1 = NX * sortedPointsList[i].x + NY * sortedPointsList[i].y + NZ * sortedPointsList[i].z;
		  t1 = (d1 - NX * posX - NY * posY - NZ * posZ) / (NX * NX + NY * NY + NZ * NZ);
		  
		  dist1 = (posX + t1 * NX - sortedPointsList[i].x) * (posX + t1 * NX - sortedPointsList[i].x) +
		    (posY + t1 * NY - sortedPointsList[i].y) * (posY + t1 * NY - sortedPointsList[i].y) +
		    (posZ + t1 * NZ - sortedPointsList[i].z) * (posZ + t1 * NZ - sortedPointsList[i].z);
		  
		  sortedPointsList[i].dist = dist1;
		  ++i;
		}
	      qsort(sortedPointsList, pointsList->index, sizeof(Cam3dPoint), camSort3dPoints);
	    }
	}
    }
  lastX = x;
  lastY = y;
}

void	processMouseMotion(int x, int y)
{
#ifdef CAM_3D_VIWER_DISPLAY_MOUSE
  printf("Mouse moved to : %i %i\n", x, y);
#endif
  if (mouseLeftDown == TRUE)
    {
      rotY += (float)(x - lastX) * ROTATI0N_ACCELERATION / (float)500;
      rotX -= (float)(y - lastY) * ROTATI0N_ACCELERATION / (float)500;
    }
  lastX = x;
  lastY = y;
}

void	processMouseWheel(int button, int dir, int x, int y)
{
  if (dir > 0)
    {
      fov /= ZOOM_FACTOR;
#ifdef CAM_3D_VIWER_DISPLAY_MOUSE
      printf("Wheel up, fov : %f\n", fov);
#endif
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      glViewport(0, 0, width, height);
      gluPerspective(fov, width / height, zNear, zFar);
      glMatrixMode(GL_MODELVIEW);
    }
  else
    {
      fov *= ZOOM_FACTOR;
#ifdef CAM_3D_VIWER_DISPLAY_MOUSE
      printf("Wheel down, fov : %f\n", fov);
#endif
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      glViewport(0, 0, width, height);
      gluPerspective(fov, width / height, zNear, zFar);
      glMatrixMode(GL_MODELVIEW);
    }
  lastX = x;
  lastY = y;
}

void	processMousePassiveMotion(int x, int y)
{
#ifdef CAM_3D_VIWER_DISPLAY_MOUSE
  printf("Mouse pos : %i %i\n", x, y);
#endif
}

void		processSpecialKeys(int key, int x, int y)
{
  CamMatrix	vect;
  CamMatrix	tmp;

  cam_3d_viewer_allocate_matrix(&vect, 1, 3);
  cam_3d_viewer_allocate_matrix(&tmp, 1, 3);
  switch (key)
    {
    case GLUT_KEY_LEFT :
#ifdef CAM_3D_VIWER_DISPLAY_KEYS
      printf("Left key pressed\n");
#endif
      cam_3d_viewer_matrix_set_value(&vect, 0, 0, (POINTS_TYPE)-1.0f * TRANSLATION_ACCELERATION);
      cam_3d_viewer_matrix_set_value(&vect, 0, 1, (POINTS_TYPE)0.0f);
      cam_3d_viewer_matrix_set_value(&vect, 0, 2, (POINTS_TYPE)0.0f);
      cam_3d_viewer_matrix_multiply(&tmp, &currentRotation, &vect);
      posX += cam_3d_viewer_matrix_get_value(&tmp, 0, 0);
      posY += cam_3d_viewer_matrix_get_value(&tmp, 0, 1);
      posZ += cam_3d_viewer_matrix_get_value(&tmp, 0, 2);
      break;
    case GLUT_KEY_RIGHT :
#ifdef CAM_3D_VIWER_DISPLAY_KEYS
      printf("Right key pressed\n");
#endif
      cam_3d_viewer_matrix_set_value(&vect, 0, 0, (POINTS_TYPE)1.0f * TRANSLATION_ACCELERATION);
      cam_3d_viewer_matrix_set_value(&vect, 0, 1, (POINTS_TYPE)0.0f);
      cam_3d_viewer_matrix_set_value(&vect, 0, 2, (POINTS_TYPE)0.0f);
      cam_3d_viewer_matrix_multiply(&tmp, &currentRotation, &vect);
      posX += cam_3d_viewer_matrix_get_value(&tmp, 0, 0);
      posY += cam_3d_viewer_matrix_get_value(&tmp, 0, 1);
      posZ += cam_3d_viewer_matrix_get_value(&tmp, 0, 2);
      break;
    case GLUT_KEY_UP :
#ifdef CAM_3D_VIWER_DISPLAY_KEYS
      printf("Up key pressed\n");
#endif
      cam_3d_viewer_matrix_set_value(&vect, 0, 0, (POINTS_TYPE)0.0f);
      cam_3d_viewer_matrix_set_value(&vect, 0, 1, (POINTS_TYPE)0.0f);
      cam_3d_viewer_matrix_set_value(&vect, 0, 2, (POINTS_TYPE)1.0f * TRANSLATION_ACCELERATION);
      cam_3d_viewer_matrix_multiply(&tmp, &currentRotation, &vect);
      posX += cam_3d_viewer_matrix_get_value(&tmp, 0, 0);
      posY += cam_3d_viewer_matrix_get_value(&tmp, 0, 1);
      posZ += cam_3d_viewer_matrix_get_value(&tmp, 0, 2);
      break;
    case GLUT_KEY_DOWN :
#ifdef CAM_3D_VIWER_DISPLAY_KEYS
      printf("Down key pressed\n");
#endif
      cam_3d_viewer_matrix_set_value(&vect, 0, 0, (POINTS_TYPE)0.0f);
      cam_3d_viewer_matrix_set_value(&vect, 0, 1, (POINTS_TYPE)0.0f);
      cam_3d_viewer_matrix_set_value(&vect, 0, 2, (POINTS_TYPE)-1.0f * TRANSLATION_ACCELERATION);
      cam_3d_viewer_matrix_multiply(&tmp, &currentRotation, &vect);
      posX += cam_3d_viewer_matrix_get_value(&tmp, 0, 0);
      posY += cam_3d_viewer_matrix_get_value(&tmp, 0, 1);
      posZ += cam_3d_viewer_matrix_get_value(&tmp, 0, 2);
      break;
    default:
      break;
    }
  cam_3d_viewer_disallocate_matrix(&vect);
  cam_3d_viewer_disallocate_matrix(&tmp);
}

/*************************/
/* End events handling */
/*************************/

void	infoChangeSize(int w, int h)
{
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, w, h, 0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glClearColor(0.0, 0.0, 0.0, 0.0);
}

void	changeSize(int w, int h)
{
  if (h == 0)
    h = 1;

  width = (float)w;
  height = (float)h;
  gluOrtho2D(0,width,height,0);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glViewport(0, 0, w, h);
  gluPerspective(fov, width / height, zNear, zFar);
  glMatrixMode(GL_MODELVIEW);
  camera();
}


void	redisplay_all(void)
{
  glutSetWindow(info);
  glutPostRedisplay();
  glutSetWindow(mainWindow);
  //  changeSize(WIDTH, HEIGHT);
  glutPostRedisplay();
}



void	drawstr(int x, int y, char* format, ...)
{
  va_list args;
  char buffer[255], *s;
    
  va_start(args, format);
  vsprintf(buffer, format, args);
  va_end(args);
    
  glRasterPos2i(x, y);
  for (s = buffer; *s; s++)
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *s);
}

void	renderInfos(void)
{
  glClearColor(0.2, 0.2, 0.2, 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glColor3f(1.0f, 0.0f, 0.0f);

  glRasterPos2f(10.0f, 10.0f);
  drawstr(10, 10, "Mouse pos :");
  drawstr(10, 25, "x: %f", posX);
  drawstr(10, 40, "y: %f", posY);
  drawstr(10, 55, "z: %f", posZ);

  if (printClosestPoint == TRUE)
    {
      drawstr(10, 85, "Selected point :");
      drawstr(10, 100, "x: %f", sortedPointsList[pointIndex].x);
      drawstr(10, 115, "y: %f", sortedPointsList[pointIndex].y);
      drawstr(10, 130, "z: %f", sortedPointsList[pointIndex].z);
      drawstr(10, 145, "distance: %f", sortedPointsList[pointIndex].dist);
    }

  glutSwapBuffers();  
}

void		renderSparsePoints(void)
{
  int		index;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  if (drawAxis)
    {
      glBegin(GL_LINES);
      // X axis
      glColor3f(1.0f, 0, 0);
      glVertex3f(0, 0, 0);
      glVertex3f(1.0f, 0, 0);
      // Y axis
      glColor3f(0, 1.0f, 0);
      glVertex3f(0, 0, 0);
      glVertex3f(0, 1.0f, 0);
      // Z axis
      glColor3f(0, 0, 1.0f);
      glVertex3f(0, 0, 0);
      glVertex3f(0, 0, 1.0f);
      glEnd();
    }

  glBegin(GL_POINTS);
  index = 0;
  glColor3f(1.0f, 1.0f, 1.0f);
  while (index < pointsList->index)
    {
      glVertex3f(sortedPointsList[index].x,
		 sortedPointsList[index].y,
		 sortedPointsList[index].z);
      ++index;
    }
  glEnd();
  
  if (printClosestPoint == TRUE)
    {
      glBegin(GL_LINES);
      glColor3f(0.5f, 1.0f, 0.0f);
      glVertex3f(rightUpMouseXPos,
		 rightUpMouseYPos,
		 rightUpMouseZPos);

      glVertex3f(sortedPointsList[pointIndex].x,
		 sortedPointsList[pointIndex].y,
		 sortedPointsList[pointIndex].z);

      if (printNormale == TRUE)
	{
	  glColor3f(1.0f, 1.0f, 0.0f);
	  glVertex3f(rightUpMouseXPos,
		     rightUpMouseYPos,
		     rightUpMouseZPos);
	  
	  glVertex3f(rightUpMouseXPos + oldX * ZFAR,
		     rightUpMouseYPos + oldY * ZFAR,
		     rightUpMouseZPos + oldZ * ZFAR);
	}

      glEnd();
    }

  camera();

  redisplay_all();
  glutSwapBuffers();
  return ;
}

void	cam_3d_viewer_init_display(int xPos, int yPos, int width, int height,
				   unsigned int displayMode, char *windowName,
				   void (*keyboardFuncCallback)(unsigned char, int, int),
				   void (*specialFuncCallback)(int, int, int),
				   void (*mouseFuncCallback)(int, int, int, int),
				   void (*mouseWheelCallback)(int, int, int , int),
 				   void	(*motionCallback)(int, int),
 				   void	(*passiveMotionCallback)(int, int),
				   void (*displayCallback)(void),
				   void (*idleCallback)(void),
				   void (*reshapeCallback)(int, int),
				   int infoWidth, int infoHeight,
				   void (*infoDisplayCallback)(void),
				   void (*infoReshapeCallback)(int, int))
{
  int	ac;

  ac = 0;
  glutInit(&ac, NULL);
  glutInitDisplayMode(displayMode);
  glutInitWindowPosition(xPos, yPos);
  glutInitWindowSize(width, height);
  mainWindow = glutCreateWindow(windowName);

  if (keyboardFuncCallback)
    glutKeyboardFunc(keyboardFuncCallback);
  if (specialFuncCallback)
    glutSpecialFunc(specialFuncCallback);
  glutDisplayFunc(displayCallback);
  if (idleCallback)
    glutIdleFunc(idleCallback);
  if (reshapeCallback)
    glutReshapeFunc(reshapeCallback);
  if (mouseFuncCallback)
    glutMouseFunc(mouseFuncCallback);
  if (mouseWheelCallback)
    glutMouseWheelFunc(mouseWheelCallback);
  if (motionCallback)
    glutMotionFunc(motionCallback);
  if (passiveMotionCallback)
    glutPassiveMotionFunc(passiveMotionCallback);

  info = glutCreateSubWindow(mainWindow, width - infoWidth, height - infoHeight, infoWidth, infoHeight);
  if (infoReshapeCallback)
    glutReshapeFunc(infoReshapeCallback);
  glutDisplayFunc(infoDisplayCallback);
  
  redisplay_all();
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

void	loadAyetPoints(char *file)
{
  FILE	*dataFd;
  int	spacingSymbol;

  dataFd = fopen(file, "r");
  if (!dataFd)
    {
      printf("Unable to open the file %s, which is containing the data\n", file);
      perror("");
      exit(-1);
    }
  while (1)
    {
      pointsList = cam_3d_viewer_add_to_linked_list(pointsList, (Cam3dPoint *)malloc(sizeof(Cam3dPoint)));
      ((Cam3dPoint *)pointsList->data)->x = extractNextValue(dataFd);
      ((Cam3dPoint *)pointsList->data)->y = extractNextValue(dataFd);
      ((Cam3dPoint *)pointsList->data)->z = extractNextValue(dataFd);
      fgetc(dataFd); // \r
      fgetc(dataFd); // \n
      spacingSymbol = fgetc(dataFd);
      if (spacingSymbol == EOF)
	break;
      else
	fseek(dataFd, -1, SEEK_CUR);
    }
  printf("Load done\n");
  fclose(dataFd);
}

void		loadSortedPointsList(CamList *list, Cam3dPoint *sortedPointsList)
{
  CamList	*ptr;
  int		index;

  ptr = list;
  index = 0;
  while (ptr)
    {
      memcpy(sortedPointsList + index, ptr->data, sizeof(Cam3dPoint));
      ptr = ptr->next;
      ++index;
    }
}

void	test_cam_3d_viewer()
{
  cam_3d_viewer_allocate_matrix(&currentRotation, 3, 3);
  cam_3d_viewer_allocate_matrix(&currentNormal, 1, 3);
  loadAyetPoints("/home/splin/manny"); // methode de chargement spécifique à chaque format de fichier
  sortedPointsList = (Cam3dPoint *)malloc(pointsList->index * sizeof(Cam3dPoint));
  loadSortedPointsList(pointsList, sortedPointsList);
  cam_3d_viewer_init_display(100, 100, WIDTH, HEIGHT, GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB, "Camellia Vizualizer",
			     processKeyboardKeys, processSpecialKeys, processMouseKeys, processMouseWheel, processMouseMotion, processMousePassiveMotion,
			     renderSparsePoints, renderSparsePoints, changeSize,
			     150, 150,
			     renderInfos, infoChangeSize);
  glutMainLoop();
  return ;
}
