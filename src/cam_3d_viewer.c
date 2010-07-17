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

#define CAM_3D_VIWER_DISPLAY_KEYS
#define CAM_3D_VIWER_DISPLAY_MOUSE
#define	CAM_3D_DEBUG

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

typedef enum
  {
    FALSE,
    TRUE
  }	BOOL;

#define FOV	60
#define ZNEAR	1
#define ZFAR	1000

/* globals */
float lastx;
float lasty;
static CamList	*pointsList = NULL;
static float	ratio;
static float	zoomFactor = 1.03f;
static float	posX=0.0f,posY=1.75f,posZ=5.0f;
static float	lx=0.0f,ly=0.0f,lz=-1.0f;
static float	fov = 60;
static float	zNear =	1;
static float	zFar = 1000;

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

void	resetView()
{
  fov = FOV;
  zNear = ZNEAR;
  zFar = ZFAR;
#ifdef CAM_3D_DEBUG
  printf("Resetting view : fov %f zNear zFar\n", fov, zNear, zFar);
#endif
  gluPerspective (fov, ratio, zNear, zFar); 
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
    case 114:
      resetView();
      break;
    default:
      break;
    }
}

void processMouseWheel(int button, int dir, int x, int y)
{
  if (dir > 0)
    {
      fov *= zoomFactor;
#ifdef CAM_3D_VIWER_DISPLAY_MOUSE
      printf("Wheel up, fov : %f\n", fov);
#endif
      gluPerspective (fov, ratio, zNear, zFar); 
    }
  else
    {
      fov /= zoomFactor;
#ifdef CAM_3D_VIWER_DISPLAY_MOUSE
      printf("Wheel down, fov : %f\n", fov);
#endif
      gluPerspective (fov, ratio, zNear, zFar); 
    }
}

void	processMouseMovement(int x, int y)
{
#ifdef CAM_3D_VIWER_DISPLAY_MOUSE
  printf("Mouse pos : %i %i\n", x, y);
#endif
}

void processSpecialKeys(int key, int x, int y)
{
  switch (key)
    {
    case GLUT_KEY_LEFT :
#ifdef CAM_3D_VIWER_DISPLAY_KEYS
      printf("Left key pressed\n");
#endif
      break;
    case GLUT_KEY_RIGHT :
#ifdef CAM_3D_VIWER_DISPLAY_KEYS
      printf("Right key pressed\n");
#endif
      break;
    case GLUT_KEY_UP :
#ifdef CAM_3D_VIWER_DISPLAY_KEYS
      printf("Up key pressed\n");
#endif
      break;
    case GLUT_KEY_DOWN :
#ifdef CAM_3D_VIWER_DISPLAY_KEYS
      printf("Down key pressed\n");
#endif
      break;
    default:
      break;
    }
}

/*************************/
/* End events handling */
/*************************/

void changeSize(int w, int h)
{
  printf("size changed\n");
  if(h == 0)
    h = 1;

  ratio = 1.0f * w / h;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glViewport(0, 0, w, h);
  gluPerspective(fov, ratio, zNear, zFar);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(posX, posY, posZ, posX + lx, posY + ly, posZ + lz, 0.0f,1.0f,0.0f);
}

void		renderSparsePoints(void)
{
  CamList	*points;

  points = pointsList;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glBegin(GL_LINES);
  // X axis
  glColor3f(1, 0, 0);
  glVertex3f(0, 0, 0);
  glVertex3f(1, 0, 0);
  // Y axis
  glColor3f(0, 1, 0);
  glVertex3f(0, 0, 0);
  glVertex3f(0, 1, 0);
  // Z axis
  glColor3f(0, 0, 1);
  glVertex3f(0, 0, 0);
  glVertex3f(0, 0, 1);
  glEnd();

  glBegin(GL_POINTS);
  while (points)
    {
      glColor3f(1, 1, 1);
      glVertex3f(((Cam3dPoint *)points->data)->x,
		 ((Cam3dPoint *)points->data)->y,
		 ((Cam3dPoint *)points->data)->z);
      points = points->next;
    }
  glEnd();
  

  glutSwapBuffers();
  return ;
}

void	cam_3d_viewer_init_display(int xPos, int yPos, int width, int height,
				   unsigned int displayMode, char *windowName,
				   void (*keyboardFuncCallback)(unsigned char, int, int),
				   void (*specialFuncCallback)(int, int, int),
				   void (*mouseWheelCallback)(int, int, int , int),
				   void	(*passiveMotionCallback)(int, int),
				   void (*displayCallback)(void),
				   void (*idleCallback)(void),
				   void (*reshapeCallback)(int, int))
{
  int	ac;

  ac = 0;
  glutInit(&ac, NULL);
  glutInitDisplayMode(displayMode);
  glutInitWindowPosition(xPos, yPos);
  glutInitWindowSize(width, height);
  glutCreateWindow(windowName);

  glClearColor(0.0,0.0,0.0,0.0); /* black background */
  gluOrtho2D(-20,20,-20,20);//(NEW)  Define our viewing area

  if (keyboardFuncCallback)
    glutKeyboardFunc(keyboardFuncCallback);
  if (specialFuncCallback)
    glutSpecialFunc(specialFuncCallback);
  glutDisplayFunc(displayCallback);
  if (idleCallback)
    glutIdleFunc(idleCallback);
  if (reshapeCallback)
    glutReshapeFunc(reshapeCallback);
  if (mouseWheelCallback)
    glutMouseWheelFunc(mouseWheelCallback);
  if (passiveMotionCallback)
    glutPassiveMotionFunc(passiveMotionCallback);
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

void	test_cam_3d_viewer()
{
  loadAyetPoints("/home/splin/manny"); // methode de chargement spécifique à chaque format de fichier
  cam_3d_viewer_init_display(100, 100, 640, 480, GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB, "Camelia Vizualizer",
			     processKeyboardKeys, processSpecialKeys, processMouseWheel, processMouseMovement,
			     renderSparsePoints, renderSparsePoints, changeSize);
  glutMainLoop();
  return ;
}
