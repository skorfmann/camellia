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
#define ZNEAR	1.0f
#define ZFAR	1000.0f
#define POSX	0.0f
#define	POSY	0.0f
#define	POSZ	15.0f
#define ROTX	0.0f
#define ROTY	0.0f
#define ROTZ	0.0f

/* globals */
float lastx;
float lasty;
static CamList	*pointsList = NULL;
static float	ratio;
static float	zoomFactor = 1.03f;
static float	posX=0.0f, posY=0.0f, posZ=0.0f;
static float	rotX=0.0f, rotY=0.0f, rotZ= PI;
static float	fov = 60;
static float	zNear =	1;
static float	zFar = 1000;
static BOOL	drawAxis = FALSE;

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
  /* Vect (0,0,1) */
  cam_3d_viewer_matrix_set_value(&vect, 0, 0, (POINTS_TYPE)0.0);
  cam_3d_viewer_matrix_set_value(&vect, 0, 1, (POINTS_TYPE)0.0);
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

  CamMatrixMultiply(&res, &tmp2, &vect);
  
  glLoadIdentity();
  gluLookAt(posX, posY, posZ,
	    posX + cam_3d_viewer_matrix_get_value(&res, 0, 0),
	    posY + cam_3d_viewer_matrix_get_value(&res, 0, 1),
	    posZ + cam_3d_viewer_matrix_get_value(&res, 0, 2),
	    0.0f, 1.0f, 0.0f);

  cam_3d_viewer_disallocate_matrix(&rotx);
  cam_3d_viewer_disallocate_matrix(&roty);
  cam_3d_viewer_disallocate_matrix(&rotz);
  cam_3d_viewer_disallocate_matrix(&tmp1);
  cam_3d_viewer_disallocate_matrix(&tmp2);
  cam_3d_viewer_disallocate_matrix(&vect);
  cam_3d_viewer_disallocate_matrix(&res);
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
      posZ -= 1;
      break;
    case 'z':
      posZ += 1;
      break;
    case 'q':
      posX += 1;
      break;
    case 'd':
      posX -= 1;
      break;
    case 'a':
      posY += 1;
      break;
    case 'e':
      posY -= 1;
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
   default:
      break;
    }
}

void processMouseWheel(int button, int dir, int x, int y)
{
  if (dir > 0)
    {
#ifdef CAM_3D_VIWER_DISPLAY_MOUSE
      printf("Wheel up, fov : %f\n", fov);
#endif
    }
  else
    {
#ifdef CAM_3D_VIWER_DISPLAY_MOUSE
      printf("Wheel down, fov : %f\n", fov);
#endif
    }
  glutPostRedisplay();
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
  camera();
}

// Here is the function
void glutPrint(float x, float y, char* text, float r, float g, float b)
{
  glColor3f(r,g,b);
  glRasterPos2f(x,y);
  while (*text)
    {
      glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, *text);
      text++;
    }
}  

void		renderSparsePoints(void)
{
  CamList	*points;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  points = pointsList;

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
  while (points)
    {
      glColor3f(1.0f, 1.0f, 1.0f);
      glVertex3f(((Cam3dPoint *)points->data)->x,
		 ((Cam3dPoint *)points->data)->y,
		 ((Cam3dPoint *)points->data)->z);
      points = points->next;
    }
  glEnd();
  camera();


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

  glClearColor(0.0,0.0,0.0,0.0);
  gluOrtho2D(0,width,height,0);

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
  cam_3d_viewer_init_display(100, 100, 800, 600, GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB, "Camellia Vizualizer",
			     processKeyboardKeys, processSpecialKeys, processMouseWheel, processMouseMovement,
			     renderSparsePoints, renderSparsePoints, changeSize);
  glutMainLoop();
  return ;
}
