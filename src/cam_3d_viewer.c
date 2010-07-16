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
#define	POINTS_TYPE float

/* globals */
static char	*dataPath;

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

inline CamList	*cam_keypoints_3d_viewer_add_to_linked_list(CamList *l, void *data)
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

void		cam_keypoints_3d_viewer_free_data_in_linked_list(CamList *l)
{
  CamList	*ptr;

  ptr = l;
  while (ptr)
    {
      free (ptr->data);
      ptr = ptr->next;
    }
}

void	cam_keypoints_3d_viewer_disallocate_linked_list(CamList *l)
{
  cam_keypoints_tracking2_free_data_in_linked_list(l);
  cam_keypoints_tracking2_free_linked_list(l);
}


void processKeyboardKeys(unsigned char key, int x, int y)
{
#ifdef CAM_3D_VIWER_DISPLAY_KEYS
  printf("keyboard key : %i\n", (int)key);
#endif
  if (key == 27) /* echap */
    exit(0);
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

void	renderSparsePoints(void)
{
  FILE	*dataFd;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glClearColor(0.0,0.0,0.0,0.0); /* black background */

  dataFd = fopen(dataPath, "r");
  if (!dataFd)
    {
      printf("Unable to open the file %s, which is containing the data\n", dataPath);
      perror("");
      exit(-1);
    }
  
  fclose(dataFd);

  glutSwapBuffers();
  return ;
}

void	cam_3d_viewer_init_display(int xPos, int yPos, int width, int height,
				   unsigned int displayMode, char *windowName,
				   void (*keyboardFuncCallback)(unsigned char, int, int),
				   void (*specialFuncCallback)(int, int, int),
				   void (*displayCallback)(void),
				   void (*idleCallback)(void),
				   void (*reshapeCallback)(int, int),
				   char *pathToData)
{
  int	ac;

  ac = 0;
  dataPath = pathToData;
  glutInit(&ac, NULL);
  glutInitDisplayMode(displayMode);
  glutInitWindowPosition(xPos, yPos);
  glutInitWindowSize(width, height);
  glutCreateWindow(windowName);

  glEnable(GL_DEPTH_TEST);
 
  if (keyboardFuncCallback)
    glutKeyboardFunc(keyboardFuncCallback);
  if (specialFuncCallback)
    glutSpecialFunc(specialFuncCallback);
  glutDisplayFunc(displayCallback);
  glutIdleFunc(idleCallback);
  if (reshapeCallback)
    glutReshapeFunc(reshapeCallback);
}

void	test_cam_3d_viewer()
{
  cam_3d_viewer_init_display(100, 100, 640, 480, GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA, "Camelia Vizualizer",
			     processKeyboardKeys, processSpecialKeys,
			     renderSparsePoints, renderSparsePoints, NULL, "/home/splin/data");
  glutMainLoop();
  return ;
}
