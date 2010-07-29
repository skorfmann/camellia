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

#include <stdio.h>
#include <stdlib.h>
#include "cam_list.h"
#include "cam_3d_points.h"

POINTS_TYPE	extractNextValue1(FILE *fd)
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

CamList		*loadPoints1(char *file)
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
      ((Cam3dPoint *)pointsList->data)->x = extractNextValue1(dataFd);
      ((Cam3dPoint *)pointsList->data)->y = extractNextValue1(dataFd);
      ((Cam3dPoint *)pointsList->data)->z = extractNextValue1(dataFd);
      fgetc(dataFd);
      fgetc(dataFd);
      spacingSymbol = fgetc(dataFd);
      if (spacingSymbol == EOF)
	break;
      else
	fseek(dataFd, -1, SEEK_CUR);
    }
  fclose(dataFd);
  return (pointsList);
}

CamList		*loadPoints2(char *file)
{
  FILE		*dataFd;
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
      fread(&((Cam3dPoint *)pointsList->data)->x, sizeof(POINTS_TYPE), 1, dataFd);
      fread(&((Cam3dPoint *)pointsList->data)->y, sizeof(POINTS_TYPE), 1, dataFd);
      fread(&((Cam3dPoint *)pointsList->data)->z, sizeof(POINTS_TYPE), 1, dataFd);
      printf("here %f %f %f\n", ((Cam3dPoint *)pointsList->data)->x,
	     ((Cam3dPoint *)pointsList->data)->y,
	     ((Cam3dPoint *)pointsList->data)->z);
    }
  fclose(dataFd);
  return (pointsList);
}
