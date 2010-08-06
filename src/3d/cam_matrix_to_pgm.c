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
#include <string.h>
#include "cam_matrix_to_pgm.h"
#include "cam_matrix.h"

void	cam_matrix_to_pgm(char *filename, CamImageMatrix *m)
{
  int	i;
  int	j;
  FILE	*file;
  char	header[20];
  
  file = fopen(filename, "w+");
  if (!file)
    {
      printf("cam_write_points_to_pgm : unable to open the destination image file\n");
      exit (-1);
    }
  sprintf(header, "P6\n%i %i\n255\n", m->r.ncols, m->r.nrows);
  fwrite(header, sizeof(char), strlen(header), file);
  for (j = 0 ; j < m->r.nrows ; ++j)
    {
      for (i = 0 ; i < m->r.ncols ; ++i)
	{
	  fprintf(file, "%c%c%c",
		  (char)cam_matrix_get_value(&m->r, i, j),
		  (char)cam_matrix_get_value(&m->g, i, j),
		  (char)cam_matrix_get_value(&m->b, i, j));
	}
    }
  fclose(file);
}

void	cam_matrix_to_pathpgm(char *dir, char *fileName, CamImageMatrix *m)
{
  char	*outputPath;

  outputPath = (char *)malloc((strlen(dir) + strlen(fileName) + strlen("pgm") + 3) * sizeof(char) );
  sprintf(outputPath,"%s/%s.%s", dir, fileName, "pgm");
  cam_matrix_to_pgm(outputPath, m);
  free(outputPath);
}
