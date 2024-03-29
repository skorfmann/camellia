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

#ifndef __CAM_PROJECT_2D_TO_3D_H__
# define __CAM_PROJECT_2D_TO_3D_H__

#include "cam_p_from_f.h"
#include "cam_2d_points.h"
#include "cam_matrix.h"

/* internals */
void		cam_compute_vector_to_3d_point(CamMatrix *v, CamMatrix *t, CamMatrix *Rt, Cam2dPoint *pt);
Cam3dPoint	*cam_vectors_intersection_perfect(CamMatrix *v1, CamMatrix *t1, CamMatrix *v2, CamMatrix *t2);
Cam3dPoint	*cam_vectors_intersection_noisy(CamMatrix *v1, CamMatrix *t1, CamMatrix *v2, CamMatrix *t2);

Cam3dPoint      *cam_triangulate_one_perfect_3d_point(CamProjectionsPair *projectionPair, CamMatrix *t1, CamMatrix *t2, CamMatrix *K, Cam2dPoint *a, Cam2dPoint *b);
Cam3dPoint      *cam_triangulate_one_noisy_3d_point(CamProjectionsPair *projectionPair, CamMatrix *t1, CamMatrix *t2, CamMatrix *K, Cam2dPoint *a, Cam2dPoint *b);

#endif /* __CAM_PROJECT_2D_TO_3D_H__*/
