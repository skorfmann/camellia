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

    Copyright (c) 2002-2006, Ecole des Mines de Paris - Centre de Robotique
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

/* Median Filtering Kernel
 * C code */

#include <stdio.h>
#include <string.h>
#include "camellia.h"
#include "camellia_internals.h"

#ifdef CAM_GENERATE_FULL_CODE

// 8 and 16 bits pixel size code generation

#undef CAM_PIXEL
#define CAM_PIXEL unsigned char
// 3x3 Kernel Median Filtering code generation
#define camMedianFilter camMedianFilter3x38
#define camMedianFilterParams camMedianFilter3x3Params
#define CAM_MF_NEIGHB 3
#include "cam_median_filtering_code.c"
#undef camMedianFilter
#undef camMedianFilterParams
#undef CAM_MF_NEIGHB

// 5x5 Kernel Median Filtering code generation
#define camMedianFilter camMedianFilter5x58
#define camMedianFilterParams camMedianFilter5x5Params
#define CAM_MF_NEIGHB 5
#include "cam_median_filtering_code.c"
#undef camMedianFilter
#undef camMedianFilterParams
#undef CAM_MF_NEIGHB

#undef CAM_PIXEL
#define CAM_PIXEL unsigned short
// 3x3 Kernel Median Filtering code generation
#define camMedianFilter camMedianFilter3x316
#define camMedianFilterParams camMedianFilter3x3Params
#define CAM_MF_NEIGHB 3
#include "cam_median_filtering_code.c"
#undef camMedianFilter
#undef camMedianFilterParams
#undef CAM_MF_NEIGHB

// 5x5 Kernel Median Filtering code generation
#define camMedianFilter camMedianFilter5x516
#define camMedianFilterParams camMedianFilter5x5Params
#define CAM_MF_NEIGHB 5
#include "cam_median_filtering_code.c"
#undef camMedianFilter
#undef camMedianFilterParams
#undef CAM_MF_NEIGHB

int camMedianFilter3x3(CamImage *source, CamImage *dest)
{
    if ((source->depth&CAM_DEPTH_MASK)>8) {
        return camMedianFilter3x316(source,dest);
    } else {
        return camMedianFilter3x38(source,dest);
    }
}

int camMedianFilter5x5(CamImage *source, CamImage *dest)
{
    if ((source->depth&CAM_DEPTH_MASK)>8) {
        return camMedianFilter5x516(source,dest);
    } else {
        return camMedianFilter5x58(source,dest);
    }
}
#else

// 3x3 Kernel Median Filtering code generation
#define camMedianFilter camMedianFilter3x3
#define camMedianFilterParams camMedianFilter3x3Params
#define CAM_MF_NEIGHB 3
#include "cam_median_filtering_code.c"
#undef camMedianFilter
#undef camMedianFilterParams
#undef CAM_MF_NEIGHB

// 5x5 Kernel Median Filtering code generation
#define camMedianFilter camMedianFilter5x5
#define camMedianFilterParams camMedianFilter5x5Params
#define CAM_MF_NEIGHB 5
#include "cam_median_filtering_code.c"
#undef camMedianFilter
#undef camMedianFilterParams
#undef CAM_MF_NEIGHB

#endif

