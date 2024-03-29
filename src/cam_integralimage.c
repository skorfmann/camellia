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

#ifndef CAM_INLINED
#define CAM_INLINED
/* Histogram computation
 * C code */
#include <stdlib.h>
#include "camellia.h"
#include "camellia_internals.h"

#ifdef CAM_VECTORIZE
typedef CAM_UINT32 v4si __attribute__ ((vector_size(16)));
union i4vector 
{
  v4si v;
  CAM_UINT32 i[4];
};
#endif

#ifdef CAM_GENERATE_FULL_CODE

// 8 and 16 bits pixel size code generation

#undef CAM_PIXEL
#define CAM_PIXEL unsigned char
#define camIntegralImage camIntegralImage8
#include "cam_integralimage.c"
#undef camIntegralImage
  
#undef CAM_PIXEL
#define CAM_PIXEL unsigned short
#define camIntegralImage camIntegralImage16
#include "cam_integralimage.c"
#undef camIntegralImage
  
int camIntegralImage(CamImage *src, CamImage *dest)
{
  if ((src->depth&CAM_DEPTH_MASK)>8) {
	return camIntegralImage16(src, dest);
    } else {
	return camIntegralImage8(src, dest);
    }
}

#else
#include "cam_integralimage.c"
#endif

#else

int camIntegralImage(CamImage *src, CamImage *dest)
{
    int x, y, i;
    int width, height;
    CAM_PIXEL *srcptr, *srctmpptr;
    CAM_UINT32 *dstptr, *dsttmpptr;
    CamInternalROIPolicyStruct iROI;
    CAM_UINT32 val;

    CAM_CHECK_ARGS2(camIntegralImage, src->imageData != NULL, "source image is not allocated");
    if (dest->imageData==NULL) {
        // Automatic allocation
	CAM_CHECK(camIntegralImage, camInternalROIPolicy(src, NULL, &iROI, 0));
        camAllocateImageEx(dest, iROI.srcroi.width, iROI.srcroi.height, CAM_DEPTH_32U, (iROI.nChannels == 1)?CAM_CHANNELSEQ_GREY:CAM_CHANNELSEQ_YUV);
	CAM_CHECK(camIntegralImage, camInternalROIPolicy(src, dest, &iROI, 0));
    }
    CAM_CHECK(camIntegralImage, camInternalROIPolicy(src, dest, &iROI, 0));
    CAM_CHECK_ARGS(camIntegralImage, (src->depth & CAM_DEPTH_MASK) >= 8);
    CAM_CHECK_ARGS(camIntegralImage, (src->depth & CAM_DEPTH_MASK) <= (sizeof(CAM_PIXEL) * 8));
    CAM_CHECK_ARGS(camIntegralImage, (dest->depth & CAM_DEPTH_MASK) == 32);
    CAM_CHECK_ARGS(camIntegralImage, !(src->depth & CAM_DEPTH_SIGN)); 

/*
    // One pass algorithm
    width = iROI.srcroi.width;
    height = iROI.srcroi.height;
    srcptr = (CAM_PIXEL*)iROI.srcptr;
    dstptr = (CAM_UINT32*)iROI.dstptr;
    srctmpptr = srcptr;
    dsttmpptr = dstptr;
    val = 0;	
    for (x = 0; x < width; x++, srcptr += iROI.srcinc, dstptr++) {
	val += *srcptr;
	*dstptr = val; 
    }
    srcptr = (CAM_PIXEL*)(((char*)srctmpptr) + src->widthStep);
    dstptr = (CAM_UINT32*)(((char*)dsttmpptr) + dest->widthStep);
    for (y = 1; y < height; y++) {
	srctmpptr = srcptr;
	dsttmpptr = dstptr;
	val = 0;	
	for (x = 0; x < width; x++, srcptr += iROI.srcinc, dstptr++) {
	    val += *srcptr;
	    *dstptr = val + *(CAM_UINT32*)(((char*)dstptr) - dest->widthStep); 
	}
	srcptr = (CAM_PIXEL*)(((char*)srctmpptr) + src->widthStep);
	dstptr = (CAM_UINT32*)(((char*)dsttmpptr) + dest->widthStep);
    }
*/

    for (i = 0; i < iROI.nChannels; i++) {
	// Two passes algorithm
	// First pass 
	width = iROI.srcroi.width;
	height = iROI.srcroi.height;
	srcptr = ((CAM_PIXEL*)iROI.srcptr) + i * iROI.srcpinc;
	dstptr = ((CAM_UINT32*)iROI.dstptr) + i * iROI.dstpinc;
	for (y = 0; y < height; y++) {
	    srctmpptr = srcptr;
	    dsttmpptr = dstptr;
	    val = 0;	
	    for (x = 0; x < width; x++, srcptr += iROI.srcinc, dstptr++) {
		val += *srcptr;
		*dstptr = val; 
	    }
	    srcptr = (CAM_PIXEL*)(((char*)srctmpptr) + src->widthStep);
	    dstptr = (CAM_UINT32*)(((char*)dsttmpptr) + dest->widthStep);
	}

	// Second pass 
	dstptr = ((CAM_UINT32*)iROI.dstptr) + i * iROI.dstpinc;
	dstptr = (CAM_UINT32*)(((char*)dstptr) + dest->widthStep);
	for (y = 1; y < height; y++) {
	    dsttmpptr = dstptr;
	    for (x = 0; x < width; x++, dstptr++) {
		*dstptr += *(CAM_UINT32*)(((char*)dstptr) - dest->widthStep); 
	    }
	    dstptr = (CAM_UINT32*)(((char*)dsttmpptr) + dest->widthStep);
	}
    }
    camInternalROIPolicyExit(&iROI);
    return 1;    
}
#endif
