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

    Copyright (c) 2002-2008, Ecole des Mines de Paris - Centre de Robotique
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

/* RLE Utilties
 * C code */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "camellia.h"
#include "camellia_internals.h"

#define parent blob

// Very useful to post-process a RLE image
// Applies a LUT (i.e. for instance thresholding), so that we can join runs
int camRLEApplyLUT(CamRLEImage *src, CamRLEImage *dest, CamLUT *LUT)
{
    CamRun *in, *out;
    int i, j=1, x, currentColor, currentLength, currentPos, value, line = 0;
    int num = src->nbRuns;

    CAM_CHECK_ARGS(camRLEApplyLUT, src != NULL);
    CAM_CHECK_ARGS(camRLEApplyLUT, dest != NULL);
    CAM_CHECK_ARGS(camRLEApplyLUT, LUT != NULL);
    // Automatic allocation
    if (dest->allocated == 0) camRLEAllocate(dest, src->nbRuns);
    CAM_CHECK_ARGS(camRLEApplyLUT, (dest->nSize == sizeof(CamRLEImage)));
    CAM_CHECK_ARGS(camRLEApplyLUT, src != dest);
    CAM_CHECK_ARGS(camRLEApplyLUT, src->options & CAM_RLEOPTS_ZERO_ENCODED);

    dest->options = CAM_RLEOPTS_ZERO_ENCODED | CAM_RLEOPTS_LINES_ENCODED;

    // Put a run at the beginning to have a reference starting point
    out = &dest->runs[0];
    out->value = 0;
    out->length = 0;
    out->parent = -1;
    out->x = 0;
 
    currentColor = LUT->t[src->runs[1].value];
    currentLength = x = src->runs[1].length;
    currentPos = 0;
    for (i = 2; i < num; i++) {
	in = &src->runs[i];
	if (in->length == 0) continue;
	value = LUT->t[in->value];
	if (x != src->width) {
	    x += in->length;
	    if (value == currentColor) {
		currentLength += in->length;
	    } else {
		out = &dest->runs[j];
		out->value = currentColor;
		out->length = currentLength;
		out->parent = -1;
		out->x = currentPos;
		j++;
		currentColor = value;
		currentPos += currentLength;
		currentLength = in->length;
	    }
	} else {
	    // Write last run of current line
	    out = &dest->runs[j];
	    out->value = currentColor;
	    out->length = currentLength;
	    out->parent = j;
            out->x = currentPos;
	    j++;
	    currentColor = value;
	    currentPos = 0;
	    currentLength = x = in->length;

	    // Write delimitating run
	    out = &dest->runs[j];
	    out->value = 0;
	    out->length = 0;
	    out->parent = j;
            out->x = ++line;
	    j++;
	}
    }

    // Write last run of current line
    out = &dest->runs[j];
    out->value = currentColor;
    out->length = currentLength;
    out->parent = j;
    out->x = currentPos;
    j++;
    
    // Write delimitating run
    out = &dest->runs[j];
    out->value = 0;
    out->length = 0;
    out->parent = j;
    out->x = ++line;
    j++;
    
    dest->nbRuns=j;
    dest->height=src->height;
    dest->width=src->width;

    assert(dest->nbRuns <= dest->allocated);
    return 1;
}

#define CAM_PIXEL unsigned char

int camRLEDecode(CamRLEImage *src, CamImage *dest, CamLUT *LUT)
{
    CAM_PIXEL *ptr;
    CamRun *in=src->runs+1;
    int i,j,k,l,m;
    int nbRuns=src->nbRuns;
    int width, height, value;
    CamInternalROIPolicyStruct iROI;

    CAM_CHECK_ARGS(camRLEDecode, src != NULL);
    CAM_CHECK_ARGS(camRLEDecode, dest != NULL);
    if ((src->runs == NULL)||(src->nbRuns == 0)) return 0;
    // Automatic allocation
    if (dest->imageData==NULL) {
        if (!camAllocateImage(dest,src->width,src->height,CAM_DEPTH_8U)) return 0;
    }

    // ROI (Region Of Interest) management
    CAM_CHECK(camRLEDecode, camInternalROIPolicy(dest, NULL, &iROI, 0));
    CAM_CHECK_ARGS(camRLEDecode, ((dest->depth&CAM_DEPTH_MASK)==(sizeof(CAM_PIXEL)*8)));
    width=iROI.srcroi.width;
    height=iROI.srcroi.height;
    if (iROI.nChannels!=1) {
        CAM_CHECK_ARGS(camRLEDecode, (dest->dataOrder==CAM_DATA_ORDER_PIXEL));
    }
    CAM_CHECK_ARGS(camRLEDecode, src->width==width);
    CAM_CHECK_ARGS(camRLEDecode, src->height==height);

    ptr=(CAM_PIXEL*)iROI.srcptr;
    for (i=1,l=0,k=0;i<nbRuns;i++,in++) {
	if (LUT) {
	    if (in->value<CAM_TABLE_SIZE) {
		value=LUT->t[in->value];
	    } else value=0;
	} else {
	    value=in->value;
	}
        for (j=0;j<in->length;j++) {
            for (m=0;m<iROI.nChannels;m++) {
                ptr[k+m]=(value>>(m<<3));
            }
            k+=iROI.srcinc;
        }
        l+=in->length;
        if (l==src->width) {
            l=0;
            k=0;
            ptr+=iROI.srclinc;
        }
    }

    return 1;
}

// TO DO : update this function according to camRLEDecode function
int camRLEDecodeBlobs(CamRLEImage *src, CamImage *dest, CamLUT *LUT)
{
    CAM_PIXEL *ptr;
    CamRun *in=src->runs+1;
    int i,j,k;
    int nbRuns=src->nbRuns;

    CAM_CHECK_ARGS(camRLEDecodeBlobs, ((dest->nChannels==1)&&((dest->depth&CAM_DEPTH_MASK)==(sizeof(CAM_PIXEL)*8))));

    if (dest->roi==NULL) {
	CAM_CHECK_ARGS(camRLEDecodeBlobs, dest->width==dest->widthStep);
	CAM_CHECK_ARGS(camRLEDecodeBlobs, src->width==dest->width);
	for (i=1,k=0;i<nbRuns;i++,in++) {
	    for (j=0;j<in->length;j++,k++) {
		if (in->value) {
		    dest->imageData[k]=LUT->t[in->parent];
		} else {
		    dest->imageData[k]=0;
		}
	    }
	}
    } else {
	ptr=(CAM_PIXEL*)(dest->imageData+dest->roi->yOffset*dest->widthStep)+dest->roi->xOffset;
	for (i=1,k=0;i<nbRuns;i++,in++) {
	    for (j=0;j<in->length;j++,k++) {
		if (in->value) {
		    *(ptr++)=LUT->t[in->parent];
		} else {
		    *(ptr++)=0;
		}
	    }
	    if (k==src->width) {
		ptr+=(dest->widthStep-src->width)/sizeof(CAM_PIXEL);
		k=0;
	    }
	}
    }
    return 1;
}

// Very useful for finding holes in blobs
int camRLEInverse(CamRLEImage *image)
{
    CamRun *in=image->runs+1;
    int nbRuns=image->nbRuns;
    int i;

    for (i=1;i<nbRuns;i++,in++) {
	if (in->length) {
	    in->value=(in->value)?0:1;
	    in->parent=i;
	}
    }
    return 1;
}

int camRLEAllocate(CamRLEImage *rle, int max_runs)
{
    rle->nSize=sizeof(CamRLEImage);
    rle->id=0;
    rle->runs=(CamRun*)malloc(sizeof(CamRun)*max_runs);
    if (rle->runs == NULL) {
        rle->allocated = 0;
	camError("camRLEAallocate","Memory allocation error");
        return 0;
    }
    rle->allocated=max_runs;
    rle->nbRuns=0;
    return 1;
}

int camRLEDeallocate(CamRLEImage *rle)
{
    if (rle->runs) free(rle->runs);
    rle->runs=NULL;
    rle->allocated=0;
    return 1;
}

int camRLEFree(CamRLEImage *rle)
{
    return camRLEDeallocate(rle);
}

int camRLEReallocate(CamRLEImage *rle, int new_max_runs)
{
    if (rle->runs==NULL) {
        return camRLEAllocate(rle,new_max_runs);
    } 
    rle->runs=(CamRun*)realloc(rle->runs,new_max_runs*sizeof(CamRun));
    if (rle->runs==NULL) {
        rle->nbRuns=0;
        rle->allocated=0;
        camError("camRLEReallocate","Memory allocation error");
        return 0;
    }
    rle->allocated=new_max_runs;
    return 1;
}

int camRLEClone(CamRLEImage *source, CamRLEImage *dest)
{
    *dest=*source;
    if (source->runs!=NULL) {
        dest->runs=(CamRun*)malloc(sizeof(CamRun)*source->allocated);
        if (dest->runs==0) {
            camError("camRLEClone","Memory allocation error");
            return 0;
        }
        memcpy(dest->runs,source->runs,source->nbRuns);
    }
    return 1;
}

int camRLEBlobSides(CamBlobInfo *blob, int *left, int *top, int *right, int *bottom)
{
    return 0;
}

int camRLEBlobROIIntersect(CamBlobInfo *blob, CamROI *roi)
{
    return 0;
}

int camRLEBlobMeasures(CamBlobInfo *blob, CamImage *original)
{
    return 0;
}


