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

    Redistribution and use in integral and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

        * Redistributions of integral code must retain the above copyright
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

/* CamKeypoints hessian computation
 * C code */

#ifndef INCLUDED
#ifdef CAM_FAST_APPROX_HESSIAN
int camFastApproxHessianDetectorFixedScale(CamImage *integral, CamImage *dest, int scale)
#else
int camFastHessianDetectorFixedScale(CamImage *integral, CamImage *dest, int scale)
#endif
{
    int i, c, x, y;
    int width, height;
    unsigned long *srcptr, *tmpsrcptr;
    unsigned short *dstptr, *tmpdstptr;
    CamInternalROIPolicyStruct iROI;
    int acc = 0, inc, det;
    int Dxx, Dyy, tmp;
#ifdef CAM_FAST_APPROX_HESSIAN
    int absDxx, absDyy;
#else
    int Dxy;
#endif
    int offset0, offset1, offset2, offset3, offset4, offset5, offset6, offset7, offset8, offset9;
    int offset10, offset11, offset12, offset13, offset14, offset15, offset16, offset17, offset18, offset19;
    int startx2, endx2;
    int sampling_mask;
    DECLARE_MASK_MANAGEMENT;
#ifdef ELIMINATE_EDGE_RESPONSE
    int r;
#endif
#ifndef POST_SCALING
    int coeff, multiplier, shift;
#endif
#if defined(__SSE2__)
    __m128i val0_sse2, val1_sse2, val2_sse2, val3_sse2, Dxx_sse2, Dyy_sse2, tmp_sse2, sse2_1, sse2_2;
    int *ptr, val[4];
#ifdef CAM_FAST_APPROX_HESSIAN
    __m128i cmp_sse2, cmp_sse2_2, absDxx_sse2, absDyy_sse2, absDxx_sse2_s2, absDyy_sse2_s2, det_sse2, str_sse2[2];
#else
    __m128i Dxy_sse2; 
    int _Dxx[4], _Dyy[4], _Dxy[4];
#endif
#endif 

    CAM_CHECK_ARGS2(camFastHessianDetectorFixedScale, integral->imageData != NULL, "source integral image is not allocated");
    CAM_CHECK_ARGS2(camFastHessianDetectorFixedScale, (scale >=0 && scale < sizeof(CamScale) / sizeof(CamScale[0])), "invalid scale parameter"); 
    if (dest->imageData==NULL) {
        // Automatic allocation
	CAM_CHECK(camFastHessianDetectorFixedScale, camInternalROIPolicy(integral, NULL, &iROI, 0));
        camAllocateImage(dest, iROI.srcroi.width >> CamSampling[scale], iROI.srcroi.height >> CamSampling[scale], CAM_DEPTH_16U);
    }
    CAM_CHECK(camFastHessianDetectorFixedScale, camInternalROIPolicy(integral, dest, &iROI, CAM_MASK_SUPPORT | CAM_NO_ROI_INTERSECTION));
    CAM_CHECK_ARGS(camFastHessianDetectorFixedScale, (integral->depth & CAM_DEPTH_MASK) == 32);
    CAM_CHECK_ARGS(camFastHessianDetectorFixedScale, (dest->depth & CAM_DEPTH_MASK) == 16);

    // ROI (Region Of Interest) management
    width = iROI.srcroi.width;
    height = iROI.srcroi.height;
    CAM_CHECK_ARGS2(camFastHessianDetectorFixedScale, (width & 7) == 0 && (height & 7) == 0, "ROI width and height must be multiple of 8");
    srcptr = (unsigned long *)iROI.srcptr;
    dstptr = (unsigned short *)iROI.dstptr;

    INIT_MASK_MANAGEMENT;
    
    // Algorithm initialization
    // Multiplier : for 9, divider should be 9x9 = 81.
    // In order to stay in 8 bits range, one should multiply by 1/81 = 0.0123
    // Which is equivalent to 809/65536, so we keep 809 as the multiplier for scale 9.
    // Formula is 65536/(s*s)
    // multiplier = 65536 / (CamScale[scale] * CamScale[scale]);  
    // It is made more accurate by using specifically the surfaces used for computation
#ifndef POST_SCALING
    multiplier = 65536 / ((CamScale[scale]*3 - 2 * CamOffset[scale][0]) * CamScale[scale]);  
    shift = 18 - 2 * scale;
    if (shift >= 0) 
    	coeff = (multiplier * multiplier) >> shift;
    else
	coeff = (multiplier * multiplier) << (-shift); 
#endif
    inc = (1 << CamSampling[scale]);
    //printf("scale #%d, multiplier = %d, coeff = %d\n", scale, multiplier, coeff);
    // Fill the offset tables
    offset0 = -CamScale[scale]*3 / 2 - 1 + CamOffset[scale][0];
    offset1 = CamScale[scale]*3 / 2 - CamOffset[scale][0];
    offset2 = (-CamScale[scale]*3 / 2 - 1) * iROI.srclinc;
    offset3 = offset2 + iROI.srclinc * CamScale[scale];
    offset4 = offset3 + iROI.srclinc * CamScale[scale];
    offset5 = offset4 + iROI.srclinc * CamScale[scale];
    offset6 = -CamScale[scale]*3 / 2 - 1;
    offset7 = offset6 + CamScale[scale];
    offset8 = offset7 + CamScale[scale];
    offset9 = offset8 + CamScale[scale];
    offset10 = iROI.srclinc * offset0;
    offset11 = iROI.srclinc * offset1;
    offset12 = -CamScale[scale]*3 / 2 - 1 + CamOffset[scale][1];
    offset13 = iROI.srclinc * offset12;
    offset14 = offset12 + CamScale[scale];
    offset15 = iROI.srclinc * offset14;
    offset16 = CamScale[scale]*3 / 2 - CamOffset[scale][1];
    offset17 = iROI.srclinc * offset16;
    offset18 = offset16 - CamScale[scale];
    offset19 = iROI.srclinc * offset18;
    // Other
    sampling_mask = (inc) - 1;

    // Main loop
    for (y = 0; y < height; y++) {
       	tmpsrcptr = srcptr;
	tmpdstptr = dstptr;
	BEGIN_MASK_MANAGEMENT(
	    srcptr = tmpsrcptr + startx + CamSamplingOffset[scale];
	    dstptr = tmpdstptr + (startx >> CamSampling[scale]);
	)
	    if ((y & sampling_mask) == CamSamplingOffset[scale]) {
		// Check border limits
		startx2 = startx;
		endx2 = endx;
		if ((y <= CamScale[scale]*3 / 2) || (y >= height - CamScale[scale]*3 / 2)) {
		    startx2 = endx;
		    endx2 = endx; // Skip the line
		} else {
		    if (startx <= CamScale[scale]*3 / 2) {
			startx2 = CamScale[scale]*3 / 2 + 1;
		    }
		    if (endx >= width - CamScale[scale]*3 / 2) {
			endx2 = width - CamScale[scale]*3 / 2 - 1;
		    }
		}

		for (x = startx; x < startx2 ; x += inc, srcptr += inc, dstptr++ ) *dstptr = 0;

#if defined(__SSE2__)
		if (inc == 1) {
#undef CAM_SSE2_LOAD
#define CAM_SSE2_LOAD(oLeft, oTop, oRight, oBottom) \
    val0_sse2 = _mm_loadu_si128((__m128i*)(srcptr + oRight + oBottom)); \
    val1_sse2 = _mm_loadu_si128((__m128i*)(srcptr + oLeft + oBottom)); \
    val2_sse2 = _mm_loadu_si128((__m128i*)(srcptr + oRight + oTop)); \
    val3_sse2 = _mm_loadu_si128((__m128i*)(srcptr + oLeft + oTop))
#define CAM_SSE2_INTEGRAL(i) \
    sse2_1 = _mm_sub_epi32(val0_sse2, val1_sse2); \
    sse2_2 = _mm_sub_epi32(val2_sse2, val3_sse2); \
    i = _mm_sub_epi32(sse2_1, sse2_2) 
#define INCLUDED
#include "cam_keypoints_hessian_code.c"
		} else {
#undef CAM_SSE2_LOAD
#define CAM_SSE2_LOAD(oLeft, oTop, oRight, oBottom) \
    for (c = 0, ptr = (int*)srcptr + oRight + oBottom; c != 4; c++, ptr+=inc) val[c] = *ptr; \
    val0_sse2 = _mm_set_epi32(val[3], val[2], val[1], val[0]); \
    for (c = 0, ptr = (int*)srcptr + oLeft + oBottom; c != 4; c++, ptr+=inc) val[c] = *ptr; \
    val1_sse2 = _mm_set_epi32(val[3], val[2], val[1], val[0]); \
    for (c = 0, ptr = (int*)srcptr + oRight + oTop; c != 4; c++, ptr+=inc) val[c] = *ptr; \
    val2_sse2 = _mm_set_epi32(val[3], val[2], val[1], val[0]); \
    for (c = 0, ptr = (int*)srcptr + oLeft + oTop; c != 4; c++, ptr+=inc) val[c] = *ptr; \
    val3_sse2 = _mm_set_epi32(val[3], val[2], val[1], val[0]);
#include "cam_keypoints_hessian_code.c"
		}
#undef INCLUDED
#endif

		for (; x < endx2; x += inc, srcptr += inc, dstptr++) {
		    
    #define CAM_INTEGRAL(ptr, oLeft, oTop, oRight, oBottom) \
    ( *(ptr + oRight + oBottom) - *(ptr + oLeft + oBottom) - *(ptr + oRight + oTop) + *(ptr + oLeft + oTop) )
		    
		    tmp = CAM_INTEGRAL(srcptr, offset0, offset3, offset1, offset4);
		    Dxx = CAM_INTEGRAL(srcptr, offset0, offset2, offset1, offset5) -
			(tmp + tmp + tmp);
		    Dxx >>= scale;
		    tmp = CAM_INTEGRAL(srcptr, offset7, offset10, offset8, offset11); 
		    Dyy = CAM_INTEGRAL(srcptr, offset6, offset10, offset9, offset11) -
			(tmp + tmp + tmp);
		    Dyy >>= scale;

#ifdef CAM_FAST_APPROX_HESSIAN
		    /* First basic algorithm
		    if (Dxx > 0) {
			if (Dyy > 0) {
			    if ((Dyy < (Dxx << 2)) && ((Dyy << 2) > Dxx)) {
				det = Dxx + Dyy;
			    } else det = 0;
			} else det = 0;
		    } else {
			if (Dyy < 0) {
			    if ((Dyy > (Dxx << 2)) && ((Dyy << 2) < Dxx)) {
				det = -(Dxx + Dyy);
			    } else det = 0;
			} else det = 0;
		    }			
		    */
		    // A much better implementation, that makes use of conditional execution 
		    if (Dxx > 0) { i = 1; absDxx = Dxx; } else { i = 0; absDxx = -Dxx; }
		    if (Dyy > 0) { i ^= 0; absDyy = Dyy; } else { i ^= 1; absDyy = -Dyy; }
		    if (i && (absDyy < (absDxx << 2)) && ((absDyy << 2) > absDxx))
			det = absDxx + absDyy;
		    else
			det = 0;
#else
		    det = Dxx * Dyy;
		    if (det <= 0) det = 0;
		    else {
			Dxy = CAM_INTEGRAL(srcptr, offset12, offset13, offset14, offset15) -
			    CAM_INTEGRAL(srcptr, offset18, offset13, offset16, offset15) -
			    CAM_INTEGRAL(srcptr, offset12, offset19, offset14, offset17) +
			    CAM_INTEGRAL(srcptr, offset18, offset19, offset16, offset17); 
			Dxy >>= scale;

			// Dxx, Dyy and Dxy should be 13 bits wide max
			// det = Dxx * Dyy - 0.81 * Dxy^2
			det = Dxx * Dyy - ((13 * Dxy * Dxy) >> 4);
			// det should then be 26 bits wide max
			if (det <= 0) det = 0;
			else {
#ifdef ELIMINATE_EDGE_RESPONSE
			    r = (Dxx + Dyy) * (Dxx + Dyy);
			    if (r > 10 * det)
				det = 0;
			    else 
#endif
#ifdef POST_SCALING
				// det is 16 bits max wide and can be stored in an unsigned short
				det >>= 10;
#else
				// coeff is a 7 bits wide max param
				// and thus the det should stay within 26 bits by shifting with 6 bits
				// 10 bits more and det is on 16 bits max
				det = (((unsigned int)det) * ((unsigned int)coeff)) >> 16;
#endif
			}
			acc += det;
		    }
#endif
		    *dstptr = (unsigned short)det;
		}
		for (; x < endx ; x += (inc), srcptr += inc, dstptr++ ) *dstptr = 0;
	    }
        END_MASK_MANAGEMENT;
	srcptr = (unsigned long*)(((char*)tmpsrcptr) + integral->widthStep);
	if ((y & sampling_mask) == CamSamplingOffset[scale])
	    dstptr = (unsigned short*)(((char*)tmpdstptr) + dest->widthStep);
    }

    camInternalROIPolicyExit(&iROI);
    return acc;    
}

#else // INCLUDED
		    for (i = 0; x < endx2 - 3 * inc; x += 4 * inc, srcptr += 4 * inc) {
			CAM_SSE2_LOAD(offset0, offset3, offset1, offset4);
			CAM_SSE2_INTEGRAL(tmp_sse2);
			CAM_SSE2_LOAD(offset0, offset2, offset1, offset5);
			CAM_SSE2_INTEGRAL(Dxx_sse2);
			sse2_1 = _mm_add_epi32(tmp_sse2, tmp_sse2);
			sse2_2 = _mm_sub_epi32(Dxx_sse2, tmp_sse2);
			Dxx_sse2 = _mm_sub_epi32(sse2_2, sse2_1);
			Dxx_sse2 = _mm_srai_epi32(Dxx_sse2, scale);
			CAM_SSE2_LOAD(offset7, offset10, offset8, offset11);
			CAM_SSE2_INTEGRAL(tmp_sse2);
			CAM_SSE2_LOAD(offset6, offset10, offset9, offset11);
			CAM_SSE2_INTEGRAL(Dyy_sse2);
			sse2_1 = _mm_add_epi32(tmp_sse2, tmp_sse2);
			sse2_2 = _mm_sub_epi32(Dyy_sse2, tmp_sse2);
			Dyy_sse2 = _mm_sub_epi32(sse2_2, sse2_1);
			Dyy_sse2 = _mm_srai_epi32(Dyy_sse2, scale);

#ifdef CAM_FAST_APPROX_HESSIAN
			// Compute absolute value of Dxx
			sse2_1 = _mm_sub_epi32(_mm_setzero_si128(), Dxx_sse2);
			cmp_sse2 = _mm_cmplt_epi32(Dxx_sse2, _mm_setzero_si128());
			sse2_1 = _mm_and_si128(cmp_sse2, sse2_1);
			sse2_2 = _mm_andnot_si128(cmp_sse2, Dxx_sse2);
			absDxx_sse2 = _mm_add_epi32(sse2_1, sse2_2);
			
			// Compute absolute value of Dyy
			sse2_1 = _mm_sub_epi32(_mm_setzero_si128(), Dyy_sse2);
			cmp_sse2_2 = _mm_cmplt_epi32(Dyy_sse2, _mm_setzero_si128());
			sse2_1 = _mm_and_si128(cmp_sse2_2, sse2_1);
			sse2_2 = _mm_andnot_si128(cmp_sse2_2, Dyy_sse2);
			absDyy_sse2 = _mm_add_epi32(sse2_1, sse2_2);
			
			cmp_sse2 = _mm_xor_si128(cmp_sse2, cmp_sse2_2);
			det_sse2 = _mm_add_epi32(absDxx_sse2, absDyy_sse2);
			det_sse2 = _mm_andnot_si128(cmp_sse2, det_sse2);
			absDxx_sse2_s2 = _mm_slli_epi32(absDxx_sse2, 2);
			absDyy_sse2_s2 = _mm_slli_epi32(absDyy_sse2, 2);
			sse2_1 = _mm_cmplt_epi32(absDyy_sse2, absDxx_sse2_s2);
			sse2_2 = _mm_cmpgt_epi32(absDyy_sse2_s2, absDxx_sse2);
			det_sse2 = _mm_and_si128(det_sse2, sse2_1);
			str_sse2[i] = _mm_and_si128(det_sse2, sse2_2);

			i++;
			if (i == 2) {
			    i = 0;
			    det_sse2 = _mm_packs_epi32(str_sse2[0], str_sse2[1]);
			    _mm_storeu_si128((__m128i*)dstptr, det_sse2);
			    dstptr += 8;
			}
		    }
		    if (i == 1) {
			_mm_storeu_si128((__m128i*)val, str_sse2[0]);
			for (i = 0; i != 4; i++, dstptr++) *dstptr = val[i];
		    }
#else
			
			CAM_SSE2_LOAD(offset12, offset13, offset14, offset15);
			CAM_SSE2_INTEGRAL(Dxy_sse2);
			CAM_SSE2_LOAD(offset18, offset13, offset16, offset15);
			CAM_SSE2_INTEGRAL(sse2_1);
			Dxy_sse2 = _mm_sub_epi32(Dxy_sse2, sse2_1);
			CAM_SSE2_LOAD(offset12, offset19, offset14, offset17);
			CAM_SSE2_INTEGRAL(sse2_1);
			Dxy_sse2 = _mm_sub_epi32(Dxy_sse2, sse2_1);
			CAM_SSE2_LOAD(offset18, offset19, offset16, offset17); 
			CAM_SSE2_INTEGRAL(sse2_1);
			Dxy_sse2 = _mm_add_epi32(Dxy_sse2, sse2_1);
			Dxy_sse2 = _mm_srai_epi32(Dxy_sse2, scale);
		    
			_mm_storeu_si128((__m128i*)_Dxy, Dxy_sse2);
			_mm_storeu_si128((__m128i*)_Dxx, Dxx_sse2);
			_mm_storeu_si128((__m128i*)_Dyy, Dyy_sse2);

			for (i = 0; i != 4; i++, dstptr++) {
			    // _Dxx[i], _Dyy[i] and _Dxy[i] should be 12 to 13 bits wide max
			    // det = _Dxx[i] * _Dyy[i] - 0.81 * _Dxy[i]^2
			    det = _Dxx[i] * _Dyy[i] - ((13 * _Dxy[i] * _Dxy[i]) >> 4);
			    // det should then be 26 bits wide max
			    if (det <= 0) det = 0;
			    else {
#ifdef ELIMINATE_EDGE_RESPONSE
				r = (_Dxx[i] + _Dyy[i]) * (_Dxx[i] + _Dyy[i]);
				if (r > 10 * det)
				    det = 0;
				else 
#endif
#ifdef POST_SCALING
				// det is 16 bits max wide and can be stored in an unsigned short
				det >>= 10;
#else
				// coeff is a 7 bits wide max param
				// and thus the det should stay within 26 bits by shifting with 6 bits
				// 10 bits more and det is on 16 bits max
				det = (((unsigned int)det) * ((unsigned int)coeff)) >> 16;
#endif
			    }
			    acc += det;
			    *dstptr = (unsigned short)det;
			}
		    }
#endif
#endif // INCLUDED

