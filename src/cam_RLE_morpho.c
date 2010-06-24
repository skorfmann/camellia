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

/* RLE morpho/erosion code
 * C code */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "camellia.h"
#include "camellia_internals.h"

#define parent blob

int camRLEErodeCross(CamRLEImage *source, CamRLEImage *dest)
{
    return 0;
}

int camRLEErode3x3(CamRLEImage *source, CamRLEImage *dest)
{
    return 0;
}

int camRLEErode3x2(CamRLEImage *source, CamRLEImage *dest)
{
    return 0;
}

#define MAX_RUNS 100000

static int camCompareRUNS (void const *a, void const *b)
{
	CamRun const *pa = (CamRun *) a;
	CamRun const *pb = (CamRun *) b;
	if(pa->blob < pb->blob) return -1;
	else if(pa->blob > pb->blob) return 1;
	else{ 
		if(pa->x < pb->x) return -1;
		else if(pa->x > pb->x) return 1;
		else return 0;
	}
}

int camRLEDilate(CamRLEImage *source, CamRLEImage *dest, CamRLEImage *strElement)
{
	CamRun primaryRuns[MAX_RUNS];
	CamRun finalRuns[MAX_RUNS];
	CamRun temp, run, last;

	int i=0, j, m, n, y_strElement = 0, y_source = 0, ind = 0, cpt_runs = 0, y_finalRuns = 0, y_indice = 0, dest_indice = 0, add_run = 0;

	//camRLEClone(source, dest);
	for(m=0;m<strElement->nbRuns;m++){
		printf("strElement Run #%d, x %d, length %d, value %d\n", m, strElement->runs[m].x, strElement->runs[m].length, strElement->runs[m].value);
		if((strElement->runs[m].value != 0)&&(strElement->runs[m].length != 0)) cpt_runs++;
	}
	printf("strElement ==> Number of Runs : %d\n", cpt_runs);
	cpt_runs = 0;
	for(n=0;n<source->nbRuns;n++){
		printf("source Run #%d, x %d, length %d, value %d\n", n, source->runs[n].x, source->runs[n].length, source->runs[n].value);
		if((source->runs[n].value != 0)&&(source->runs[n].length != 0)) cpt_runs++;
	}
	printf("source ==> Number of Runs : %d\n", cpt_runs);
	for(m=0;m<strElement->nbRuns;m++){
		for(n=0;n<source->nbRuns;n++){
			if(source->runs[n].length == 0) y_source = source->runs[n].x;
			else if(strElement->runs[m].length == 0) y_strElement = strElement->runs[m].x;
			else{
				if((source->runs[n].value != 0) && (strElement->runs[m].value != 0)){
					temp.blob = y_strElement + y_source;
					temp.x = strElement->runs[m].x + source->runs[n].x;
					if(temp.x > source->width) temp.x = source->width;
					if(temp.x < 0) temp.x = 0;
					temp.length = strElement->runs[m].length + source->runs[n].length;
					if((temp.x + temp.length) > source->width) temp.length = source->width - temp.x;
					temp.value = source->runs[n].value*255;
					primaryRuns[ind] = temp;
					ind++;
				}
			}
		}
	}
	qsort(primaryRuns, ind, sizeof(temp), camCompareRUNS);
	//y_indice = primaryRuns[0].blob;
	temp = primaryRuns[0];
	while(i<ind){
		y_indice = temp.blob;
		add_run = 1;
		while(y_indice == primaryRuns[i].blob){
			printf("primaryRuns[%d] => xs : %d, xe : %d, y : %d, value : %d\n", i, primaryRuns[i].x, primaryRuns[i].x + primaryRuns[i].length, primaryRuns[i].blob, primaryRuns[i].value);
			printf("temp => xs : %d, xe : %d, y : %d, value : %d\n", temp.x, temp.x + temp.length, temp.blob, temp.value);
			//i++;
			if(primaryRuns[i].x > (temp.x + temp.length)){
				//c'est un nouveau run, passer a un autre run et enregistrer le temp dans le tableau final
				printf("got a new run!! ");
				printf("added #%d xs : %d, xe : %d, y : %d\n", y_finalRuns, temp.x, temp.x + temp.length, y_indice);
				finalRuns[y_finalRuns] = temp;
				y_finalRuns++;
				temp = primaryRuns[i];
				i++;
				add_run = 0;
			}
			else{
				//le debut du run qu'on teste est avant le fin du temp, il y a deux cas
				if((primaryRuns[i].x + primaryRuns[i].length) < (temp.x + temp.length)){
					//le run qu'on teste est inclus dans le temp donc on le neglige et on avance
					printf("Run ignored\n");
					i++;
					add_run = 1;
				}
				else{
					//le run qu'on teste a une partie commune avec le temp, il faut les fusionner et remplacer le temp par le nouveau run resultant
					temp.length = primaryRuns[i].x + primaryRuns[i].length - temp.x;
					printf("Run expansion. New Run : xs : %d, xe : %d\n", temp.x, temp.x + temp.length);
					i++;
					add_run = 1;
				}
			}
		}
		if(add_run){
			printf("New run added #%d xs : %d, xe : %d\n", y_finalRuns, temp.x, temp.x + temp.length);
			finalRuns[y_finalRuns] = temp;
			y_finalRuns++;
			temp = primaryRuns[i];
		}
		//printf("current indice : %d, i : %d\n", y_indice, i);
	}
	printf("final runs : %d\n", y_finalRuns);
	for(i=0;i<y_finalRuns; i++){
		//printf("finalRun #%d xs : %d, xe : %d, y : %d, value : %d\n", i, finalRuns[i].x , finalRuns[i].x + finalRuns[i].length, finalRuns[i].blob, finalRuns[i].value);
	}
	printf("--------------------------\n");
	//creer l'image dilatee
	last.x = 0;
	last.blob = 0;
	last.length = 0;
	last.value = 0;
	dest->runs[dest_indice] = last;
	printf("dest #%d xs : %d, length : %d, end : %d, y : %d, value : %d\n", dest_indice, dest->runs[dest_indice].x , dest->runs[dest_indice].length, dest->runs[dest_indice].x + dest->runs[dest_indice].length, dest->runs[dest_indice].blob, dest->runs[dest_indice].value);
	dest_indice++;
	for(i=0; i<y_finalRuns; i++){
		run = finalRuns[i];
		if(run.blob != last.blob){
			temp.x = last.x + last.length;
			temp.length = source->width - (last.x + last.length);
			temp.blob = last.blob;
			temp.value = 0;
			if(temp.x < source->width){
				dest->runs[dest_indice] = temp;
				printf("dest #%d xs : %d, length : %d, end : %d, y : %d, value : %d\n", dest_indice, dest->runs[dest_indice].x , dest->runs[dest_indice].length, dest->runs[dest_indice].x + dest->runs[dest_indice].length, dest->runs[dest_indice].blob, dest->runs[dest_indice].value);
				dest_indice++;
			}
			for(j=1;j<(run.blob - last.blob); j++){
				temp.x = j;
				temp.length = 0;
				temp.blob = j;
				temp.value = 0;
				dest->runs[dest_indice] = temp;
				printf("dest #%d xs : %d, length : %d, end : %d, y : %d, value : %d\n", dest_indice, dest->runs[dest_indice].x , dest->runs[dest_indice].length, dest->runs[dest_indice].x + dest->runs[dest_indice].length, dest->runs[dest_indice].blob, dest->runs[dest_indice].value);
				dest_indice++;
				temp.x = 0;
				temp.length = source->width;
				temp.blob = j;
				temp.value = 0;
				dest->runs[dest_indice] = temp;
				printf("dest #%d xs : %d, length : %d, end : %d, y : %d, value : %d\n", dest_indice, dest->runs[dest_indice].x , dest->runs[dest_indice].length, dest->runs[dest_indice].x + dest->runs[dest_indice].length, dest->runs[dest_indice].blob, dest->runs[dest_indice].value);
				dest_indice++;
			}
			temp.x = run.blob;
			temp.length = 0;
			temp.blob = run.blob;
			temp.value = 0;
			dest->runs[dest_indice] = temp;
			printf("dest #%d xs : %d, length : %d, end : %d, y : %d, value : %d\n", dest_indice, dest->runs[dest_indice].x , dest->runs[dest_indice].length, dest->runs[dest_indice].x + dest->runs[dest_indice].length, dest->runs[dest_indice].blob, dest->runs[dest_indice].value);
			dest_indice++;
			if(run.x != 0){
				temp.x = 0;
				temp.length = run.x;
				temp.blob = run.blob;
				temp.value = 0;
				dest->runs[dest_indice] = temp;
				printf("dest #%d xs : %d, length : %d, end : %d, y : %d, value : %d\n", dest_indice, dest->runs[dest_indice].x , dest->runs[dest_indice].length, dest->runs[dest_indice].x + dest->runs[dest_indice].length, dest->runs[dest_indice].blob, dest->runs[dest_indice].value);
				dest_indice++;
			}
			dest->runs[dest_indice] = run;
			printf("dest #%d xs : %d, length : %d, end : %d, y : %d, value : %d\n", dest_indice, dest->runs[dest_indice].x , dest->runs[dest_indice].length, dest->runs[dest_indice].x + dest->runs[dest_indice].length, dest->runs[dest_indice].blob, dest->runs[dest_indice].value);
			dest_indice++;
			last = run;
		}
		else{
			temp.x = last.x + last.length;
			temp.length = run.x - (last.x + last.length);
			temp.blob = run.blob;
			temp.value = 0;
			dest->runs[dest_indice] = temp;
			printf("dest #%d xs : %d, length : %d, end : %d, y : %d, value : %d\n", dest_indice, dest->runs[dest_indice].x , dest->runs[dest_indice].length, dest->runs[dest_indice].x + dest->runs[dest_indice].length, dest->runs[dest_indice].blob, dest->runs[dest_indice].value);
			dest_indice++;
			dest->runs[dest_indice] = run;
			printf("dest #%d xs : %d, length : %d, end : %d, y : %d, value : %d\n", dest_indice, dest->runs[dest_indice].x , dest->runs[dest_indice].length, dest->runs[dest_indice].x + dest->runs[dest_indice].length, dest->runs[dest_indice].blob, dest->runs[dest_indice].value);
			dest_indice++;
			last = run;
		}
	}
	if((dest->runs[dest_indice-1].x + dest->runs[dest_indice-1].length) < source->width){
		temp.x = dest->runs[dest_indice-1].x + dest->runs[dest_indice-1].length;
		temp.length = source->width - dest->runs[dest_indice-1].x + dest->runs[dest_indice-1].length;
		temp.blob = dest->runs[dest_indice-1].blob;
		temp.value = 0;
		dest->runs[dest_indice] = temp;
		printf("dest #%d xs : %d, length : %d, end : %d, y : %d, value : %d\n", dest_indice, dest->runs[dest_indice].x , dest->runs[dest_indice].length, dest->runs[dest_indice].x + dest->runs[dest_indice].length, dest->runs[dest_indice].blob, dest->runs[dest_indice].value);
		dest_indice++;
	}
	/*
	for(i=0;i<dest_indice; i++){
		printf("dest #%d xs : %d, length : %d, y : %d, value : %d\n", i, dest->runs[i].x , dest->runs[i].length, dest->runs[i].blob, dest->runs[i].value);
	}*/

	if(dest->runs[dest_indice-1].blob < source->height){
		printf("adding needed lines\n");
		printf("dest->runs[dest_indice].blob : %d, source->height : %d\n", dest->runs[dest_indice - 1].blob, source->height);
		for(j=dest->runs[dest_indice-1].blob; j<source->height-1;j++){
				temp.x = j;
				temp.length = 0;
				temp.blob = j;
				temp.value = 0;
				dest->runs[dest_indice] = temp;
				printf("destRuns #%d xs : %d, length : %d, y : %d, value : %d\n", dest_indice, dest->runs[dest_indice].x , dest->runs[dest_indice].length, dest->runs[dest_indice].blob, dest->runs[dest_indice].value);
				dest_indice++;
				temp.x = 0;
				temp.length = source->width;
				temp.blob = j;
				temp.value = 0;
				dest->runs[dest_indice] = temp;
				printf("destRuns #%d xs : %d, length : %d, y : %d, value : %d\n", dest_indice, dest->runs[dest_indice].x ,dest->runs[dest_indice].length, dest->runs[dest_indice].blob, dest->runs[dest_indice].value);
				dest_indice++;

		}
		temp.x = j;
		temp.length = 0;
		temp.blob = j;
		temp.value = 0;
		dest->runs[dest_indice] = temp;
		printf("destRuns(last) #%d xs : %d, length : %d, y : %d, value : %d\n", dest_indice, dest->runs[dest_indice].x , dest->runs[dest_indice].length, dest->runs[dest_indice].blob, dest->runs[dest_indice].value);
	}
	printf("src width : %d, height : %d\n", source->width, source->height);
	dest->nbRuns = dest_indice;
	return 0;
}

// Generic algorithm
int camRLEErode(CamRLEImage *source, CamRLEImage *dest, CamRLEImage *elt)
{
    return 0;
}
