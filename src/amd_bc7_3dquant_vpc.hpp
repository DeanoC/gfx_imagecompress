//===============================================================================
// Copyright (c) 2007-2016  Advanced Micro Devices, Inc. All rights reserved.
// Copyright (c) 2004-2006 ATI Technologies Inc.
//===============================================================================
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files(the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions :
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef _3DQUANT_H_INCLUDED
#define _3DQUANT_H_INCLUDED

#include "amd_bc7_partitions.hpp"
typedef struct {
	float d;
	int i;
} a;

void sortProjection(float projection[MAX_ENTRIES], int order[MAX_ENTRIES], int numEntries);
void covariance(float data[][DIMENSION], int numEntries, float cov[DIMENSION][DIMENSION]);
void centerInPlace(float data[][DIMENSION], int numEntries, float mean[DIMENSION]);
void project(float data[][DIMENSION], int numEntries, float vector[DIMENSION], float projection[MAX_ENTRIES]);
void sortProjection(float projection[MAX_ENTRIES], int order[MAX_ENTRIES], int numEntries);
void eigenVector(float cov[DIMENSION][DIMENSION], float vector[DIMENSION]);
float partition2 (float data[][DIMENSION], int numEntries,int index[]);

float optQuantEven(
		float data[MAX_ENTRIES][DIMENSION],
		int numEntries, int numClusters, int index[MAX_ENTRIES],
		float out[MAX_ENTRIES][DIMENSION],
		float direction [DIMENSION],float *step
) ;

float totalError(float data[MAX_ENTRIES][DIMENSION],float data2[MAX_ENTRIES][DIMENSION],int numEntries);
float totalError(float data[MAX_ENTRIES][MAX_DIMENSION_BIG],float data2[MAX_ENTRIES][MAX_DIMENSION_BIG],int numEntries, int dimension);

/****************************************************/
/****************************************************/
//
// Trace dirven continious quntizer
// traceBuilder/ allTrace builder should be called once first to initilize static trace (file static)
//
//    returns resulting error
//
float optQuantTrace(
		float data[MAX_ENTRIES][DIMENSION],    // input data
		int numEntries,                            // number of input points above (not clear about 1, better to avoid)
		int numClusters,                        // number of clusters on the ramp, max 8 (not clear about 1, better to avoid)
		int index[MAX_ENTRIES],                    // output index, if not all points of the ramp used, 0 may not be assigned
		float out[MAX_ENTRIES][DIMENSION],        // resulting quantization
		float direction [DIMENSION],            // direction vector of the ramp (check normalization)
		float *step                            // step size (check normalization)
);

float optQuantTrace_d(
		float data[MAX_ENTRIES][MAX_DIMENSION_BIG],    // input data
		int numEntries,                            // number of input points above (not clear about 1, better to avoid)
		int numClusters,                        // number of clusters on the ramp, max 8 (not clear about 1, better to avoid)
		int index[MAX_ENTRIES],                    // output index, if not all points of the ramp used, 0 may not be assigned
		float out[MAX_ENTRIES][MAX_DIMENSION_BIG],        // resulting quantization
		float direction [MAX_DIMENSION_BIG],            // direction vector of the ramp (check normalization)
		float *step,                            // step size (check normalization)
		int dimension);

/****************************************************/
/****************************************************/
/****************************************************/

//
//    ping-pong style  KATC type (but for RGB x 2)  trace (with fallback to quantAnD (?) driven "continius" quantizer
//
float optQuant2Trace(
		float data[MAX_ENTRIES][DIMENSION],
		int numEntries, int numClusters, int index_[MAX_ENTRIES],
		float out_[MAX_ENTRIES][DIMENSION],
		float direction [DIMENSION],float *step
);

void printStep (void);

/********************************************/
// continious quantizer for 16 clusters, could use bette testing
float optQuantAnD(
		float data[MAX_ENTRIES][DIMENSION],  // 0-255
		int numEntries, int numClusters, int index[MAX_ENTRIES],
		float out[MAX_ENTRIES][DIMENSION],
		float direction [DIMENSION],float *step
);


float optQuantAnD_d(
		float data[MAX_ENTRIES][MAX_DIMENSION_BIG],  // 0-255
		int numEntries, int numClusters, int index[MAX_ENTRIES],
		float out[MAX_ENTRIES][MAX_DIMENSION_BIG],
		float direction [MAX_DIMENSION_BIG],float *step,
		int dimension
);

/********************************************/


float superQuantAnD(
		float data[MAX_ENTRIES][DIMENSION],
		int numEntries, int numClusters, int index[MAX_ENTRIES],
		float out[MAX_ENTRIES][DIMENSION],
		float direction [DIMENSION],float *step
);

int reconstructGetDirConstr (
		float data[MAX_ENTRIES][DIMENSION],
		float mean[DIMENSION],
		int numEntries,
		int numClusters,
		int index[MAX_ENTRIES],
		float direction [DIMENSION],
		float *step,
		float *idxmean,
		int clump
);

// New
void Quant_Init();
void Quant_DeInit();

#endif



