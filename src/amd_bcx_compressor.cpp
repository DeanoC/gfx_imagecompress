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
//
//
//  File Name:   ATIXCodec.cpp
//  Description: Performs the DXT-style block compression
//
//
//////////////////////////////////////////////////////////////////////////////

#include "al2o3_platform/platform.h"
#include "al2o3_cmath/scalar.h"
#include <float.h>
#include "amd_bcx_compressor.hpp"

#define ALIGN_16(x) AL2O3_DEFINE_ALIGNED(x, 16)

#define BLOCK_SIZE MAX_BLOCK

#define EPS        (2.f / 255.f) * (2.f / 255.f)
#define EPS2    3.f * (2.f / 255.f) * (2.f / 255.f)
#ifndef MAX_ERROR
#define MAX_ERROR 128000.f
#endif

#ifndef GBL_SCH_STEP
#define GBL_SCH_STEP_MXS 0.018f
#define GBL_SCH_EXT_MXS 0.1f
#define LCL_SCH_STEP_MXS 0.6f
#define GBL_SCH_STEP_MXQ 0.0175f
#define GBL_SCH_EXT_MXQ 0.154f
#define LCL_SCH_STEP_MXQ 0.45f

#define GBL_SCH_STEP GBL_SCH_STEP_MXS
#define GBL_SCH_EXT  GBL_SCH_EXT_MXS
#define LCL_SCH_STEP LCL_SCH_STEP_MXS
#endif

// Channel indexes
#define AC 3
#define RC 2
#define GC 1
#define BC 0

// Grid precision
#define PIX_GRID 8

#define ConstructColour(r, g, b)  (((r) << 11) | ((g) << 5) | (b))

uint8_t nByteBitsMask[] =
		{
				0x00,
				0x80,
				0xc0,
				0xe0,
				0xf0,
				0xf8,
				0xfc,
				0xfe,
				0xff,
		};

uint32_t ConstructColor(uint8_t R, uint8_t nRedBits, uint8_t G, uint8_t nGreenBits, uint8_t B, uint8_t nBlueBits)
{
	return (    ((R & nByteBitsMask[nRedBits])    << (nGreenBits + nBlueBits - (PIX_GRID - nRedBits))) |
			((G & nByteBitsMask[nGreenBits])<< (nBlueBits - (PIX_GRID - nGreenBits))) |
			((B & nByteBitsMask[nBlueBits]) >> ((PIX_GRID - nBlueBits))));
}

/*--------------------------------------------------------------------------
 3 DIM VECTOR CASE
---------------------------------------------------------------------------*/


/*------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------*/
static int QSortIntCmp(const void * Elem1, const void * Elem2)
{
	return (*(int*)Elem1 - *(int*)Elem2);
}

static int QSortFloatCmp(const void * Elem1, const void * Elem2)
{
	uint32_t* pdwElem1 = (uint32_t*) Elem1;
	uint32_t* pdwElem2 = (uint32_t*) Elem2;
	if(    (pdwElem1[2] > pdwElem2[2]) ||
			(pdwElem1[2] == pdwElem2[2] && pdwElem1[1] > pdwElem2[1]) ||
			(pdwElem1[2] == pdwElem2[2] && pdwElem1[1] == pdwElem2[1] && pdwElem1[0] > pdwElem2[0]))
		return 1;
	else if((pdwElem1[2] < pdwElem2[2]) ||
			(pdwElem1[2] == pdwElem2[2] && pdwElem1[1] < pdwElem2[1]) ||
			(pdwElem1[2] == pdwElem2[2] && pdwElem1[1] == pdwElem2[1] && pdwElem1[0] < pdwElem2[0]))
		return -1;
	else
		return 0;
}

/*------------------------------------------------------------------------------------------------
// this is how the end points is going to be rounded in compressed format
------------------------------------------------------------------------------------------------*/
static void MkRmpOnGrid(float _RmpF[NUM_CHANNELS][NUM_ENDPOINTS], float _MnMx[NUM_CHANNELS][NUM_ENDPOINTS],
												float _Min, float _Max, uint8_t nRedBits, uint8_t nGreenBits, uint8_t nBlueBits)
{
	float Fctrs0[3];
	float Fctrs1[3];

	Fctrs1[RC] = (float)(1 << nRedBits);
	Fctrs1[GC] = (float)(1 << nGreenBits);
	Fctrs1[BC] = (float)(1 << nBlueBits);
	Fctrs0[RC] = (float)(1 << (PIX_GRID-nRedBits));
	Fctrs0[GC] = (float)(1 << (PIX_GRID-nGreenBits));
	Fctrs0[BC] = (float)(1 << (PIX_GRID-nBlueBits));

	for(int j = 0; j < 3; j++)
	{
		for(int k = 0; k < 2; k++)
		{
			_RmpF[j][k] = floor(_MnMx[j][k]);
			if(_RmpF[j][k] <= _Min)
				_RmpF[j][k] = _Min;
			else
			{
				_RmpF[j][k] += floor(128.f / Fctrs1[j]) - floor(_RmpF[j][k] / Fctrs1[j]);
				_RmpF[j][k] = Math_MinF(_RmpF[j][k], _Max);
			}

			_RmpF[j][k] = floor(_RmpF[j][k] / Fctrs0[j]) * Fctrs0[j];
		}
	}
}


/*------------------------------------------------------------------------------------------------
// this is how the end points is going to be look like when decompressed
------------------------------------------------------------------------------------------------*/
inline void MkWkRmpPts(bool *_bEq, float _OutRmpPts[NUM_CHANNELS][NUM_ENDPOINTS],
											 float _InpRmpPts[NUM_CHANNELS][NUM_ENDPOINTS], uint8_t nRedBits, uint8_t nGreenBits, uint8_t nBlueBits)
{
	float Fctrs[3];
	Fctrs[RC] = (float)(1 << nRedBits);
	Fctrs[GC] = (float)(1 << nGreenBits);
	Fctrs[BC] = (float)(1 << nBlueBits);

	*_bEq = true;
	// find whether input ramp is flat
	for(int j = 0; j < 3; j++)
		*_bEq  &= (_InpRmpPts[j][0] == _InpRmpPts[j][1]);

	// end points on the integer grid
	for(int j = 0; j <3; j++)
	{
		for(int k = 0; k <2; k++)
		{
			// Apply the lower bit replication to give full dynamic range
			_OutRmpPts[j][k] = _InpRmpPts[j][k] + floor(_InpRmpPts[j][k] / Fctrs[j]);
			_OutRmpPts[j][k] = Math_MaxF(_OutRmpPts[j][k], 0.f);
			_OutRmpPts[j][k] = Math_MinF(_OutRmpPts[j][k], 255.f);
		}
	}
}

uint32_t dwRndAmount[] = {0, 0, 0, 0, 1, 1, 2, 2, 3};

/*------------------------------------------------------------------------------------------------
1 DIM ramp 
------------------------------------------------------------------------------------------------*/
inline void BldClrRmp(float _Rmp[MAX_POINTS], float _InpRmp[NUM_ENDPOINTS], uint8_t dwNumPoints)
{
	// linear interpolate end points to get the ramp 
	_Rmp[0] = _InpRmp[0];
	_Rmp[dwNumPoints - 1] = _InpRmp[1];
	if(dwNumPoints % 2)
		_Rmp[dwNumPoints] = 1000000.f; // for 3 point ramp; not to select the 4th point as min
	for(int e = 1; e < dwNumPoints - 1; e++)
		_Rmp[e] = floor((_Rmp[0] * (dwNumPoints - 1 - e) + _Rmp[dwNumPoints - 1] * e + dwRndAmount[dwNumPoints])/ (float)(dwNumPoints - 1));
}

/*------------------------------------------------------------------------------------------------
// build 3D ramp
------------------------------------------------------------------------------------------------*/
inline void BldRmp(float _Rmp[NUM_CHANNELS][MAX_POINTS], float _InpRmp[NUM_CHANNELS][NUM_ENDPOINTS],
									 uint8_t dwNumPoints)
{
	for(int j = 0; j < 3; j++)
		BldClrRmp(_Rmp[j], _InpRmp[j], dwNumPoints);
}



/*------------------------------------------------------------------------------------------------
Compute cumulative error for the current cluster
------------------------------------------------------------------------------------------------*/
static float ClstrErr(float _Blk[MAX_BLOCK][NUM_CHANNELS], float _Rpt[MAX_BLOCK],
													 float _Rmp[NUM_CHANNELS][MAX_POINTS], int _NmbClrs, int _blcktp,
													 bool _ConstRamp, float* _pfWeights)
{
	float fError = 0.f;
	int rmp_l = (_ConstRamp) ? 1 : _blcktp;

	// For each colour in the original block, find the closest cluster
	// and compute the comulative error
	for(int i=0; i<_NmbClrs; i++)
	{
		float fShortest = 99999999999.f;

		if(_pfWeights)
			for(int r=0; r < rmp_l; r++)
			{
				// calculate the distance for each component
				float fDistance =    (_Blk[i][RC] - _Rmp[RC][r]) * (_Blk[i][RC] - _Rmp[RC][r]) * _pfWeights[0] +
						(_Blk[i][GC] - _Rmp[GC][r]) * (_Blk[i][GC] - _Rmp[GC][r]) * _pfWeights[1] +
						(_Blk[i][BC] - _Rmp[BC][r]) * (_Blk[i][BC] - _Rmp[BC][r]) * _pfWeights[2];

				if(fDistance < fShortest)
					fShortest = fDistance;
			}
		else
			for(int r=0; r < rmp_l; r++)
			{
				// calculate the distance for each component
				float fDistance =    (_Blk[i][RC] - _Rmp[RC][r]) * (_Blk[i][RC] - _Rmp[RC][r]) +
						(_Blk[i][GC] - _Rmp[GC][r]) * (_Blk[i][GC] - _Rmp[GC][r]) +
						(_Blk[i][BC] - _Rmp[BC][r]) * (_Blk[i][BC] - _Rmp[BC][r]);

				if(fDistance < fShortest)
					fShortest = fDistance;
			}

		// accumulate the error
		fError += fShortest * _Rpt[i];
	}

	return fError;
}

// Compute error and find DXTC indexes for the current cluster
static float ClstrIntnl(float _Blk[MAX_BLOCK][NUM_CHANNELS], uint8_t* _Indxs,
														 float _Rmp[NUM_CHANNELS][MAX_POINTS], int dwBlockSize, uint8_t dwNumPoints,
														 bool _ConstRamp, float* _pfWeights, bool _bUseAlpha)
{
	float Err = 0.f;
	uint8_t rmp_l = (_ConstRamp) ? 1 : dwNumPoints;

	// For each colour in the original block assign it
	// to the closest cluster and compute the cumulative error
	for(int i=0; i< dwBlockSize; i++)
	{
		if(_bUseAlpha && *((uint32_t*) &_Blk[i][AC]) == 0)
			_Indxs[i] = dwNumPoints;
		else
		{
			float shortest = 99999999999.f;
			uint8_t shortestIndex = 0;
			if(_pfWeights)
				for(uint8_t r=0; r < rmp_l; r++)
				{
					// calculate the distance for each component
					float distance =    (_Blk[i][RC] - _Rmp[RC][r]) * (_Blk[i][RC] - _Rmp[RC][r]) * _pfWeights[0] +
							(_Blk[i][GC] - _Rmp[GC][r]) * (_Blk[i][GC] - _Rmp[GC][r]) * _pfWeights[1] +
							(_Blk[i][BC] - _Rmp[BC][r]) * (_Blk[i][BC] - _Rmp[BC][r]) * _pfWeights[2];

					if(distance < shortest)
					{
						shortest = distance;
						shortestIndex = r;
					}
				}
			else
				for(uint8_t r=0; r < rmp_l; r++)
				{
					// calculate the distance for each component
					float distance =    (_Blk[i][RC] - _Rmp[RC][r]) * (_Blk[i][RC] - _Rmp[RC][r]) +
							(_Blk[i][GC] - _Rmp[GC][r]) * (_Blk[i][GC] - _Rmp[GC][r]) +
							(_Blk[i][BC] - _Rmp[BC][r]) * (_Blk[i][BC] - _Rmp[BC][r]);

					if(distance < shortest)
					{
						shortest = distance;
						shortestIndex = r;
					}
				}

			Err += shortest;

			// We have the index of the best cluster, so assign this in the block
			// Reorder indices to match correct DXTC ordering
			if(shortestIndex == dwNumPoints - 1)
				shortestIndex = 1;
			else if(shortestIndex)
				shortestIndex++;
			_Indxs[i] = shortestIndex;
		}
	}

	return Err;
}

/*------------------------------------------------------------------------------------------------
// input ramp is on the coarse grid
------------------------------------------------------------------------------------------------*/
static float ClstrBas(uint8_t* _Indxs, float _Blk[MAX_BLOCK][NUM_CHANNELS],
													 float _InpRmp[NUM_CHANNELS][NUM_ENDPOINTS], int dwBlockSize, uint8_t dwNumPoints, float* _pfWeights,
													 bool _bUseAlpha, uint8_t nRedBits, uint8_t nGreenBits, uint8_t nBlueBits)
{
	// make ramp endpoints the way they'll going to be decompressed
	bool Eq = true;
	float InpRmp[NUM_CHANNELS][NUM_ENDPOINTS];
	MkWkRmpPts(&Eq, InpRmp, _InpRmp, nRedBits, nGreenBits, nBlueBits);

	// build ramp as it would be built by decompressor
	float Rmp[NUM_CHANNELS][MAX_POINTS];
	BldRmp(Rmp, InpRmp, dwNumPoints);

	// clusterize and find a cumulative error
	return ClstrIntnl(_Blk, _Indxs, Rmp, dwBlockSize, dwNumPoints, Eq,  _pfWeights, _bUseAlpha);
}

/*------------------------------------------------------------------------------------------------
Clusterization the way it looks from the DXTC decompressor
------------------------------------------------------------------------------------------------*/

float Clstr(uint32_t block_32[MAX_BLOCK], uint16_t dwBlockSize,
								 uint8_t nEndpoints[3][NUM_ENDPOINTS],
								 uint8_t* pcIndices, uint8_t dwNumPoints,
								 float* _pfWeights, bool _bUseAlpha, uint8_t _nAlphaThreshold,
								 uint8_t nRedBits, uint8_t nGreenBits, uint8_t nBlueBits
)
{
	unsigned int c0 = ConstructColor(nEndpoints[RC][0], nRedBits, nEndpoints[GC][0], nGreenBits, nEndpoints[BC][0], nBlueBits);
	unsigned int c1 = ConstructColor(nEndpoints[RC][1], nRedBits, nEndpoints[GC][1], nGreenBits, nEndpoints[BC][1], nBlueBits);
	unsigned int nEndpointIndex0 = 0;
	unsigned int nEndpointIndex1 = 1;
	if((!(dwNumPoints & 0x1) && c0 <= c1) || ((dwNumPoints & 0x1) && c0 > c1))
	{
		nEndpointIndex0 = 1;
		nEndpointIndex1 = 0;
	}

	float InpRmp[NUM_CHANNELS][NUM_ENDPOINTS];
	InpRmp[RC][0] = (float)nEndpoints[RC][nEndpointIndex0];
	InpRmp[RC][1] = (float)nEndpoints[RC][nEndpointIndex1];
	InpRmp[GC][0] = (float)nEndpoints[GC][nEndpointIndex0];
	InpRmp[GC][1] = (float)nEndpoints[GC][nEndpointIndex1];
	InpRmp[BC][0] = (float)nEndpoints[BC][nEndpointIndex0];
	InpRmp[BC][1] = (float)nEndpoints[BC][nEndpointIndex1];

	uint32_t dwAlphaThreshold = _nAlphaThreshold << 24;
	float Blk[MAX_BLOCK][NUM_CHANNELS];
	for(int i = 0; i < dwBlockSize; i++)
	{
		Blk[i][RC] = (float)((block_32[i] & 0xff0000) >> 16);
		Blk[i][GC] = (float)((block_32[i] & 0xff00) >> 8);
		Blk[i][BC] = (float)(block_32[i] & 0xff);
		if(_bUseAlpha)
			Blk[i][AC] = ((block_32[i] & 0xff000000) >= dwAlphaThreshold) ? 1.f : 0.f;
	}

	return ClstrBas(pcIndices, Blk, InpRmp, dwBlockSize, dwNumPoints, _pfWeights, _bUseAlpha, nRedBits, nGreenBits, nBlueBits);
}

float Clstr(float block_32[MAX_BLOCK*4], uint16_t dwBlockSize,
								 uint8_t nEndpoints[3][NUM_ENDPOINTS],
								 uint8_t* pcIndices, uint8_t dwNumPoints,
								 float* _pfWeights, bool _bUseAlpha, float _fAlphaThreshold,
								 uint8_t nRedBits, uint8_t nGreenBits, uint8_t nBlueBits)
{
	unsigned int c0 = ConstructColor(nEndpoints[RC][0], nRedBits, nEndpoints[GC][0], nGreenBits, nEndpoints[BC][0], nBlueBits);
	unsigned int c1 = ConstructColor(nEndpoints[RC][1], nRedBits, nEndpoints[GC][1], nGreenBits, nEndpoints[BC][1], nBlueBits);
	unsigned int nEndpointIndex0 = 0;
	unsigned int nEndpointIndex1 = 1;
	if((!(dwNumPoints & 0x1) && c0 <= c1) || ((dwNumPoints & 0x1) && c0 > c1))
	{
		nEndpointIndex0 = 1;
		nEndpointIndex1 = 0;
	}

	float InpRmp[NUM_CHANNELS][NUM_ENDPOINTS];
	InpRmp[RC][0] = (float)nEndpoints[RC][nEndpointIndex0];
	InpRmp[RC][1] = (float)nEndpoints[RC][nEndpointIndex1];
	InpRmp[GC][0] = (float)nEndpoints[GC][nEndpointIndex0];
	InpRmp[GC][1] = (float)nEndpoints[GC][nEndpointIndex1];
	InpRmp[BC][0] = (float)nEndpoints[BC][nEndpointIndex0];
	InpRmp[BC][1] = (float)nEndpoints[BC][nEndpointIndex1];

	float fAlphaThreshold = _fAlphaThreshold * 255.f;
	float Blk[MAX_BLOCK][NUM_CHANNELS];
	for(int i = 0; i < dwBlockSize; i++)
	{
		Blk[i][RC] = block_32[(i * 4) + 2];
		Blk[i][GC] = block_32[(i * 4) + 1];
		Blk[i][BC] = block_32[(i * 4)];
		if(_bUseAlpha)
			Blk[i][AC] = (block_32[(i * 4) + 3] >= fAlphaThreshold) ? 1.f : 0.f;
	}

	return ClstrBas(pcIndices, Blk, InpRmp, dwBlockSize, dwNumPoints, _pfWeights, _bUseAlpha, nRedBits, nGreenBits, nBlueBits);
}

static float Refine(float _OutRmpPnts[NUM_CHANNELS][NUM_ENDPOINTS],
												 float _InpRmpPnts[NUM_CHANNELS][NUM_ENDPOINTS],
												 float _Blk[MAX_BLOCK][NUM_CHANNELS], float _Rpt[MAX_BLOCK],
												 int _NmrClrs, uint8_t dwNumPoints, float* _pfWeights,
												 uint8_t nRedBits, uint8_t nGreenBits, uint8_t nBlueBits, uint8_t nRefineSteps);

static float Refine3D(float _OutRmpPnts[NUM_CHANNELS][NUM_ENDPOINTS],
													 float _InpRmpPnts[NUM_CHANNELS][NUM_ENDPOINTS],
													 float _Blk[MAX_BLOCK][NUM_CHANNELS], float _Rpt[MAX_BLOCK],
													 int _NmrClrs, uint8_t dwNumPoints, float* _pfWeights,
													 uint8_t nRedBits, uint8_t nGreenBits, uint8_t nBlueBits, uint8_t nRefineSteps);
/*------------------------------------------------------------------------------------------------
1 dim error
------------------------------------------------------------------------------------------------*/
static float RampSrchW(float _Blck[MAX_BLOCK],
														float _BlckErr[MAX_BLOCK],
														float _Rpt[MAX_BLOCK],
														float _maxerror, float _min_ex, float _max_ex,
														int _NmbClrs,
														int _block)
{
	float error = 0;
	float step = (_max_ex - _min_ex) / (_block - 1);
	float step_h = step * (float)0.5;
	float rstep = (float)1.0f / step;

	for(int i=0; i < _NmbClrs; i++)
	{
		float v;
		// Work out which value in the block this select
		float del;

		if((del = _Blck[i] - _min_ex) <= 0)
			v = _min_ex;
		else if(_Blck[i] -  _max_ex >= 0)
			v = _max_ex;
		else
			v = floor((del + step_h) * rstep) * step + _min_ex;

		// And accumulate the error
		float d = (_Blck[i] - v);
		d *= d;
		float err = _Rpt[i] * d + _BlckErr[i];
		error += err;
		if(_maxerror < error)
		{
			error  = _maxerror;
			break;
		}
	}
	return error;
}

// Find the first approximation of the line
// Assume there is a linear relation 
//   Z = a * X_In
//   Z = b * Y_In
// Find a,b to minimize MSE between Z and Z_In
static void FindAxis(float _outBlk[MAX_BLOCK][NUM_CHANNELS], float fLineDirection[NUM_CHANNELS],
										 float fBlockCenter[NUM_CHANNELS], bool * _pbSmall, float _inpBlk[MAX_BLOCK][NUM_CHANNELS],
										 float _inpRpt[MAX_BLOCK], int nDimensions, int nNumColors)
{
	float Crrl[NUM_CHANNELS];
	float RGB2[NUM_CHANNELS];

	fLineDirection[0] = fLineDirection[1] = fLineDirection[2] = RGB2[0] = RGB2[1] = RGB2[2] =
	Crrl[0] = Crrl[1] = Crrl[2] = fBlockCenter[0] = fBlockCenter[1] = fBlockCenter[2] = 0.f;

	// sum position of all points
	float fNumPoints = 0.f;
	for(int i=0; i < nNumColors; i++)
	{
		fBlockCenter[0] += _inpBlk[i][0] * _inpRpt[i];
		fBlockCenter[1] += _inpBlk[i][1] * _inpRpt[i];
		fBlockCenter[2] += _inpBlk[i][2] * _inpRpt[i];
		fNumPoints += _inpRpt[i];
	}

	// and then average to calculate center coordinate of block
	fBlockCenter[0] /= fNumPoints;
	fBlockCenter[1] /= fNumPoints;
	fBlockCenter[2] /= fNumPoints;

	for(int i = 0; i < nNumColors; i++)
	{
		// calculate output block as offsets around block center
		_outBlk[i][0] = _inpBlk[i][0] - fBlockCenter[0];
		_outBlk[i][1] = _inpBlk[i][1] - fBlockCenter[1];
		_outBlk[i][2] = _inpBlk[i][2] - fBlockCenter[2];

		// compute correlation matrix
		// RGB2 = sum of ((distance from point from center) squared)
		// Crrl = ???????. Seems to be be some calculation based on distance from point center in two dimensions
		for(int j = 0; j < nDimensions; j++)
		{
			RGB2[j] += _outBlk[i][j] * _outBlk[i][j] * _inpRpt[i];
			Crrl[j] += _outBlk[i][j] * _outBlk[i][(j+1)%3] * _inpRpt[i];
		}
	}

	// if set's diameter is small 
	int i0 = 0, i1 = 1;
	float mxRGB2 = 0.f;
	int k = 0, j = 0;
	float fEPS = fNumPoints * EPS;
	for(k = 0, j = 0; j < 3; j++)
	{
		if(RGB2[j] >= fEPS)
			k++;
		else
			RGB2[j] = 0.f;

		if(mxRGB2 < RGB2[j])
		{
			mxRGB2 = RGB2[j];
			i0 = j;
		}
	}

	float fEPS2 = fNumPoints * EPS2;
	*_pbSmall = true;
	for(j = 0; j < 3; j++)
		*_pbSmall &= (RGB2[j] < fEPS2);

	if(*_pbSmall) // all are very small to avoid division on the small determinant
		return;

	if(k == 1) // really only 1 dimension
		fLineDirection[i0]= 1.;
	else if(k == 2) // really only 2 dimensions
	{
		i1 = (RGB2[(i0+1)%3] > 0.f) ? (i0+1)%3 : (i0+2)%3;
		float Crl = (i1 == (i0+1)%3) ? Crrl[i0] : Crrl[(i0+2)%3];
		fLineDirection[i1] = Crl/ RGB2[i0];
		fLineDirection[i0]= 1.;
	}
	else
	{
		float maxDet = 100000.f;
		float Cs[3];
		// select max det for precision
		for(j = 0; j < nDimensions; j++)
		{
			float Det = RGB2[j] * RGB2[(j+1)%3] - Crrl[j] * Crrl[j];
			Cs[j] = abs(Crrl[j]/sqrt(RGB2[j] * RGB2[(j+1)%3]));
			if(maxDet < Det)
			{
				maxDet = Det;
				i0 = j;
			}
		}

		// inverse correl matrix
		//  --      --       --      --
		//  |  A   B |       |  C  -B |
		//  |  B   C |  =>   | -B   A |
		//  --      --       --     --
		float mtrx1[2][2];
		float vc1[2];
		float vc[2];
		vc1[0] = Crrl[(i0 + 2) %3];
		vc1[1] = Crrl[(i0 + 1) %3];
		// C
		mtrx1[0][0] = RGB2[(i0+1)%3];
		// A
		mtrx1[1][1] = RGB2[i0];
		// -B
		mtrx1[1][0] = mtrx1[0][1] = -Crrl[i0];
		// find a solution
		vc[0] = mtrx1[0][0] * vc1[0] + mtrx1[0][1] * vc1[1];
		vc[1] = mtrx1[1][0] * vc1[0] + mtrx1[1][1] * vc1[1];
		// normalize
		vc[0] /= maxDet;
		vc[1] /= maxDet;
		// find a line direction vector
		fLineDirection[i0] = 1.;
		fLineDirection[(i0 + 1) %3] = 1.;
		fLineDirection[(i0 + 2) %3] = vc[0] + vc[1];
	}

	// normalize direction vector
	float Len = fLineDirection[0] * fLineDirection[0] + fLineDirection[1] * fLineDirection[1] + fLineDirection[2] * fLineDirection[2];
	Len = sqrt(Len);

	for(j = 0; j < 3; j++)
		fLineDirection[j] = (Len > 0.f) ? fLineDirection[j] / Len : 0.f;
}

/*-------------------------------------------------------------------------------------------

This C version. For more performance, please, take SSE2 route.
--------------------------------------------------------------------------------------------*/
#define SCH_STPS 3 // number of search steps to make at each end of interval

// coefficient defining a number of grid units in one search step 

static const float sMvF[] = { 0.f, -1.f, 1.f, -2.f, 2.f, -3.f, 3.f, -4.f, 4.f, -5.f, 5.f, -6.f, 6.f, -7.f, 7.f, -8.f, 8.f};

float Refine(float _OutRmpPnts[NUM_CHANNELS][NUM_ENDPOINTS],
									float _InpRmpPnts[NUM_CHANNELS][NUM_ENDPOINTS],
									float _Blk[MAX_BLOCK][NUM_CHANNELS], float _Rpt[MAX_BLOCK],
									int _NmrClrs, uint8_t dwNumPoints, float* _pfWeights,
									uint8_t nRedBits, uint8_t nGreenBits, uint8_t nBlueBits, uint8_t nRefineSteps)
{
	ALIGN_16(float) Rmp[NUM_CHANNELS][MAX_POINTS];

	float Blk[MAX_BLOCK][NUM_CHANNELS];
	for(int i = 0; i < _NmrClrs; i++)
		for(int j = 0; j < 3; j++)
			Blk[i][j] = _Blk[i][j];

	float fWeightRed = _pfWeights ? _pfWeights[0] : 1.f;
	float fWeightGreen = _pfWeights ? _pfWeights[1] : 1.f;
	float fWeightBlue = _pfWeights ? _pfWeights[2] : 1.f;

	// here is our grid
	float Fctrs[3];
	Fctrs[RC] = (float)(1 << (PIX_GRID-nRedBits));
	Fctrs[GC] = (float)(1 << (PIX_GRID-nGreenBits));
	Fctrs[BC] = (float)(1 << (PIX_GRID-nBlueBits));

	float InpRmp0[NUM_CHANNELS][NUM_ENDPOINTS];
	float InpRmp[NUM_CHANNELS][NUM_ENDPOINTS];
	for(int k = 0; k < 2; k++)
		for(int j = 0; j < 3; j++)
			InpRmp0[j][k] = InpRmp[j][k] = _OutRmpPnts[j][k] = _InpRmpPnts[j][k];

	// make ramp endpoints the way they'll going to be decompressed
	// plus check whether the ramp is flat
	bool Eq;
	float WkRmpPts[NUM_CHANNELS][NUM_ENDPOINTS];
	MkWkRmpPts(&Eq, WkRmpPts, InpRmp, nRedBits, nGreenBits, nBlueBits);

	// build ramp for all 3 colors
	BldRmp(Rmp, WkRmpPts, dwNumPoints);

	// clusterize for the current ramp
	float bestE = ClstrErr(Blk, _Rpt, Rmp, _NmrClrs, dwNumPoints, Eq, _pfWeights);
	if(bestE == 0.f || !nRefineSteps)    // if exact, we've done
		return bestE;

	// Tweak each component in isolation and get the best values

	// precompute ramp errors for Green and Blue
	float RmpErr[MAX_POINTS][MAX_BLOCK];
	for(int i=0; i < _NmrClrs; i++)
	{
		for(int r = 0; r < dwNumPoints; r++)
		{
			float DistG = (Rmp[GC][r] - Blk[i][GC]);
			float DistB = (Rmp[BC][r] - Blk[i][BC]);
			RmpErr[r][i] = DistG * DistG * fWeightGreen + DistB * DistB * fWeightBlue;
		}
	}

	// First Red
	float bstC0 = InpRmp0[RC][0];
	float bstC1 = InpRmp0[RC][1];
	int nRefineStart = 0 - (Math_MinF(nRefineSteps, (uint8_t)8));
	int nRefineEnd = Math_MinF(nRefineSteps, (uint8_t)8);
	for(int i = nRefineStart; i <= nRefineEnd; i++)
	{
		for(int j = nRefineStart; j <= nRefineEnd; j++)
		{
			// make a move; both sides of interval.        
			InpRmp[RC][0] = Math_MinF(Math_MaxF(InpRmp0[RC][0] + i * Fctrs[RC], 0.f), 255.f);
			InpRmp[RC][1] = Math_MinF(Math_MaxF(InpRmp0[RC][1] + j * Fctrs[RC], 0.f), 255.f);

			// make ramp endpoints the way they'll going to be decompressed
			// plus check whether the ramp is flat
			MkWkRmpPts(&Eq, WkRmpPts, InpRmp, nRedBits, nGreenBits, nBlueBits);

			// build ramp only for red
			BldClrRmp(Rmp[RC], WkRmpPts[RC], dwNumPoints);

			// compute cumulative error
			float mse = 0.f;
			int rmp_l = (Eq) ? 1 : dwNumPoints;
			for(int k = 0; k < _NmrClrs; k++)
			{
				float MinErr = 10000000.f;
				for(int r = 0; r < rmp_l; r++)
				{
					float Dist = (Rmp[RC][r] - Blk[k][RC]);
					float Err = RmpErr[r][k] + Dist * Dist * fWeightRed;
					MinErr = Math_MinF(MinErr, Err);
				}
				mse += MinErr * _Rpt[k];
			}

			// save if we achieve better result
			if(mse < bestE)
			{
				bstC0 = InpRmp[RC][0];
				bstC1 = InpRmp[RC][1];
				bestE = mse;
			}
		}
	}

	// our best REDs
	InpRmp[RC][0] = bstC0;
	InpRmp[RC][1] = bstC1;

	// make ramp endpoints the way they'll going to be decompressed
	// plus check whether the ramp is flat
	MkWkRmpPts(&Eq, WkRmpPts, InpRmp, nRedBits, nGreenBits, nBlueBits);

	// build ramp only for green
	BldRmp(Rmp, WkRmpPts, dwNumPoints);

	// precompute ramp errors for Red and Blue
	for(int i=0; i < _NmrClrs; i++)
	{
		for(int r = 0; r < dwNumPoints; r++)
		{
			float DistR = (Rmp[RC][r] - Blk[i][RC]);
			float DistB = (Rmp[BC][r] - Blk[i][BC]);
			RmpErr[r][i] = DistR * DistR * fWeightRed + DistB * DistB * fWeightBlue;
		}
	}

	// Now green
	bstC0 = InpRmp0[GC][0];
	bstC1 = InpRmp0[GC][1];
	for(int i = nRefineStart; i <= nRefineEnd; i++)
	{
		for(int j = nRefineStart; j <= nRefineEnd; j++)
		{
			InpRmp[GC][0] = Math_MinF(Math_MaxF(InpRmp0[GC][0] + i * Fctrs[GC], 0.f), 255.f);
			InpRmp[GC][1] = Math_MinF(Math_MaxF(InpRmp0[GC][1] + j * Fctrs[GC], 0.f), 255.f);

			MkWkRmpPts(&Eq, WkRmpPts, InpRmp, nRedBits, nGreenBits, nBlueBits);
			BldClrRmp(Rmp[GC], WkRmpPts[GC], dwNumPoints);

			float mse = 0.f;
			int rmp_l = (Eq) ? 1 : dwNumPoints;
			for(int k = 0; k < _NmrClrs; k++)
			{
				float MinErr = 10000000.f;
				for(int r = 0; r < rmp_l; r++)
				{
					float Dist = (Rmp[GC][r] - Blk[k][GC]);
					float Err = RmpErr[r][k] +  Dist * Dist * fWeightGreen;
					MinErr = Math_MinF(MinErr, Err);
				}
				mse += MinErr * _Rpt[k];
			}

			if(mse < bestE)
			{
				bstC0 = InpRmp[GC][0];
				bstC1 = InpRmp[GC][1];
				bestE = mse;
			}
		}
	}

	// our best GREENs
	InpRmp[GC][0] = bstC0;
	InpRmp[GC][1] = bstC1;

	MkWkRmpPts(&Eq, WkRmpPts, InpRmp, nRedBits, nGreenBits, nBlueBits);
	BldRmp(Rmp, WkRmpPts, dwNumPoints);

	// ramp err for Red and Green
	for(int i=0; i < _NmrClrs; i++)
	{
		for(int r = 0; r < dwNumPoints; r++)
		{
			float DistR = (Rmp[RC][r] - Blk[i][RC]);
			float DistG = (Rmp[GC][r] - Blk[i][GC]);
			RmpErr[r][i] = DistR * DistR * fWeightRed + DistG * DistG * fWeightGreen;
		}
	}

	bstC0 = InpRmp0[BC][0];
	bstC1 = InpRmp0[BC][1];
	// Now blue
	for(int i = nRefineStart; i <= nRefineEnd; i++)
	{
		for(int j = nRefineStart; j <= nRefineEnd; j++)
		{
			InpRmp[BC][0] = Math_MinF(Math_MaxF(InpRmp0[BC][0] + i * Fctrs[BC], 0.f), 255.f);
			InpRmp[BC][1] = Math_MinF(Math_MaxF(InpRmp0[BC][1] + j * Fctrs[BC], 0.f), 255.f);

			MkWkRmpPts(&Eq, WkRmpPts, InpRmp, nRedBits, nGreenBits, nBlueBits);
			BldClrRmp(Rmp[BC], WkRmpPts[BC], dwNumPoints);

			float mse = 0.f;
			int rmp_l = (Eq) ? 1 : dwNumPoints;
			for(int k = 0; k < _NmrClrs; k++)
			{
				float MinErr = 10000000.f;
				for(int r = 0; r < rmp_l; r++)
				{
					float Dist = (Rmp[BC][r] - Blk[k][BC]);
					float Err = RmpErr[r][k] +  Dist * Dist * fWeightBlue;
					MinErr = Math_MinF(MinErr, Err);
				}
				mse += MinErr * _Rpt[k];
			}

			if(mse < bestE)
			{
				bstC0 = InpRmp[BC][0];
				bstC1 = InpRmp[BC][1];
				bestE = mse;
			}
		}
	}

	// our best BLUEs
	InpRmp[BC][0] = bstC0;
	InpRmp[BC][1] = bstC1;

	// return our best choice
	for(int j = 0; j < 3; j++)
		for(int k = 0; k < 2; k++)
			_OutRmpPnts[j][k] = InpRmp[j][k];

	return bestE;
}

float Refine3D(float _OutRmpPnts[NUM_CHANNELS][NUM_ENDPOINTS],
										float _InpRmpPnts[NUM_CHANNELS][NUM_ENDPOINTS],
										float _Blk[MAX_BLOCK][NUM_CHANNELS], float _Rpt[MAX_BLOCK],
										int _NmrClrs, uint8_t dwNumPoints, float* _pfWeights,
										uint8_t nRedBits, uint8_t nGreenBits, uint8_t nBlueBits, uint8_t nRefineSteps)
{
	ALIGN_16(float) Rmp[NUM_CHANNELS][MAX_POINTS];

	float Blk[MAX_BLOCK][NUM_CHANNELS];
	for(int i = 0; i < _NmrClrs; i++)
		for(int j = 0; j < 3; j++)
			Blk[i][j] = _Blk[i][j];

	float fWeightRed = _pfWeights ? _pfWeights[0] : 1.f;
	float fWeightGreen = _pfWeights ? _pfWeights[1] : 1.f;
	float fWeightBlue = _pfWeights ? _pfWeights[2] : 1.f;

	// here is our grid
	float Fctrs[3];
	Fctrs[RC] = (float)(1 << (PIX_GRID-nRedBits));
	Fctrs[GC] = (float)(1 << (PIX_GRID-nGreenBits));
	Fctrs[BC] = (float)(1 << (PIX_GRID-nBlueBits));

	float InpRmp0[NUM_CHANNELS][NUM_ENDPOINTS];
	float InpRmp[NUM_CHANNELS][NUM_ENDPOINTS];
	for(int k = 0; k < 2; k++)
		for(int j = 0; j < 3; j++)
			InpRmp0[j][k] = InpRmp[j][k] = _OutRmpPnts[j][k] = _InpRmpPnts[j][k];

	// make ramp endpoints the way they'll going to be decompressed
	// plus check whether the ramp is flat
	bool Eq;
	float WkRmpPts[NUM_CHANNELS][NUM_ENDPOINTS];
	MkWkRmpPts(&Eq, WkRmpPts, InpRmp, nRedBits, nGreenBits, nBlueBits);

	// build ramp for all 3 colors
	BldRmp(Rmp, WkRmpPts, dwNumPoints);

	// clusterize for the current ramp
	float bestE = ClstrErr(Blk, _Rpt, Rmp, _NmrClrs, dwNumPoints, Eq, _pfWeights);
	if(bestE == 0.f || !nRefineSteps)    // if exact, we've done
		return bestE;

	// Jitter endpoints in each direction
	int nRefineStart = 0 - (Math_MinF(nRefineSteps, (uint8_t)8));
	int nRefineEnd = Math_MinF(nRefineSteps, (uint8_t)8);
	for(int nJitterG0 = nRefineStart; nJitterG0 <= nRefineEnd; nJitterG0++)
	{
		InpRmp[GC][0] = Math_MinF(Math_MaxF(InpRmp0[GC][0] + nJitterG0 * Fctrs[GC], 0.f), 255.f);
		for(int nJitterG1 = nRefineStart; nJitterG1 <= nRefineEnd; nJitterG1++)
		{
			InpRmp[GC][1] = Math_MinF(Math_MaxF(InpRmp0[GC][1] + nJitterG1 * Fctrs[GC], 0.f), 255.f);
			MkWkRmpPts(&Eq, WkRmpPts, InpRmp, nRedBits, nGreenBits, nBlueBits);
			BldClrRmp(Rmp[GC], WkRmpPts[GC], dwNumPoints);

			float RmpErrG[MAX_POINTS][MAX_BLOCK];
			for(int i=0; i < _NmrClrs; i++)
			{
				for(int r = 0; r < dwNumPoints; r++)
				{
					float DistG = (Rmp[GC][r] - Blk[i][GC]);
					RmpErrG[r][i] = DistG * DistG * fWeightGreen;
				}
			}

			for(int nJitterB0 = nRefineStart; nJitterB0 <= nRefineEnd; nJitterB0++)
			{
				InpRmp[BC][0] = Math_MinF(Math_MaxF(InpRmp0[BC][0] + nJitterB0 * Fctrs[BC], 0.f), 255.f);
				for(int nJitterB1 = nRefineStart; nJitterB1 <= nRefineEnd; nJitterB1++)
				{
					InpRmp[BC][1] = Math_MinF(Math_MaxF(InpRmp0[BC][1] + nJitterB1 * Fctrs[BC], 0.f), 255.f);
					MkWkRmpPts(&Eq, WkRmpPts, InpRmp, nRedBits, nGreenBits, nBlueBits);
					BldClrRmp(Rmp[BC], WkRmpPts[BC], dwNumPoints);

					float RmpErr[MAX_POINTS][MAX_BLOCK];
					for(int i=0; i < _NmrClrs; i++)
					{
						for(int r = 0; r < dwNumPoints; r++)
						{
							float DistB = (Rmp[BC][r] - Blk[i][BC]);
							RmpErr[r][i] = RmpErrG[r][i] + DistB * DistB * fWeightBlue;
						}
					}

					for(int nJitterR0 = nRefineStart; nJitterR0 <= nRefineEnd; nJitterR0++)
					{
						InpRmp[RC][0] = Math_MinF(Math_MaxF(InpRmp0[RC][0] + nJitterR0 * Fctrs[RC], 0.f), 255.f);
						for(int nJitterR1 = nRefineStart; nJitterR1 <= nRefineEnd; nJitterR1++)
						{
							InpRmp[RC][1] = Math_MinF(Math_MaxF(InpRmp0[RC][1] + nJitterR1 * Fctrs[RC], 0.f), 255.f);
							MkWkRmpPts(&Eq, WkRmpPts, InpRmp, nRedBits, nGreenBits, nBlueBits);
							BldClrRmp(Rmp[RC], WkRmpPts[RC], dwNumPoints);

							// compute cumulative error
							float mse = 0.f;
							int rmp_l = (Eq) ? 1 : dwNumPoints;
							for(int k = 0; k < _NmrClrs; k++)
							{
								float MinErr = 10000000.f;
								for(int r = 0; r < rmp_l; r++)
								{
									float Dist = (Rmp[RC][r] - Blk[k][RC]);
									float Err = RmpErr[r][k] + Dist * Dist * fWeightRed;
									MinErr = Math_MinF(MinErr, Err);
								}
								mse += MinErr * _Rpt[k];
							}

							// save if we achieve better result
							if(mse < bestE)
							{
								bestE = mse;
								for(int k = 0; k < 2; k++)
									for(int j = 0; j < 3; j++)
										_OutRmpPnts[j][k] = InpRmp[j][k];
							}
						}
					}
				}
			}
		}
	}

	return bestE;
}


//    This is a float point-based compression
//    it assumes that the number of unique colors is already known; input is in [0., 255.] range.
//    This is C version.
static void CompressRGBBlockX(float _RsltRmpPnts[NUM_CHANNELS][NUM_ENDPOINTS],
															float _BlkIn[MAX_BLOCK][NUM_CHANNELS],
															float _Rpt[MAX_BLOCK],
															int _UniqClrs,
															uint8_t dwNumPoints, bool b3DRefinement, uint8_t nRefinementSteps,
															float* _pfWeights,
															uint8_t nRedBits, uint8_t nGreenBits, uint8_t nBlueBits)
{
	ALIGN_16(float) Prj0[MAX_BLOCK];
	ALIGN_16(float) Prj[MAX_BLOCK];
	ALIGN_16(float) PrjErr[MAX_BLOCK];
	ALIGN_16(float) LineDir[NUM_CHANNELS];
	ALIGN_16(float) RmpIndxs[MAX_BLOCK];

	float LineDirG[NUM_CHANNELS];
	float PosG[NUM_ENDPOINTS];
	float Blk[MAX_BLOCK][NUM_CHANNELS];
	float BlkSh[MAX_BLOCK][NUM_CHANNELS];
	float LineDir0[NUM_CHANNELS];
	float Mdl[NUM_CHANNELS];

	float rsltC[NUM_CHANNELS][NUM_ENDPOINTS];
	int i, j, k;

	// down to [0., 1.]
	for(i = 0; i < _UniqClrs; i++)
		for(j = 0; j < 3; j++)
			Blk[i][j] = _BlkIn[i][j] / 255.f;

	bool isDONE = false;

	// as usual if not more then 2 different colors, we've done 
	if(_UniqClrs <= 2)
	{
		for(j = 0; j < 3; j++)
		{
			rsltC[j][0] = _BlkIn[0][j];
			rsltC[j][1] = _BlkIn[_UniqClrs - 1][j];
		}
		isDONE = true;
	}

	if ( !isDONE )
	{
		//    This is our first attempt to find an axis we will go along.
		//    The cumulation is done to find a line minimizing the MSE from the input 3D points.
		bool bSmall = true;
		FindAxis(BlkSh, LineDir0, Mdl, &bSmall, Blk, _Rpt, 3, _UniqClrs);

		//    While trying to find the axis we found that the diameter of the input set is quite small.
		//    Do not bother.
		if(bSmall)
		{
			for(j = 0; j < 3; j++)
			{
				rsltC[j][0] = _BlkIn[0][j];
				rsltC[j][1] = _BlkIn[_UniqClrs - 1][j];
			}
			isDONE = true;
		}
	}

	// GCC is being an awful being when it comes to goto-jumps.
	// So please bear with this.
	if ( !isDONE )
	{
		float ErrG = 10000000.f;
		float PrjBnd[NUM_ENDPOINTS];
		ALIGN_16(float) PreMRep[MAX_BLOCK];
		for(j =0; j < 3; j++)
			LineDir[j] = LineDir0[j];

		//    Here is the main loop.
		//    1. Project input set on the axis in consideration.
		//    2. Run 1 dimensional search (see scalar case) to find an (sub) optimal pair of end points.
		//    3. Compute the vector of indexes (or clusters) for the current approximate ramp.
		//    4. Present our color channels as 3 16DIM vectors.
		//    5. Find closest approximation of each of 16DIM color vector with the projection of the 16DIM index vector.
		//    6. Plug the projections as a new directional vector for the axis.
		//    7. Goto 1.

		//    D - is 16 dim "index" vector (or 16 DIM vector of indexes - {0, 1/3, 2/3, 0, ...,}, but shifted and normalized).
		//    Ci - is a 16 dim vector of color i.
		//    for each Ci find a scalar Ai such that
		//    (Ai * D - Ci) (Ai * D - Ci) -> min , i.e distance between vector AiD and C is min.
		//    You can think of D as a unit interval(vector) "clusterizer",
		//    and Ai is a scale you need to apply to the clusterizer to 
		//    approximate the Ci vector instead of the unit vector.

		//    Solution is 

		//    Ai = (D . Ci) / (D . D); . - is a dot product.

		//    in 3 dim space Ai(s) represent a line direction, along which
		//    we again try to find (sub)optimal quantizer.

		//    That's what our for(;;) loop is about.
		for(;;)
		{
			//  1. Project input set on the axis in consideration.
			// From Foley & Van Dam: Closest point of approach of a line (P + v) to a point (R) is
			//                            P + ((R-P).v) / (v.v))v
			// The distance along v is therefore (R-P).v / (v.v)
			// (v.v) is 1 if v is a unit vector.
			//
			PrjBnd[0] = 1000.;
			PrjBnd[1] = -1000.;
			for(i = 0; i < MAX_BLOCK; i++)
				Prj0[i] = Prj[i] = PrjErr[i] = PreMRep[i] = 0.f;

			for(i = 0; i < _UniqClrs; i++)
			{
				Prj0[i] = Prj[i] = BlkSh[i][0] * LineDir[0] + BlkSh[i][1] * LineDir[1] + BlkSh[i][2] * LineDir[2];

				PrjErr[i] = (BlkSh[i][0] - LineDir[0] * Prj[i]) * (BlkSh[i][0] - LineDir[0] * Prj[i])
						+ (BlkSh[i][1] - LineDir[1] * Prj[i]) * (BlkSh[i][1] - LineDir[1] * Prj[i])
						+ (BlkSh[i][2] - LineDir[2] * Prj[i]) * (BlkSh[i][2] - LineDir[2] * Prj[i]);

				PrjBnd[0] = Math_MinF(PrjBnd[0], Prj[i]);
				PrjBnd[1] = Math_MaxF(PrjBnd[1], Prj[i]);
			}

			//  2. Run 1 dimensional search (see scalar case) to find an (sub) optimal pair of end points.

			// min and max of the search interval
			float Scl[NUM_ENDPOINTS];
			Scl[0] = PrjBnd[0] - (PrjBnd[1] - PrjBnd[0]) * 0.125f;;
			Scl[1] = PrjBnd[1] + (PrjBnd[1] - PrjBnd[0]) * 0.125f;;

			// compute scaling factor to scale down the search interval to [0.,1] 
			const float Scl2 = (Scl[1] - Scl[0]) * (Scl[1] - Scl[0]);
			const float overScl = 1.f/(Scl[1] - Scl[0]);

			for(i = 0; i < _UniqClrs; i++)
			{
				// scale them
				Prj[i] = (Prj[i] - Scl[0]) * overScl;
				// premultiply the scale squire to plug into error computation later
				PreMRep[i] = _Rpt[i] * Scl2;
			}

			// scale first approximation of end points
			for(k = 0; k <2; k++)
				PrjBnd[k] = (PrjBnd[k] - Scl[0]) * overScl;

			float Err = MAX_ERROR;

			// search step
			static const float stp = 0.025f;

			// low Start/End; high Start/End
			const float lS = (PrjBnd[0] - 2.f * stp > 0.f) ?  PrjBnd[0] - 2.f * stp : 0.f;
			const float hE = (PrjBnd[1] + 2.f * stp < 1.f) ?  PrjBnd[1] + 2.f * stp : 1.f;

			// find the best endpoints 
			float Pos[NUM_ENDPOINTS];
			float lP, hP;
			int l, h;
			for(l = 0, lP = lS; l < 8; l++, lP += stp)
			{
				for(h = 0, hP = hE; h < 8; h++, hP -= stp)
				{
					float err = Err;
					// compute an error for the current pair of end points.
					err = RampSrchW(Prj, PrjErr, PreMRep, err, lP, hP, _UniqClrs, dwNumPoints);

					if(err < Err)
					{
						// save better result
						Err = err;
						Pos[0] = lP;
						Pos[1] = hP;
					}
				}
			}

			// inverse the scaling
			for(k = 0; k < 2; k++)
				Pos[k] = Pos[k] * (Scl[1] - Scl[0])+ Scl[0];

			// did we find somthing better from the previous run?
			if(Err + 0.001 < ErrG)
			{
				// yes, remember it
				ErrG = Err;
				LineDirG[0] =  LineDir[0];
				LineDirG[1] =  LineDir[1];
				LineDirG[2] =  LineDir[2];
				PosG[0] = Pos[0];
				PosG[1] = Pos[1];
				//  3. Compute the vector of indexes (or clusters) for the current approximate ramp.
				// indexes
				const float step = (Pos[1] - Pos[0]) / (float)(dwNumPoints - 1);
				const float step_h = step * (float)0.5;
				const float rstep = (float)1.0f / step;
				const float overBlkTp = 1.f/  (float)(dwNumPoints - 1) ;

				// here the index vector is computed, 
				// shifted and normalized
				float indxAvrg = (float)(dwNumPoints - 1) / 2.f;

				for(i=0; i < _UniqClrs; i++)
				{
					float del;
					//int n = (int)((b - _min_ex + (step*0.5f)) * rstep);
					if((del = Prj0[i] - Pos[0]) <= 0)
						RmpIndxs[i] = 0.f;
					else if(Prj0[i] -  Pos[1] >= 0)
						RmpIndxs[i] = (float)(dwNumPoints - 1);
					else
						RmpIndxs[i] = floor((del + step_h) * rstep);
					// shift and normalization
					RmpIndxs[i] = (RmpIndxs[i] - indxAvrg) * overBlkTp;
				}

				//  4. Present our color channels as 3 16DIM vectors.
				//  5. Find closest aproximation of each of 16DIM color vector with the pojection of the 16DIM index vector.
				float Crs[3], Len, Len2;
				for(i = 0, Crs[0] = Crs[1] = Crs[2] = Len = 0.f; i < _UniqClrs; i++)
				{
					const float PreMlt = RmpIndxs[i] * _Rpt[i];
					Len += RmpIndxs[i] * PreMlt;
					for(j = 0; j < 3; j++)
						Crs[j] += BlkSh[i][j] * PreMlt;
				}

				LineDir[0] = LineDir[1] = LineDir[2] = 0.f;
				if(Len > 0.f)
				{
					LineDir[0] = Crs[0]/ Len;
					LineDir[1] = Crs[1]/ Len;
					LineDir[2] = Crs[2]/ Len;

					//  6. Plug the projections as a new directional vector for the axis.
					//  7. Goto 1.
					Len2 = LineDir[0] * LineDir[0] + LineDir[1] * LineDir[1] + LineDir[2] * LineDir[2];
					Len2 = sqrt(Len2);

					LineDir[0] /= Len2;
					LineDir[1] /= Len2;
					LineDir[2] /= Len2;
				}
			}
			else // We was not able to find anything better.  Drop dead.
				break;
		}

		// inverse transform to find end-points of 3-color ramp
		for(k = 0; k < 2; k++)
			for(j = 0; j < 3; j++)
				rsltC[j][k] = (PosG[k] * LineDirG[j]  + Mdl[j]) * 255.f;
	}

	// We've dealt with (almost) unrestricted full precision realm.
	// Now back to the dirty digital world.

	// round the end points to make them look like compressed ones
	float inpRmpEndPts[NUM_CHANNELS][NUM_ENDPOINTS];
	MkRmpOnGrid(inpRmpEndPts, rsltC, 0.f, 255.f, nRedBits, nGreenBits, nBlueBits);

	//    This not a small procedure squeezes and stretches the ramp along each axis (R,G,B) separately while other 2 are fixed.
	//    It does it only over coarse grid - 565 that is. It tries to squeeze more precision for the real world ramp.
	if(b3DRefinement)
		Refine3D(_RsltRmpPnts, inpRmpEndPts, _BlkIn, _Rpt, _UniqClrs, dwNumPoints, _pfWeights, nRedBits, nGreenBits, nBlueBits, nRefinementSteps);
	else
		Refine(_RsltRmpPnts, inpRmpEndPts, _BlkIn, _Rpt, _UniqClrs, dwNumPoints, _pfWeights, nRedBits, nGreenBits, nBlueBits, nRefinementSteps);
}

/*--------------------------------------------------------------------------------------------------------

--------------------------------------------------------------------------------------------------------*/

float CompRGBBlock(float* block_32, uint16_t dwBlockSize,
												uint8_t nRedBits, uint8_t nGreenBits, uint8_t nBlueBits,
												uint8_t nEndpoints[3][NUM_ENDPOINTS], uint8_t* pcIndices, uint8_t dwNumPoints,
												bool b3DRefinement, uint8_t nRefinementSteps, float* _pfChannelWeights,
												bool _bUseAlpha, float _fAlphaThreshold)
{
	ALIGN_16(float) Rpt[MAX_BLOCK];
	ALIGN_16(float) BlkIn[MAX_BLOCK][NUM_CHANNELS];

	memset(Rpt, 0, sizeof(Rpt));
	memset(BlkIn, 0, sizeof(BlkIn));

	uint32_t dwColors = 0;
	float fBlk[BLOCK_SIZE*4];
	for(uint32_t i = 0; i < dwBlockSize; i++)
		if(!_bUseAlpha || (block_32[(i* 4) + 3] >= _fAlphaThreshold))
		{
			fBlk[(dwColors* 4) + 0] = block_32[(i*4) + 2];
			fBlk[(dwColors* 4) + 1] = block_32[(i*4) + 1];
			fBlk[(dwColors* 4) + 2] = block_32[(i*4) + 0];
			fBlk[(dwColors* 4) + 3] = 0.f;
			dwColors++;
		}

	// Do we have any colors ?
	if(dwColors)
	{
		bool bHasAlpha = (dwColors != dwBlockSize);
		if(bHasAlpha && _bUseAlpha && !(dwNumPoints & 0x1))
			return FLT_MAX;

		//  Here we are computing an uniq number of colors.
		//  For each uniq value we compute the number of it appearences.
		qsort((void *)fBlk, (size_t)dwColors, 4 * sizeof(float), QSortFloatCmp);

		float new_p[NUM_CHANNELS];
		uint32_t dwUniqueColors = 0;
		memcpy(&BlkIn[0], &fBlk[0], 4 * sizeof(float));
		memcpy(&new_p, &fBlk[0], 4 * sizeof(float));
		Rpt[dwUniqueColors] = 1.f;
		for(uint32_t i = 1; i < dwColors; i++)
		{
			if(memcmp(&new_p, &fBlk[i*4], 4 * sizeof(float)) != 0)
			{
				dwUniqueColors++;
				memcpy(&BlkIn[dwUniqueColors], &fBlk[i*4], 4 * sizeof(float));
				memcpy(&new_p, &fBlk[i*4], 4 * sizeof(float));
				Rpt[dwUniqueColors] = 1.f;
			}
			else
				Rpt[dwUniqueColors] += 1.f;
		}
		dwUniqueColors++;

		// switch to float
		for(uint32_t i=0; i < dwUniqueColors; i++)
			for(uint32_t j=0; j < 4; j++)
				BlkIn[i][j] *= 255.0;

		float rsltC[NUM_CHANNELS][NUM_ENDPOINTS];
		CompressRGBBlockX(rsltC, BlkIn, Rpt, dwUniqueColors, dwNumPoints, b3DRefinement, nRefinementSteps,
											_pfChannelWeights, nRedBits, nGreenBits, nBlueBits);

		// return to integer realm
		for(int i = 0; i < 3; i++)
			for(int j = 0; j < 2; j++)
				nEndpoints[i][j] =  (uint8_t)rsltC[i][j];

		for(uint32_t i = 0; i < dwBlockSize; i++)
		{
			fBlk[(i* 4) + 0] = block_32[(i*4) + 2] * 255.0f;
			fBlk[(i* 4) + 1] = block_32[(i*4) + 1] * 255.0f;
			fBlk[(i* 4) + 2] = block_32[(i*4) + 0] * 255.0f;
			fBlk[(i* 4) + 3] = block_32[(i*4) + 3] * 255.0f;
		}

		return Clstr(fBlk, dwBlockSize, nEndpoints, pcIndices, dwNumPoints, _pfChannelWeights, _bUseAlpha, _fAlphaThreshold,
								 nRedBits, nGreenBits, nBlueBits);
	}
	else
	{
		// All colors transparent
		nEndpoints[0][0] = nEndpoints[1][0] = nEndpoints[2][0] = 0;
		nEndpoints[0][1] = nEndpoints[1][1] = nEndpoints[2][1] = 0xff;
		memset(pcIndices, 0xff, dwBlockSize);
		return 0.0;
	}
}

/*--------------------------------------------------------------------------------------------------------

--------------------------------------------------------------------------------------------------------*/

float CompRGBBlock(uint32_t* block_32, uint16_t dwBlockSize,
												uint8_t nRedBits, uint8_t nGreenBits, uint8_t nBlueBits,
												uint8_t nEndpoints[3][NUM_ENDPOINTS], uint8_t* pcIndices, uint8_t dwNumPoints,
												bool _bUseSSE2, bool b3DRefinement, uint8_t nRefinementSteps, float* _pfChannelWeights,
												bool _bUseAlpha, uint8_t _nAlphaThreshold)
{
	ALIGN_16(float) Rpt[BLOCK_SIZE];
	ALIGN_16(float) BlkIn[BLOCK_SIZE][NUM_CHANNELS];

	memset(Rpt, 0, sizeof(Rpt));
	memset(BlkIn, 0, sizeof(BlkIn));

	uint32_t dwAlphaThreshold = _nAlphaThreshold << 24;
	uint32_t dwColors = 0;
	uint32_t dwBlk[BLOCK_SIZE];
	for(uint32_t i = 0; i < dwBlockSize; i++)
		if(!_bUseAlpha || (block_32[i] & 0xff000000) >= dwAlphaThreshold)
			dwBlk[dwColors++] = block_32[i] | 0xff000000;

	// Do we have any colors ?
	if(dwColors)
	{
		bool bHasAlpha = (dwColors != dwBlockSize);
		if(bHasAlpha && _bUseAlpha && !(dwNumPoints & 0x1))
			return FLT_MAX;

		// Here we are computing an unique number of colors.
		// For each unique value we compute the number of it appearences.
		qsort((void *)dwBlk, (size_t)dwColors, sizeof(uint32_t), QSortIntCmp);

		uint32_t new_p;
		uint32_t dwBlkU[BLOCK_SIZE];
		uint32_t dwUniqueColors = 0;
		new_p = dwBlkU[0] = dwBlk[0];
		Rpt[dwUniqueColors] = 1.f;
		for(uint32_t i = 1; i < dwColors; i++)
		{
			if(new_p != dwBlk[i])
			{
				dwUniqueColors++;
				new_p = dwBlkU[dwUniqueColors] = dwBlk[i];
				Rpt[dwUniqueColors] = 1.f;
			}
			else
				Rpt[dwUniqueColors] += 1.f;
		}
		dwUniqueColors++;

		// switch to float
		for(uint32_t i=0; i<dwUniqueColors; i++)
		{
			BlkIn[i][RC] = (float)((dwBlkU[i] >> 16) & 0xff); // R
			BlkIn[i][GC] = (float)((dwBlkU[i] >> 8)  & 0xff); // G
			BlkIn[i][BC] = (float)((dwBlkU[i] >> 0)  & 0xff); // B
			BlkIn[i][AC] =  255.f; // A
		}

		float rsltC[NUM_CHANNELS][NUM_ENDPOINTS];
		CompressRGBBlockX(rsltC, BlkIn, Rpt, dwUniqueColors, dwNumPoints, b3DRefinement, nRefinementSteps,
											_pfChannelWeights, nRedBits, nGreenBits, nBlueBits);

		// return to integer realm
		for(int i = 0; i < 3; i++)
			for(int j = 0; j < 2; j++)
				nEndpoints[i][j] =  (uint8_t)rsltC[i][j];

		return Clstr(block_32, dwBlockSize, nEndpoints, pcIndices, dwNumPoints, _pfChannelWeights, _bUseAlpha,
								 _nAlphaThreshold, nRedBits, nGreenBits, nBlueBits);
	}
	else
	{
		// All colors transparent
		nEndpoints[0][0] = nEndpoints[1][0] = nEndpoints[2][0] = 0;
		nEndpoints[0][1] = nEndpoints[1][1] = nEndpoints[2][1] = 0xff;
		memset(pcIndices, 0xff, dwBlockSize);
		return 0.0;
	}
}

/*----------------------------------------------------------------------------
  END of 3 COLOR case
-----------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------
SCALAR CASE
---------------------------------------------------------------------------*/

static float CompBlock1(float _RmpPnts[NUM_ENDPOINTS],
														 float _Blk[MAX_BLOCK], int _Nmbr,
														 uint8_t dwNumPoints, bool bFixedRampPoints,
														 int _IntPrc = 8,
														 int _FracPrc = 0,
														 bool _bFixedRamp = true,
														 bool _bUseSSE2 = true
);

static  float   Clstr1(uint8_t* pcIndices,
														float _blockIn[MAX_BLOCK],
														float _ramp[NUM_ENDPOINTS],
														int _NmbrClrs,
														int nNumPoints,
														bool bFixedRampPoints,
														int _intPrec = 8,
														int _fracPrec = 0,
														bool _bFixedRamp = true
);

/*------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------*/
static void BldRmp1(float _Rmp[MAX_POINTS], float _InpRmp[NUM_ENDPOINTS], int nNumPoints)
{
	// for 3 point ramp; not to select the 4th point in min
	for(int e = nNumPoints; e < MAX_POINTS; e++)
		_Rmp[e] = 100000.f;

	_Rmp[0] = _InpRmp[0];
	_Rmp[1] = _InpRmp[1];
	for(int e = 1; e < nNumPoints - 1; e++)
		_Rmp[e + 1] = (_Rmp[0] * (nNumPoints - 1 - e) + _Rmp[1] * e)/(float)(nNumPoints - 1);
}
/*--------------------------------------------------------------------------------------------

---------------------------------------------------------------------------------------------*/
static void GetRmp1(float _rampDat[MAX_POINTS], float _ramp[NUM_ENDPOINTS], int nNumPoints,
										bool bFixedRampPoints, int _intPrec, int _fracPrec, bool _bFixedRamp)
{
	if(_ramp[0] == _ramp[1])
		return;

	if((!bFixedRampPoints  && _ramp[0] <= _ramp[1]) || (bFixedRampPoints && _ramp[0] > _ramp[1]))
	{
		float t = _ramp[0];
		_ramp[0] = _ramp[1];
		_ramp[1] = t;
	}

	_rampDat[0] = _ramp[0];
	_rampDat[1] = _ramp[1];

	float IntFctr = (float) (1 << _intPrec);
	float FracFctr = (float) (1 << _fracPrec);

	float ramp[NUM_ENDPOINTS];
	ramp[0] = _ramp[0] * FracFctr;
	ramp[1] = _ramp[1] * FracFctr;

	BldRmp1(_rampDat, ramp, nNumPoints);
	if(bFixedRampPoints)
	{
		_rampDat[nNumPoints] = 0.;
		_rampDat[nNumPoints+1] = FracFctr * IntFctr - 1.f;
	}

	if(_bFixedRamp)
	{
		for(int i = 0; i < nNumPoints; i++)
		{
			_rampDat[i] = floor(_rampDat[i] + 0.5f);
			_rampDat[i] /= FracFctr;
		}
	}
}

/*--------------------------------------------------------------------------------------------

---------------------------------------------------------------------------------------------*/
static float Clstr1(uint8_t* pcIndices, float _blockIn[MAX_BLOCK], float _ramp[NUM_ENDPOINTS],
												 int _NmbrClrs, int nNumPoints, bool bFixedRampPoints, int _intPrec, int _fracPrec, bool _bFixedRamp)
{
	float Err = 0.f;
	float alpha[MAX_POINTS];

	for(int i = 0; i < _NmbrClrs; i++)
		pcIndices[i] = 0;

	if(_ramp[0] == _ramp[1])
		return Err;

	if(!_bFixedRamp)
	{
		_intPrec = 8;
		_fracPrec = 0;
	}

	GetRmp1(alpha, _ramp, nNumPoints, bFixedRampPoints, _intPrec, _fracPrec, _bFixedRamp);

	if(bFixedRampPoints)
		nNumPoints += 2;

	const float OverIntFctr = 1.f / ((float) (1 << _intPrec) - 1.f);
	for(int i = 0; i < nNumPoints; i++)
		alpha[i] *= OverIntFctr;

	// For each colour in the original block, calculate its weighted
	// distance from each point in the original and assign it
	// to the closest cluster
	for(int i = 0; i < _NmbrClrs; i++)
	{
		float shortest = 10000000.f;

		// Get the original alpha
		float acur = _blockIn[i];

		for(uint8_t j = 0; j < nNumPoints; j++)
		{
			float adist = (acur - alpha[j]);
			adist *= adist;

			if(adist < shortest)
			{
				shortest = adist;
				pcIndices[i] = j;
			}
		}

		Err += shortest;
	}

	return Err;
}

/*--------------------------------------------------------------------------------------------

---------------------------------------------------------------------------------------------*/
static float RmpSrch1(float _Blk[MAX_BLOCK],
													 float _Rpt[MAX_BLOCK],
													 float _maxerror,
													 float _min_ex,
													 float _max_ex,
													 int _NmbrClrs,
													 uint8_t nNumPoints)
{
	float error = 0;
	const float step = (_max_ex - _min_ex) / (float)(nNumPoints - 1);
	const float step_h = step * 0.5f;
	const float rstep = 1.0f / step;

	for(int i=0; i< _NmbrClrs; i++)
	{
		float v;
		// Work out which value in the block this select
		float del;

		if((del = _Blk[i] - _min_ex) <= 0)
			v = _min_ex;
		else if(_Blk[i] -  _max_ex >= 0)
			v = _max_ex;
		else
			v = (floor((del + step_h) * rstep) * step) + _min_ex;

		// And accumulate the error
		float del2 = (_Blk[i] - v);
		error += del2 * del2 * _Rpt[i];

		// if we've already lost to the previous step bail out
		if(_maxerror < error)
		{
			error  = _maxerror;
			break;
		}
	}
	return error;
}


/*--------------------------------------------------------------------------------------------

---------------------------------------------------------------------------------------------*/

static float Refine1(float _Blk[MAX_BLOCK], float _Rpt[MAX_BLOCK],
													float _MaxError, float& _min_ex, float& _max_ex, float _m_step,
													float _min_bnd, float _max_bnd, int _NmbrClrs,
													uint8_t dwNumPoints, bool _bUseSSE2)
{
	// Start out assuming our endpoints are the min and max values we've determined

	// Attempt a (simple) progressive refinement step to reduce noise in the
	// output image by trying to find a better overall match for the endpoints.

	float maxerror = _MaxError;
	float min_ex = _min_ex;
	float max_ex = _max_ex;

	int mode, bestmode;
	do
	{
		float cr_min0 = min_ex;
		float cr_max0 = max_ex;
		for(bestmode = -1, mode = 0; mode < SCH_STPS * SCH_STPS; mode++)
		{
			// check each move (see sStep for direction)
			float cr_min =  min_ex + _m_step * sMvF[mode / SCH_STPS];
			float cr_max =  max_ex + _m_step * sMvF[mode % SCH_STPS];

			cr_min = Math_MaxF(cr_min, _min_bnd);
			cr_max = Math_MinF(cr_max, _max_bnd);

			float error;
			error = RmpSrch1(_Blk, _Rpt, maxerror, cr_min, cr_max, _NmbrClrs, dwNumPoints);

			if(error < maxerror)
			{
				maxerror = error;
				bestmode = mode;
				cr_min0 = cr_min;
				cr_max0 = cr_max;
			}
		}

		if(bestmode != -1)
		{
			// make move (see sStep for direction)
			min_ex = cr_min0;
			max_ex = cr_max0;
		}
	} while(bestmode != -1);

	_min_ex = min_ex;
	_max_ex = max_ex;

	return maxerror;
}

static int QSortFCmp(const void * Elem1, const void * Elem2)
{
	int ret = 0;

	if(*(float*)Elem1 - *(float*)Elem2 < 0.)
		ret = -1;
	else if(*(float*)Elem1 - *(float*)Elem2 > 0.)
		ret = 1;
	return ret;
}


/*--------------------------------------------------------------------------------------------
// input [0,1]
static float CompBlock1(float _RmpPnts[NUM_ENDPOINTS], [OUT] Min amd Max value of the ramp in float
                                                     format in [0., (1 << _IntPrc) - 1] range


---------------------------------------------------------------------------------------------*/
/*
   this is the case when the input data and possible ramp values are all on integer grid.
*/
#define _INT_GRID (_bFixedRamp && _FracPrc == 0)

static float CompBlock1(float _RmpPnts[NUM_ENDPOINTS], float _Blk[MAX_BLOCK], int _Nmbr,
														 uint8_t dwNumPoints, bool bFixedRampPoints,
														 int _IntPrc, int _FracPrc, bool _bFixedRamp, bool _bUseSSE2)
{
	float fMaxError = 0.f;

	float Ramp[NUM_ENDPOINTS];

	float IntFctr = (float)(1 << _IntPrc);
	//    float FracFctr = (float)(1 << _FracPrc);

	ALIGN_16(float) afUniqueValues[MAX_BLOCK];
	ALIGN_16(float) afValueRepeats[MAX_BLOCK];
	for(int i = 0; i < MAX_BLOCK; i++)
		afUniqueValues[i] = afValueRepeats[i] = 0.f;

	// For each unique value we compute the number of it appearances.
	float fBlk[MAX_BLOCK];
	memcpy(fBlk, _Blk, _Nmbr * sizeof(float));

	// sort the input
	qsort((void *)fBlk, (size_t)_Nmbr, sizeof(float), QSortFCmp);

	float new_p = -2.;

	int N0s = 0, N1s = 0;
	uint32_t dwUniqueValues = 0;
	afUniqueValues[0] = 0.f;

	bool requiresCalculation = true;

	if(bFixedRampPoints)
	{
		for(int i = 0; i < _Nmbr; i++)
		{
			if(new_p != fBlk[i])
			{
				new_p = fBlk[i];
				if(new_p <= 1.5/255.)
					N0s++;
				else if(new_p >= 253.5/255.)
					N1s++;
				else
				{
					afUniqueValues[dwUniqueValues] = fBlk[i];
					afValueRepeats[dwUniqueValues] = 1.f;
					dwUniqueValues++;
				}
			}
			else if(dwUniqueValues > 0 && afUniqueValues[dwUniqueValues - 1] == new_p)
				afValueRepeats[dwUniqueValues - 1] += 1.f;
		}

		// if number of unique colors is less or eq 2 we've done either, but we know that we may have 0s and/or 1s as well. 
		// To avoid for the ramp to be considered flat we invented couple entries on the way.
		if(dwUniqueValues <= 2)
		{
			if(dwUniqueValues == 2) // if 2, take them
			{
				Ramp[0]  = floor(afUniqueValues[0] * (IntFctr - 1) + 0.5f);
				Ramp[1]  = floor(afUniqueValues[1] * (IntFctr - 1) + 0.5f);
			}
			else if(dwUniqueValues == 1) // if 1, add another one
			{
				Ramp[0]  = floor(afUniqueValues[0] * (IntFctr - 1) + 0.5f);
				Ramp[1] = Ramp[0] + 1.f;
			}
			else // if 0, invent them 
			{
				Ramp[0]  = 128.f;
				Ramp[1] = Ramp[0] + 1.f;
			}

			fMaxError = 0.f;
			requiresCalculation = false;
		}
	}
	else
	{
		for(int i = 0; i < _Nmbr; i++)
		{
			if(new_p != fBlk[i])
			{
				afUniqueValues[dwUniqueValues] = new_p = fBlk[i];
				afValueRepeats[dwUniqueValues] = 1.f;
				dwUniqueValues++;
			}
			else
				afValueRepeats[dwUniqueValues - 1] += 1.f;
		}

		// if number of unique colors is less or eq 2, we've done 
		if(dwUniqueValues <= 2)
		{
			Ramp[0] = floor(afUniqueValues[0] * (IntFctr - 1) + 0.5f);
			if(dwUniqueValues == 1)
				Ramp[1] = Ramp[0] + 1.f;
			else
				Ramp[1] = floor(afUniqueValues[1] * (IntFctr - 1) + 0.5f);
			fMaxError = 0.f;
			requiresCalculation = false;
		}
	}

	if ( requiresCalculation )
	{
		float min_ex  = afUniqueValues[0];
		float max_ex  = afUniqueValues[dwUniqueValues - 1];
		float min_bnd = 0, max_bnd = 1.;
		float min_r = min_ex, max_r = max_ex;
		float gbl_l = 0, gbl_r = 0;
		float cntr = (min_r + max_r)/2;

		float gbl_err = MAX_ERROR;
		// Trying to avoid unnecessary calculations. Heuristics: after some analisis it appears 
		// that in integer case, if the input interval not more then 48 we won't get much better

		bool wantsSearch = !( _INT_GRID && max_ex - min_ex <= 48.f / IntFctr );

		if ( wantsSearch )
		{
			// Search.
			// 1. take the vicinities of both low and high bound of the input interval.
			// 2. setup some search step
			// 3. find the new low and high bound which provides an (sub) optimal (infinite precision) clusterization.
			float gbl_llb = (min_bnd >  min_r - GBL_SCH_EXT) ? min_bnd : min_r - GBL_SCH_EXT;
			float gbl_rrb = (max_bnd <  max_r + GBL_SCH_EXT) ? max_bnd : max_r + GBL_SCH_EXT;
			float gbl_lrb = (cntr <  min_r + GBL_SCH_EXT) ? cntr : min_r + GBL_SCH_EXT;
			float gbl_rlb = (cntr >  max_r - GBL_SCH_EXT) ? cntr : max_r - GBL_SCH_EXT;
			for(float step_l = gbl_llb; step_l < gbl_lrb ; step_l+= GBL_SCH_STEP)
			{
				for(float step_r = gbl_rrb; gbl_rlb <= step_r; step_r-=GBL_SCH_STEP)
				{
					float sch_err;
					sch_err = RmpSrch1(afUniqueValues, afValueRepeats, gbl_err, step_l, step_r, dwUniqueValues, dwNumPoints);
					if(sch_err < gbl_err)
					{
						gbl_err = sch_err;
						gbl_l = step_l;
						gbl_r = step_r;
					}
				}
			}

			min_r = gbl_l;
			max_r = gbl_r;
		}

		// This is a refinement call. The function tries to make several small stretches or squashes to 
		// minimize quantization error.
		float m_step = LCL_SCH_STEP/ IntFctr;
		fMaxError = Refine1(afUniqueValues, afValueRepeats, gbl_err, min_r, max_r, m_step, min_bnd, max_bnd, dwUniqueValues,
												dwNumPoints, _bUseSSE2);

		min_ex = min_r;
		max_ex = max_r;

		max_ex *= (IntFctr - 1);
		min_ex *= (IntFctr - 1);
		/*
		this one is tricky. for the float or high fractional precision ramp it tries to avoid
		for the ramp to be collapsed into one integer number after rounding.
		Notice the condition. There is a difference between max_ex and min_ex but after rounding 
		they may collapse into the same integer.
		
		So we try to run the same refinement procedure but with starting position on the integer grid
		and step equal 1.
		*/
		if(!_INT_GRID && max_ex - min_ex > 0. && floor(min_ex + 0.5f) == floor(max_ex + 0.5f))
		{
			m_step = 1.;
			gbl_err = MAX_ERROR;
			for(uint32_t i = 0; i < dwUniqueValues; i++)
				afUniqueValues[i] *= (IntFctr - 1);

			max_ex = min_ex = floor(min_ex + 0.5f);

			gbl_err = Refine1(afUniqueValues, afValueRepeats, gbl_err, min_ex, max_ex, m_step, 0.f, 255.f, dwUniqueValues, dwNumPoints, _bUseSSE2);

			fMaxError = gbl_err;

		}
		Ramp[1] = floor(max_ex + 0.5f);
		Ramp[0] = floor(min_ex + 0.5f);
	}

	// Ensure that the two endpoints are not the same
	// This is legal but serves no need & can break some optimizations in the compressor
	if(Ramp[0] == Ramp[1])
	{
		if(Ramp[1] < 255.f)
			Ramp[1]++;
		else
			Ramp[1]--;
	}
	_RmpPnts[0] = Ramp[0];
	_RmpPnts[1] = Ramp[1];

	return fMaxError;
}

/*--------------------------------------------------------------------------------------------
// input [0,1]
void CompBlock1X(float* _Blk, [IN] scalar data block (alphas or normals) in float format 
                 CMP_uint32_t blockCompressed[NUM_ENDPOINTS],  [OUT] compressed data in DXT5 alpha foramt
                 int _NbrClrs,              [IN] actual number of elements in the block
                 int _intPrec,              [IN} integer precision; it applies both to the input data and
                                                 to the ramp points
                 int _fracPrec,             [IN] fractional precision of the ramp points
                 bool _bFixedRamp,          [IN] non-fixed ramp means we have input and generate
                                                 output as float. fixed ramp means that they are fractional numbers.
                 bool _bUseSSE2             [IN] forces to switch to the SSE2 implementation
               )
---------------------------------------------------------------------------------------------*/

float CompBlock1X(float* _Blk, uint16_t dwBlockSize, uint8_t nEndpoints[2], uint8_t* pcIndices,
											 uint8_t dwNumPoints, bool bFixedRampPoints, bool _bUseSSE2, int _intPrec, int _fracPrec, bool _bFixedRamp)
{
	// just to make them initialized
	if(!_bFixedRamp)
	{
		_intPrec = 8;
		_fracPrec = 0;
	}

	// this one makes the bulk of the work
	float Ramp[NUM_ENDPOINTS];
	CompBlock1(Ramp, _Blk, dwBlockSize, dwNumPoints, bFixedRampPoints, _intPrec, _fracPrec, _bFixedRamp, _bUseSSE2);

	// final clusterization applied
	float fError = Clstr1(pcIndices, _Blk, Ramp, dwBlockSize, dwNumPoints, bFixedRampPoints, _intPrec, _fracPrec, _bFixedRamp);
	nEndpoints[0] = (uint8_t)Ramp[0];
	nEndpoints[1] = (uint8_t)Ramp[1];

	return fError;
}

/*--------------------------------------------------------------------------------------------
// input [0,255]
void CompBlock1X(uint8_t* _Blk, [IN] scalar data block (alphas or normals) in 8 bits format 
                 CMP_uint32_t blockCompressed[NUM_ENDPOINTS],  [OUT] compressed data in DXT5 alpha foramt
                 int _NbrClrs,              [IN] actual number of elements in the block
                 int _intPrec,              [IN} integer precision; it applies both to the input data and
                                                 to the ramp points
                 int _fracPrec,             [IN] fractional precision of the ramp points
                 bool _bFixedRamp,          [IN] always true at this point
                 bool _bUseSSE2             [IN] forces to switch to the SSE2-based assembler implementation
               )
---------------------------------------------------------------------------------------------*/

float CompBlock1X(uint8_t* _Blk, uint16_t dwBlockSize, uint8_t nEndpoints[2], uint8_t* pcIndices,
											 uint8_t dwNumPoints, bool bFixedRampPoints, bool _bUseSSE2, int _intPrec, int _fracPrec, bool _bFixedRamp)
{
	// convert the input and call the float equivalent.
	float fBlk[MAX_BLOCK];
	for(int i = 0; i < dwBlockSize; i++)
		fBlk[i] = (float)_Blk[i] / 255.f;

	return CompBlock1X(fBlk, dwBlockSize, nEndpoints, pcIndices, dwNumPoints, bFixedRampPoints, _bUseSSE2, _intPrec, _fracPrec, _bFixedRamp);
}
