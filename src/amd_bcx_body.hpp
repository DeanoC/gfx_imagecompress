#pragma once
#ifndef GFX_IMAGECOMPRESS_SRC_AMD_BCX_BODY_HPP_
#define GFX_IMAGECOMPRESS_SRC_AMD_BCX_BODY_HPP_

#define MAX_BLOCK 64
#define MAX_POINTS 16

#define NUM_CHANNELS 4
#define NUM_ENDPOINTS 2

float CompRGBBlock(float const * block_32, uint16_t dwBlockSize,
									 uint8_t nRedBits, uint8_t nGreenBits, uint8_t nBlueBits,
									 uint8_t nEndpoints[3][NUM_ENDPOINTS], uint8_t* pcIndices, uint8_t dwNumPoints,
									 bool b3DRefinement, uint8_t nRefinementSteps, float const * _pfChannelWeights,
									 bool _bUseAlpha, float _fAlphaThreshold);

float CompBlock1X(float const * _Blk, uint16_t dwBlockSize, uint8_t nEndpoints[2], uint8_t* pcIndices,
									uint8_t dwNumPoints, bool bFixedRampPoints, int _intPrec, int _fracPrec, bool _bFixedRamp);


#endif // end GFX_IMAGECOMPRESS_SRC_AMD_BCX_BODY_HPP_
