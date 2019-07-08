#pragma once
#ifndef GFX_IMAGECOMPRESS_SRC_AMD_BCX_COMPRESSOR_HPP_
#define GFX_IMAGECOMPRESS_SRC_AMD_BCX_COMPRESSOR_HPP_

#define MAX_BLOCK 64
#define MAX_POINTS 16

#define NUM_CHANNELS 4
#define NUM_ENDPOINTS 2

float CompRGBBlock(float* block_32, uint16_t dwBlockSize,
									 uint8_t nRedBits, uint8_t nGreenBits, uint8_t nBlueBits,
									 uint8_t nEndpoints[3][NUM_ENDPOINTS], uint8_t* pcIndices, uint8_t dwNumPoints,
									 bool b3DRefinement, uint8_t nRefinementSteps, float* _pfChannelWeights,
									 bool _bUseAlpha, float _fAlphaThreshold);


#endif // end GFX_IMAGECOMPRESS_SRC_AMD_BCX_COMPRESSOR_HPP_
