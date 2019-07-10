#pragma once
#ifndef GFX_IMAGECOMPRESS_SRC_AMD_BCX_HELPER_HPP_
#define GFX_IMAGECOMPRESS_SRC_AMD_BCX_HELPER_HPP_

Image_CompressAMDBackendOptions const* Image_CompressDefaultAmdOptions();

void CompressRGBBlockBC1(float rgbBlock[4 * 4 * 4],
											uint32_t compressedBlock[2],
											float pfChannelWeights[3],
											bool threeDRefinement = false,
											uint8_t refinementSteps = 1,
											bool bc1UseAlpha = false,
											float bc1AlphaThreshold = 0.5f);
void CompressRGBBlock(float rgbBlock[4 * 4 * 4],
												 uint32_t compressedBlock[2],
												 float pfChannelWeights[3],
												 bool threeDRefinement = false,
												 uint8_t refinementSteps = 1);
void CompressRGBABlock(float rgbaBlock[4 * 4 * 4],
											 uint32_t compressedBlock[4],
											 float *pfChannelWeights,
											 bool threeDRefinement = false,
											 uint8_t refinementSteps = 1);

void CompressRGBABlock_ExplicitAlpha(float rgbaBlock[4 * 4 * 4],
																		 uint32_t compressedBlock[4],
																		 float *pfChannelWeights,
																		 bool threeDRefinement = false,
																		 uint8_t refinementSteps = 1);
void CompressAlphaBlock(float alphaBlock[4*4*4], uint32_t compressedBlock[2]);
void CompressExplicitAlphaBlock(float const alphaBlock[4 * 4],
																uint32_t compressedBlock[2]);

#endif // end GFX_IMAGECOMPRESS_SRC_AMD_BCX_COMPRESSOR_HPP_
