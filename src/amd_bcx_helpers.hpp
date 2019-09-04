#pragma once
#ifndef GFX_IMAGECOMPRESS_SRC_AMD_BCX_HELPER_HPP_
#define GFX_IMAGECOMPRESS_SRC_AMD_BCX_HELPER_HPP_

Image_CompressAMDBackendOptions const* Image_CompressDefaultAmdOptions();

void CompressRGBBlock(float const rgbBlock[4 * 4 * 4],
												 uint32_t compressedBlock[2],
												 float const pfChannelWeights[3],
												 bool threeDRefinement = false,
												 uint8_t refinementSteps = 1);


#endif // end GFX_IMAGECOMPRESS_SRC_AMD_BCX_COMPRESSOR_HPP_
