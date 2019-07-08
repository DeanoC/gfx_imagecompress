#pragma once
#ifndef GFX_IMAGECOMPRESS_SRC_BLOCK_UTILS_HPP_
#define GFX_IMAGECOMPRESS_SRC_BLOCK_UTILS_HPP_

#include "gfx_image/image.hpp"

namespace ImageCompress {
void ReadNxNBlock(Image_ImageHeader const *src,
													uint32_t blockWidth,
													uint32_t blockHeight,
													float *dst,
													uint32_t x,
													uint32_t y,
													uint32_t w);

void WriteNxNBlock(Image_ImageHeader const *dst,
									uint32_t blockWidth,
									uint32_t blockHeight,
									void const *blockData,
									uint32_t blockByteCount,
									uint32_t x,
									uint32_t y,
									uint32_t w);

void CalculateColourWeightings(float block[4 * 4 * 4], float weights[3]);
}

#endif // end GFX_IMAGECOMPRESS_SRC_BLOCK_UTILS_HPP_
