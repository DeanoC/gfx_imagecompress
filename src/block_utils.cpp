#include "gfx_imagecompress/imagecompress.h"
#include "block_utils.hpp"

namespace ImageCompress {

void ReadNxNBlock(Image_ImageHeader const *src,
									uint32_t blockWidth,
									uint32_t blockHeight,
									float *dst,
									uint32_t sx,
									uint32_t sy,
									uint32_t sw) {

	uint32_t cy = sy;
	for(uint32_t y = 0; y < blockHeight;++y) {
		uint32_t cx = sx;
		for(uint32_t x = 0; x < blockHeight;++x) {
			if(sy + y >= src->height) cy = src->height-1;
			if(sx + x >= src->width) cx = src->width-1;

			size_t const index = Image_CalculateIndex(src, cx, cy, 0, sw);
			Image_PixelD pixel;
			Image_GetPixelAt(src, &pixel, index);
			dst[(((y*4) + x) * 4) + 0] = (float) pixel.r;
			dst[(((y*4) + x) * 4) + 1] = (float) pixel.g;
			dst[(((y*4) + x) * 4) + 2] = (float) pixel.b;
			dst[(((y*4) + x) * 4) + 3] = (float) pixel.a;
			cx++;
		}
		cy++;
	}
}

void WriteNxNBlock(Image_ImageHeader const *dst,
									 uint32_t blockWidth,
									 uint32_t blockHeight,
									 void const *blockData,
									 uint32_t blockByteCount,
									 uint32_t x,
									 uint32_t y,
									 uint32_t w) {
	ASSERT(blockByteCount == (ImageFormat_BitSizeOfBlock(dst->format)/8));
	uint8_t * rawData = (uint8_t *) Image_RawDataPtr(dst);

	size_t const blockIndex = Image_GetBlockIndex(dst, x, y, 0, w);
	uint8_t * const dstPtr = rawData + (blockIndex * blockByteCount);
	memcpy(dstPtr, blockData, blockByteCount);
}

void CalculateColourWeightings(float block[4 * 4 * 4], float weights[3]) {
	static const float baseWeights[3] = {
			0.3086f,
			0.6094f,
			0.0820f
	};
	float medianR = 0.0f, medianG = 0.0f, medianB = 0.0f;

	for (size_t k = 0; k < (4 * 4); k++) {
		block++;
		medianB += *block++;
		medianG += *block++;
		medianR += *block++;
	}

	medianR /= (4 * 4);
	medianG /= (4 * 4);
	medianB /= (4 * 4);

	// Now skew the colour weightings based on the gravity center of the block
	float largest = Math_MaxF(Math_MaxF(medianR, medianG), medianB);

	if (largest > 0) {
		medianR /= largest;
		medianG /= largest;
		medianB /= largest;
	} else
		medianR = medianG = medianB = 1.0f;

	// Scale weightings back up to 1.0f
	float fWeightScale = 1.0f / (baseWeights[0] + baseWeights[1] + baseWeights[2]);
	weights[0] *= baseWeights[0] * fWeightScale;
	weights[1] *= baseWeights[1] * fWeightScale;
	weights[2] *= baseWeights[2] * fWeightScale;
	weights[0] = ((weights[0] * 3 * medianR) + weights[0]) * 0.25f;
	weights[1] = ((weights[1] * 3 * medianG) + weights[1]) * 0.25f;
	weights[2] = ((weights[2] * 3 * medianB) + weights[2]) * 0.25f;
	fWeightScale = 1.0f / (weights[0] + weights[1] + weights[2]);
	weights[0] *= fWeightScale;
	weights[1] *= fWeightScale;
	weights[2] *= fWeightScale;
}

}