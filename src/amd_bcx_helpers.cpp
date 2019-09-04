#include "al2o3_platform/platform.h"
#include "al2o3_cmath/scalar.h"
#include "gfx_imagecompress/imagecompress.h"
#include "amd_bcx_body.hpp"
#include "amd_bcx_helpers.hpp"
#include "block_utils.hpp"
#include <float.h>

//ARGB Channel indexes
#define AC 3
#define RC 2
#define GC 1
#define BC 0

// bits per pixel
#define RG 5
#define GG 6
#define BG 5

// construct a 565 colour
#define ConstructColour(r, g, b)  (((r) << 11) | ((g) << 5) | (b))

Image_CompressAMDBackendOptions const *Image_CompressDefaultAmdOptions() {
	static Image_CompressAMDBackendOptions const defaultAmdOptions{
			false,
			false,
			1,
			0xFF
	};
	return &defaultAmdOptions;
}
static void EncodeAlphaBlock(uint32_t compressedBlock[2], uint8_t nEndpoints[2], uint8_t nIndices[4 * 4]) {
	compressedBlock[0] = ((int) nEndpoints[0]) | (((int) nEndpoints[1]) << 8);
	compressedBlock[1] = 0;

	for (int i = 0; i < 4 * 4; i++) {
		if (i < 5)
			compressedBlock[0] |= (nIndices[i] & 0x7) << (16 + (i * 3));
		else if (i > 5)
			compressedBlock[1] |= (nIndices[i] & 0x7) << (2 + (i - 6) * 3);
		else {
			compressedBlock[0] |= (nIndices[i] & 0x1) << 31;
			compressedBlock[1] |= (nIndices[i] & 0x6) >> 1;
		}
	}
}

void CompressAlphaBlock(float alphaBlock[4 * 4], uint32_t compressedBlock[2]) {
}

AL2O3_EXTERN_C void Image_CompressAMDBC1Block(float const *input,
																							bool adaptiveColourWeights,
																							bool threeDRefinement,
																							uint8_t refinementSteps,
																							float alphaThreshold,
																							void *outv) {

	float weights[3];
	ImageCompress::CalculateColourWeightings(input, weights, adaptiveColourWeights);

	auto out = (uint32_t *) outv;

	uint8_t nEndpoints[2][3][2];
	uint8_t nIndices[2][4 * 4];

	double fError3 = CompRGBBlock(input, 4 * 4,
																RG, GG, BG,
																nEndpoints[0],
																nIndices[0],
																3,
																threeDRefinement,
																refinementSteps,
																weights,
																alphaThreshold > 0.0f,
																alphaThreshold);

	double fError4 = (fError3 == 0.0) ? FLT_MAX :
									 CompRGBBlock(input,
																4 * 4,
																RG, GG, BG,
																nEndpoints[1],
																nIndices[1],
																4,
																threeDRefinement,
																refinementSteps,
																weights,
																alphaThreshold > 0.0f,
																alphaThreshold);

	unsigned int nMethod = (fError3 <= fError4) ? 0 : 1;
	unsigned int c0 = ConstructColour((nEndpoints[nMethod][RC][0] >> (8 - RG)),
																		(nEndpoints[nMethod][GC][0] >> (8 - GG)),
																		(nEndpoints[nMethod][BC][0] >> (8 - BG)));
	unsigned int c1 = ConstructColour((nEndpoints[nMethod][RC][1] >> (8 - RG)),
																		(nEndpoints[nMethod][GC][1] >> (8 - GG)),
																		(nEndpoints[nMethod][BC][1] >> (8 - BG)));
	if ((nMethod == 1 && c0 <= c1) || (nMethod == 0 && c0 > c1))
		out[0] = c1 | (c0 << 16);
	else
		out[0] = c0 | (c1 << 16);

	out[1] = 0;
	for (int i = 0; i < 16; i++)
		out[1] |= (nIndices[nMethod][i] << (2 * i));
}

AL2O3_EXTERN_C void Image_CompressAMDExplictAlphaSingleModeBlock(float const input[ 4 * 4 ], void *out) {
	static const uint8_t EXPLICIT_ALPHA_PIXEL_MASK = 0xf;
	static const uint8_t EXPLICIT_ALPHA_PIXEL_BPP = 4;

	uint32_t* compressedBlock = (uint32_t*)out;
	compressedBlock[0] = compressedBlock[1] = 0;
	for (int i = 0; i < 16; i++) {
		int nBlock = i < 8 ? 0 : 1;
		uint8_t cAlpha = (uint8_t) (input[i] * 255.0f);
		cAlpha = (uint8_t) (
				(cAlpha + ((cAlpha >> EXPLICIT_ALPHA_PIXEL_BPP) < 0x8 ? 7 : 8) - (cAlpha >> EXPLICIT_ALPHA_PIXEL_BPP))
						>> EXPLICIT_ALPHA_PIXEL_BPP);
		if (cAlpha > EXPLICIT_ALPHA_PIXEL_MASK)
			cAlpha = EXPLICIT_ALPHA_PIXEL_MASK;
		compressedBlock[nBlock] |= (cAlpha << ((i % 8) * EXPLICIT_ALPHA_PIXEL_BPP));
	}
}

AL2O3_EXTERN_C void Image_CompressAMDAlphaSingleModeBlock(float const input[ 4 * 4 ], void *out) {
	uint8_t nEndpoints[2][2];
	uint8_t nIndices[2][4 * 4];

	uint32_t* compressedBlock = (uint32_t*)out;
	compressedBlock[0] = compressedBlock[1] = 0;

	float fError8 = CompBlock1X(input, 4 * 4, nEndpoints[0], nIndices[0], 8, false, 8, 0, true);
	float fError6 =
			(fError8 == 0.f) ? FLT_MAX : CompBlock1X(input, 4 * 4, nEndpoints[1], nIndices[1], 6, true, 8, 0, true);
	if (fError8 <= fError6)
		EncodeAlphaBlock(compressedBlock, nEndpoints[0], nIndices[0]);
	else
		EncodeAlphaBlock(compressedBlock, nEndpoints[1], nIndices[1]);

}

void Image_CompressAMDRGBSingleModeBlock(float rgbBlock[4 * 4 * 4],
											uint32_t compressedBlock[2],
											float *pfChannelWeights,
											bool threeDRefinement,
											uint8_t refinementSteps) {

	uint8_t nEndpoints[3][2];
	uint8_t nIndices[4 * 4];

	CompRGBBlock(rgbBlock, 4 * 4,
							 RG, GG, BG,
							 nEndpoints,
							 nIndices,
							 4,
							 threeDRefinement,
							 refinementSteps,
							 pfChannelWeights,
							 false,
							 0.0f);

	unsigned int c0 = ConstructColour((nEndpoints[RC][0] >> (8 - RG)),
																		(nEndpoints[GC][0] >> (8 - GG)),
																		(nEndpoints[BC][0] >> (8 - BG)));
	unsigned int c1 = ConstructColour((nEndpoints[RC][1] >> (8 - RG)),
																		(nEndpoints[GC][1] >> (8 - GG)),
																		(nEndpoints[BC][1] >> (8 - BG)));
	if (c0 <= c1)
		compressedBlock[0] = c1 | (c0 << 16);
	else
		compressedBlock[0] = c0 | (c1 << 16);

	compressedBlock[1] = 0;
	for (int i = 0; i < 16; i++)
		compressedBlock[1] |= (nIndices[i] << (2 * i));
}

void CompressRGBABlock(float rgbaBlock[4 * 4 * 4],
											 uint32_t compressedBlock[4],
											 float *pfChannelWeights,
											 bool threeDRefinement,
											 uint8_t refinementSteps) {
	float alphaBlock[4 * 4];
	for (uint32_t i = 0; i < 16; i++) {
		alphaBlock[i] = rgbaBlock[(i * 4) + AC];
	}

	CompressAlphaBlock(alphaBlock, &compressedBlock[0]);
	CompressRGBBlock(rgbaBlock, &compressedBlock[2], pfChannelWeights, threeDRefinement, refinementSteps);
}

