// this is just a remix of AMD Code_DXT1.cpp compressor. All the hard bits are
// from AMD. See there original license at the end of this file


#include "al2o3_platform/platform.h"
#include "al2o3_cmath/scalar.h"
#include "gfx_image/image.hpp"
#include "gfx_imagecompress/imagecompress.h"
#include "gfx_image/utils.h"
#include "amd_bcx_compressor.hpp"
#include "block_utils.hpp"

void CompressRGBBlock(float rgbBlock[4 * 4 * 4],
											uint32_t compressedBlock[2],
											float pfChannelWeights[3],
											Image_CompressBC1Options const *options) {
	/*
	ARGB Channel indexes
	*/

	int RC = 2, GC = 1, BC = 0;
	/*
	Channel Bits
	*/
#define RG 5
#define GG 6
#define BG 5
#define ConstructColour(r, g, b)  (((r) << 11) | ((g) << 5) | (b))

	uint8_t nEndpoints[2][3][2];
	uint8_t nIndices[2][4 * 4];

	double fError3 = CompRGBBlock(rgbBlock, 4 * 4,
																RG, GG, BG,
																nEndpoints[0],
																nIndices[0], 3,
																options->b3DRefinement,
																options->RefinementSteps,
																pfChannelWeights,
																options->UseAlpha,
																options->AlphaThreshold);
	double fError4 = (fError3 == 0.0) ? FLT_MAX : CompRGBBlock(rgbBlock,
																														 4 * 4,
																														 RG,GG,BG,
																														 nEndpoints[1],
																														 nIndices[1],
																														 4,
																														 options->b3DRefinement,
																														 options->RefinementSteps,
																														 pfChannelWeights,
																														 options->UseAlpha,
																														 options->AlphaThreshold);

	unsigned int nMethod = (fError3 <= fError4) ? 0 : 1;
	unsigned int c0 = ConstructColour((nEndpoints[nMethod][RC][0] >> (8 - RG)),
																		(nEndpoints[nMethod][GC][0] >> (8 - GG)),
																		(nEndpoints[nMethod][BC][0] >> (8 - BG)));
	unsigned int c1 = ConstructColour((nEndpoints[nMethod][RC][1] >> (8 - RG)),
																		(nEndpoints[nMethod][GC][1] >> (8 - GG)),
																		(nEndpoints[nMethod][BC][1] >> (8 - BG)));
	if ((nMethod == 1 && c0 <= c1) || (nMethod == 0 && c0 > c1))
		compressedBlock[0] = c1 | (c0 << 16);
	else
		compressedBlock[0] = c0 | (c1 << 16);

	compressedBlock[1] = 0;
	for (int i = 0; i < 16; i++)
		compressedBlock[1] |= (nIndices[nMethod][i] << (2 * i));
}

AL2O3_EXTERN_C Image_ImageHeader const *Image_CompressAMDBC1(Image_ImageHeader const *src, Image_CompressBC1Options const *options) {
	if(src->depth > 1) { return nullptr; }

	if (options == nullptr) {
		static Image_CompressBC1Options const defaultOptions
				{
						false,
						false,
						128,
						1,
				};
		options = &defaultOptions;
	};
	bool const sRGB = ImageFormat_IsSRGB((src->format));
	bool const isAlpha = options->UseAlpha;
	ImageFormat dstFmt = sRGB ?
												isAlpha ? ImageFormat_BC1_RGBA_SRGB_BLOCK : ImageFormat_BC1_RGB_UNORM_BLOCK :
												isAlpha ? ImageFormat_BC1_RGBA_UNORM_BLOCK : ImageFormat_BC1_RGB_UNORM_BLOCK;
	Image_ImageHeader const * dst = Image_CreateNoClear(src->width, src->height, 1, src->slices, dstFmt);
	if(!dst) return nullptr;

	// get block size round up to 4
	size_t const blocksX = (src->width + 3) / 4;
	size_t const blocksY = (src->height + 3) / 4;

	float fAlphaThreshold = ((float) options->AlphaThreshold) / 255.f;

	for (size_t w = 0; w < src->slices; ++w) {
		for (size_t y = 0; y < blocksY; ++y) {
			for (size_t x = 0; x < blocksX; ++x) {

				uint32_t compressedBlock[2];

				float srcBlock[4 * 4 * 4];
				float weights[3];

				ImageCompress::ReadNxNBlock(src, 4, 4, srcBlock, x * 4, y * 4, w);
				ImageCompress::CalculateColourWeightings(srcBlock, weights);

				CompressRGBBlock(srcBlock,
												 compressedBlock,
												 weights,
												 options);
				ImageCompress::WriteNxNBlock(dst, 4, 4, compressedBlock, 2 * sizeof(uint32_t), x * 4, y * 4, w);

				//			bufferOut.WriteBlock(i * 4, j * 4, compressedBlock, 2);
			}
		}
		//		if (pFeedbackProc) {
		//			float fProgress = 100.f * (j * dwBlocksX) / (dwBlocksX * dwBlocksY);
		//			if (pFeedbackProc(fProgress, pUser1, pUser2))
		//				return CE_Aborted;
		//		}
	}
	return dst;
}


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
//  File Name:   Codec_DXT1.cpp
//  Description: implementation of the CCodec_DXT1 class
//
//////////////////////////////////////////////////////////////////////////////
