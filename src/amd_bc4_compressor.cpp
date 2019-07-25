// this is just a remix of AMD Code_DXT1.cpp compressor. All the hard bits are
// from AMD. See there original license at the end of this file


#include "al2o3_platform/platform.h"
#include "gfx_image/image.hpp"
#include "gfx_imagecompress/imagecompress.h"
#include "block_utils.hpp"
#include "amd_bcx_helpers.hpp"

AL2O3_EXTERN_C Image_ImageHeader const *Image_CompressAMDBC4(Image_ImageHeader const *src,
																														 Image_CompressProgressFunc progressCallback,
																														 void *userCallbackData) {
	if (src->depth > 1) return nullptr;

	bool const srcIsSigned = TinyImageFormat_IsSigned(src->format);
	bool const srcHasAlpha = TinyImageFormat_ChannelCount(src->format) > 3;

	TinyImageFormat dstFmt = srcIsSigned ? TinyImageFormat_BC4_SNORM_BLOCK : TinyImageFormat_BC4_UNORM_BLOCK;
	Image_ImageHeader const *dst = Image_CreateNoClear(src->width, src->height, 1, src->slices, dstFmt);
	if (!dst) return nullptr;

	// get block size round up to 4
	size_t const blocksX = (src->width + 3) / 4;
	size_t const blocksY = (src->height + 3) / 4;

	for (uint32_t w = 0; w < src->slices; ++w) {
		for (uint32_t y = 0; y < blocksY; ++y) {
			for (uint32_t x = 0; x < blocksX; ++x) {
				uint32_t compressedBlock[2];

				float srcBlock[4 * 4 * 4];

				ImageCompress::ReadNxNBlock(src, 4, 4,
																		!srcHasAlpha, srcBlock, x * 4, y * 4, w);
				// extract just the red
				float redBlock[4 * 4];
				for (uint32_t i = 0; i < 16; i++) {
					redBlock[i] = srcBlock[(i * 4) + 0];
				}

				CompressAlphaBlock(redBlock, compressedBlock);

				ImageCompress::WriteNxNBlock(dst, 4, 4,
																		 compressedBlock, sizeof(compressedBlock),
																		 x * 4, y * 4, w);
			}
			if (progressCallback) {
				float fProgress = 100.f * (y * blocksX) / (blocksX * blocksY);
				if (progressCallback(userCallbackData, fProgress))
					return nullptr;
			}
		}
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
//  File Name:   Codec_DXT3.cpp
//  Description: implementation of the CCodec_DXT3 class
//
//////////////////////////////////////////////////////////////////////////////


