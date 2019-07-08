#pragma once
#ifndef GFX_IMAGECOMPRESS_IMAGECOMPRESS_H_
#define GFX_IMAGECOMPRESS_IMAGECOMPRESS_H_ 1

#include "gfx_image/image.h"

typedef struct Image_Compress *Image_CompressHandle;
typedef Image_ImageHeader const* (*Image_CompressFunc)(void const * options, Image_ImageHeader const * src);

typedef enum Image_CompressType {
	Image_CT_None = 0,

	Image_CT_BC1,
	Image_CT_BC2,
	Image_CT_BC3,
	Image_CT_BC4,
	Image_CT_BC5,
	Image_CT_BC6H,
	Image_CT_BC7,

	Image_CT_ETC_RGB,
	Image_CT_ETC2_RGB,
	Image_CT_ETC_RGBA_Explicit,
	Image_CT_ETC_RGBA_Interpolated,

	Image_CT_ASTC,

	Image_CT_MAX
} Image_CompressType;

typedef struct Image_CompressBC1Options {
	bool UseAlpha;
	bool b3DRefinement;
	uint8_t AlphaThreshold;
	uint8_t RefinementSteps;

} Image_CompressBC1Options;

AL2O3_EXTERN_C Image_CompressHandle Image_CompressCreateContext();
AL2O3_EXTERN_C void Image_CompressDestroyContext(Image_CompressHandle handle);

AL2O3_EXTERN_C Image_ImageHeader const* ImageCompress_Compress(	Image_CompressHandle handle,
																																Image_CompressType type,
																																void const * options, // null or option structur for type
																																Image_ImageHeader const * src);

AL2O3_EXTERN_C bool Image_CompressAddCompressor(Image_CompressHandle handle, Image_CompressType type, Image_CompressFunc func);


// direct access to each compressor
AL2O3_EXTERN_C Image_ImageHeader const *Image_CompressAMDBC1(Image_ImageHeader const *src, Image_CompressBC1Options const *options);

#endif // end GFX_IMAGECOMPRESS_IMAGECOMPRESS_H_