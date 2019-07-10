#pragma once
#ifndef GFX_IMAGECOMPRESS_IMAGECOMPRESS_H_
#define GFX_IMAGECOMPRESS_IMAGECOMPRESS_H_ 1

#include "gfx_image/image.h"

typedef struct Image_Compress *Image_CompressHandle;
typedef Image_ImageHeader const* (*Image_CompressFunc)(void const * options, Image_ImageHeader const * src);

typedef bool (*Image_CompressProgressFunc)(void* user, float percentage);

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
	bool UseAlpha;						// default false
	uint8_t AlphaThreshold;		// default 128
} Image_CompressBC1Options;

typedef struct Image_CompressAMDBackendOptions {
	bool b3DRefinement;				// default false
	bool AdaptiveColourWeights;	// default false (NOTE not sure this is working yet)
	uint8_t RefinementSteps; 	// default 1
	uint8_t ModeMask;					// default 0xFF used by BC6H and BC7
} Image_CompressAMDBackendOptions;

AL2O3_EXTERN_C Image_CompressHandle Image_CompressCreateContext();
AL2O3_EXTERN_C void Image_CompressDestroyContext(Image_CompressHandle handle);

AL2O3_EXTERN_C Image_ImageHeader const* ImageCompress_Compress(	Image_CompressHandle handle,
																																Image_CompressType type,
																																void const * options, // null or option structure
																																Image_ImageHeader const * src);

AL2O3_EXTERN_C bool Image_CompressAddCompressor(Image_CompressHandle handle, Image_CompressType type, Image_CompressFunc func);


// direct access to each compressor options can be null
AL2O3_EXTERN_C Image_ImageHeader const *Image_CompressAMDBC1(Image_ImageHeader const *src, Image_CompressAMDBackendOptions const* amdOptions, Image_CompressBC1Options const *options, Image_CompressProgressFunc progressCallback, void* userCallbackData);
AL2O3_EXTERN_C Image_ImageHeader const *Image_CompressAMDBC2(Image_ImageHeader const *src, Image_CompressAMDBackendOptions const* amdOptions, Image_CompressProgressFunc progressCallback, void* userCallbackData);
AL2O3_EXTERN_C Image_ImageHeader const *Image_CompressAMDBC3(Image_ImageHeader const *src, Image_CompressAMDBackendOptions const* amdOptions, Image_CompressProgressFunc progressCallback, void* userCallbackData);
AL2O3_EXTERN_C Image_ImageHeader const *Image_CompressAMDBC4(Image_ImageHeader const *src, Image_CompressProgressFunc progressCallback, void* userCallbackData);
AL2O3_EXTERN_C Image_ImageHeader const *Image_CompressAMDBC5(Image_ImageHeader const *src, Image_CompressProgressFunc progressCallback, void* userCallbackData);
AL2O3_EXTERN_C Image_ImageHeader const *Image_CompressAMDBC6H(Image_ImageHeader const *src, Image_CompressAMDBackendOptions const* amdOptions, Image_CompressProgressFunc progressCallback, void* userCallbackData);
AL2O3_EXTERN_C Image_ImageHeader const *Image_CompressAMDBC7(Image_ImageHeader const *src, Image_CompressAMDBackendOptions const* amdOptions, Image_CompressProgressFunc progressCallback, void* userCallbackData);

#endif // end GFX_IMAGECOMPRESS_IMAGECOMPRESS_H_