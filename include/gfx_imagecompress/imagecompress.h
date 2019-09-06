#pragma once

#include "gfx_image/image.h"

typedef bool (*Image_CompressProgressFunc)(void *user, float percentage);

typedef enum Image_CompressType {
	Image_CT_None = 0,

	Image_CT_DXBC1,
	Image_CT_DXBC2,
	Image_CT_DXBC3,
	Image_CT_DXBC4,
	Image_CT_DXBC5,
	Image_CT_DXBC6H,
	Image_CT_DXBC7,

	Image_CT_ETC_RGB,
	Image_CT_ETC2_RGB,
	Image_CT_ETC_RGBA_Explicit,
	Image_CT_ETC_RGBA_Interpolated,

	Image_CT_ASTC,

	Image_CT_MAX
} Image_CompressType;

typedef enum Image_CompressPickFlags {
	Image_CPF_AllowDXBC1to5 = 0x1,
	Image_CPF_AllowASTC = 0x2,
	Image_CPF_AllowETC = 0x8,
	Image_CPF_AllowDXBC6and7 = 0x10
} Image_CompressPickFlags;

typedef struct Image_CompressBC1Options {
	bool UseAlpha;            // default false
	uint8_t AlphaThreshold;    // default 128
} Image_CompressBC1Options;

typedef struct Image_CompressAMDBackendOptions {
	bool b3DRefinement;        // default false
	bool AdaptiveColourWeights;  // default false (NOTE not sure this is working yet)
	uint8_t RefinementSteps;  // default 1
	uint8_t ModeMask;          // default 0xFF used by BC6H and BC7
} Image_CompressAMDBackendOptions;

typedef struct Image_CompressRichGel99BackendOptions {
	bool perceptual;
	bool fast;
} Image_CompressRichGel999BackendOptions;

// some compression system (BC7 at the moment) create global tables
// these work as ref counters for these tables. Generally only need to be
// called if using the loweset level block API. Higher apis call them around
// each compress image call. However if doing many may be more performant to
// init/deinit around a block of compress calls
AL2O3_EXTERN_C void Image_CompressInit();
AL2O3_EXTERN_C void Image_CompressDeinit();

AL2O3_EXTERN_C Image_ImageHeader const *ImageCompress_Compress(Image_CompressType type,
																															 bool fast,
																															 Image_ImageHeader const *src);

AL2O3_EXTERN_C Image_CompressType ImageCompress_PickCompressionType(Image_CompressPickFlags flags,
																																		Image_ImageHeader const *src);


// direct access to each compressor options can be null
AL2O3_EXTERN_C Image_ImageHeader const *Image_CompressAMDBC1(Image_ImageHeader const *src,
																														 Image_CompressAMDBackendOptions const *amdOptions,
																														 Image_CompressBC1Options const *options,
																														 Image_CompressProgressFunc progressCallback,
																														 void *userCallbackData);
AL2O3_EXTERN_C Image_ImageHeader const *Image_CompressAMDBC2(Image_ImageHeader const *src,
																														 Image_CompressAMDBackendOptions const *amdOptions,
																														 Image_CompressProgressFunc progressCallback,
																														 void *userCallbackData);
AL2O3_EXTERN_C Image_ImageHeader const *Image_CompressAMDBC3(Image_ImageHeader const *src,
																														 Image_CompressAMDBackendOptions const *amdOptions,
																														 Image_CompressProgressFunc progressCallback,
																														 void *userCallbackData);
AL2O3_EXTERN_C Image_ImageHeader const *Image_CompressAMDBC4(Image_ImageHeader const *src,
																														 Image_CompressProgressFunc progressCallback,
																														 void *userCallbackData);
AL2O3_EXTERN_C Image_ImageHeader const *Image_CompressAMDBC5(Image_ImageHeader const *src,
																														 Image_CompressProgressFunc progressCallback,
																														 void *userCallbackData);
AL2O3_EXTERN_C Image_ImageHeader const *Image_CompressAMDBC6H(Image_ImageHeader const *src,
																															Image_CompressAMDBackendOptions const *amdOptions,
																															Image_CompressProgressFunc progressCallback,
																															void *userCallbackData);
AL2O3_EXTERN_C Image_ImageHeader const *Image_CompressAMDBC7(Image_ImageHeader const *src,
																														 Image_CompressAMDBackendOptions const *amdOptions,
																														 Image_CompressProgressFunc progressCallback,
																														 void *userCallbackData);

AL2O3_EXTERN_C Image_ImageHeader const *Image_CompressRichGel999BC7(Image_ImageHeader const *src,
																																		Image_CompressRichGel999BackendOptions const *richOptions,
																																		Image_CompressProgressFunc progressCallback,
																																		void *userCallbackData);


// lowest level interface block compression API
// inputs are normalised float values (0-1 for LDR) in most cases
// exceptions being HDR and SNORM modes

// block using by 'single mode' DXTC (BC1-5) RGB compression
// these 4 block make up BC 1 to 5

// takes 16 x RGB (3) float input, outputs 8 byte RGB block
AL2O3_EXTERN_C void Image_CompressAMDRGBSingleModeBlock(float const input[4 * 4 * 3],
																												bool adaptiveColourWeights,
																												bool b3DRefinement,
																												uint8_t refinementSteps,
																												void *out);
// takes 16 x A float input, outputs 8 byte Alpha block
AL2O3_EXTERN_C void Image_CompressAMDAlphaSingleModeBlock(float const input[4 * 4], void *out);
// takes 16 x A float input, outputs 8 byte Explicit 4 bit Alpha block
AL2O3_EXTERN_C void Image_CompressAMDExplictAlphaSingleModeBlock(float const input[4 * 4], void *out);

// special case BC1 compress 4x4 (4 channel) to 8 byte BC1 RGB block threshold is 0 - 1 if <= 0 no alpha
AL2O3_EXTERN_C void Image_CompressAMDBC1Block(float const input[4 * 4 * 4],
																							bool adaptiveColourWeight,
																							bool b3DRefinement,
																							uint8_t refinementSteps,
																							float alphaThreshold,
																							void *out);
// Takes 16 x RGBA (4) float inputs and outputs 16 byte multi mode block
AL2O3_EXTERN_C void Image_CompressAMDMultiModeLDRBlock(float const input[4 * 4 * 4],
																											 uint8_t modeMask,
																											 bool srcHasAlpha,
																											 float quality,
																											 bool colourRestrict,
																											 bool alphaRestrict,
																											 float performance,
																											 void *out);

// takes 16 RGBA uint32_t inputs an outputs 16 byte multi mode (only 1 or 6) block
AL2O3_EXTERN_C void Image_CompressRichGel999BC7enc16(uint32_t const input[4 * 4],
																										 bool fast,
																										 bool perceptual,
																										 void *out);