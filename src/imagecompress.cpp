#include "al2o3_platform/platform.h"
#include "al2o3_memory/memory.h"
#include "gfx_image/image.hpp"
#include "gfx_imagecompress/imagecompress.h"

struct Image_Compress {
	Image_CompressFunc compressFunctions[Image_CT_MAX];

};

AL2O3_EXTERN_C Image_CompressHandle Image_CompressCreateContext() {
	auto ic = (Image_Compress *) MEMORY_CALLOC(1, sizeof(Image_Compress));
	if (!ic)
		return nullptr;

	return ic;
}

AL2O3_EXTERN_C void Image_CompressDestroyContext(Image_CompressHandle handle) {
	auto ic = (Image_Compress *) handle;
	if (!ic)
		return;
	MEMORY_FREE(ic);

}

AL2O3_EXTERN_C Image_ImageHeader const *ImageCompress_Compress(Image_CompressHandle handle,
																															 Image_CompressType type,
																															 void const *options,
																															 Image_ImageHeader const *src) {
	auto ic = (Image_Compress *) handle;
	if (!ic)
		return nullptr;

	if(ic->compressFunctions[type] == nullptr) return nullptr;

	return ic->compressFunctions[type](options, src);
}



