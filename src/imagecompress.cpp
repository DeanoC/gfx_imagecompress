#include "al2o3_platform/platform.h"
#include "al2o3_memory/memory.h"
#include "gfx_image/image.hpp"
#include "gfx_imagecompress/imagecompress.h"

AL2O3_EXTERN_C Image_ImageHeader const *ImageCompress_Compress(Image_CompressType type,
																															 bool fast,
																															 Image_ImageHeader const *src) {
	switch (type) {
	case Image_CT_None: return src;

	case Image_CT_DXBC1: return Image_CompressAMDBC1(src, nullptr, nullptr, nullptr, nullptr);
	case Image_CT_DXBC2: return Image_CompressAMDBC2(src, nullptr, nullptr, nullptr);
	case Image_CT_DXBC3: return Image_CompressAMDBC3(src, nullptr, nullptr, nullptr);
	case Image_CT_DXBC4: return Image_CompressAMDBC4(src, nullptr, nullptr);
	case Image_CT_DXBC5: return Image_CompressAMDBC5(src, nullptr, nullptr);
	case Image_CT_DXBC6H: return Image_CompressAMDBC6H(src, nullptr, nullptr, nullptr);
	case Image_CT_DXBC7:
		if (fast) {
			return Image_CompressRichGel999BC7(src, nullptr, nullptr, nullptr);
		} else {
			return Image_CompressAMDBC7(src, nullptr, nullptr, nullptr);
		}
	case Image_CT_ETC_RGB:
	case Image_CT_ETC2_RGB:
	case Image_CT_ETC_RGBA_Explicit:
	case Image_CT_ETC_RGBA_Interpolated:
	case Image_CT_ASTC: return nullptr;
	default: ASSERT(false);
		return nullptr;
	}
}

AL2O3_EXTERN_C Image_CompressType ImageCompress_PickCompressionType(Image_CompressPickFlags flags,
																																		Image_ImageHeader const *src) {

	// only HDR codec we have is BC6H so advise no compression for floats unless thats supported
	if (TinyImageFormat_IsFloat(src->format)) {
		if ((flags & Image_CPF_AllowDXBC6and7) == 0) {
			return Image_CT_None;
		}
	} else if (!TinyImageFormat_IsNormalised(src->format)) {
		// except for floats only normalised formats have compression types
		return Image_CT_None;
	}

	bool hasAlpha = false;

	// how many channels
	// TODO SNORM vs UNORM
	switch (TinyImageFormat_ChannelCount(src->format)) {
	case 1:
		if (flags & Image_CPF_AllowDXBC1to5) {
			return Image_CT_DXBC4;
		};
		// TODO ETC2 single channel format
		break;
	case 2:
		if (flags & Image_CPF_AllowDXBC1to5) {
			return Image_CT_DXBC5;
		}
		// TODO ETC2 dual channel format
		break;
	case 3: break;
	case 4: hasAlpha = true;
		break;
	}

	// which RGB or RGBA formats should we use?
	if (flags & Image_CPF_AllowDXBC6and7) {
		// TODO in some cases BC1 might be preferable (half size)...
		return Image_CT_DXBC7; // if bc7 is avail use it
	}

	if (flags & Image_CPF_AllowASTC) {
		return Image_CT_ASTC; // if astc is avail use it
	}
	if (hasAlpha) {
		if (flags & Image_CPF_AllowDXBC1to5) {
			return Image_CT_DXBC3;
		}
		if (flags & Image_CPF_AllowETC) {
			return Image_CT_None;
			//			return Image_CT_ETC_RGBA_Interpolated;
		}
	} else {
		if (flags & Image_CPF_AllowDXBC1to5) {
			return Image_CT_DXBC1;
		}
		if (flags & Image_CPF_AllowETC) {
			return Image_CT_None;
			//			return Image_CT_ETC2_RGB;
		}
	}

	// no available compression format to pick from
	return Image_CT_None;
}


