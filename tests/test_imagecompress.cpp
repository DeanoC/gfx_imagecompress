#include <al2o3_vfile/vfile.hpp>
#include "al2o3_platform/platform.h"
#include "al2o3_catch2/catch2.hpp"
#include "gfx_image/create.h"
#include "gfx_imagecompress/imagecompress.h"
#include "gfx_imageio/io.h"


TEST_CASE("AMD BC1 Direct RGB", "[ImageCompress]") {
	auto image = Image_Create2D(256, 256, ImageFormat_R8G8B8_UNORM);
	REQUIRE(image);

	auto dst = Image_CompressAMDBC1(image, nullptr);

	REQUIRE(dst);
	REQUIRE(dst->width == image->width);
	REQUIRE(dst->height == image->height);
	REQUIRE(dst->depth == image->depth);
	REQUIRE(dst->slices == image->slices);
	REQUIRE(dst->format == ImageFormat_BC1_RGB_UNORM_BLOCK);
	Image_Destroy(image);
	Image_Destroy(dst);
}
TEST_CASE("AMD BC1 Direct RGB npot", "[ImageCompress]") {
	auto image = Image_Create2D(257, 257, ImageFormat_R8G8B8_UNORM);
	REQUIRE(image);

	Image_PixelD const red = { 1, 0, 0, 1 };
	Image_PixelD const green = { 0, 1, 0, 1 };

	for(auto y = 0u;y < 257;++y) {
		for(auto x = 0u;x < 257;++x) {
			size_t index = Image_CalculateIndex(image, x, y, 0, 0);
			if((x / 2) & 2) Image_SetPixelAt(image, &red, index);
			else Image_SetPixelAt(image, &green, index);
		}
	}

	VFile::ScopedFile fhi = VFile::File::FromFile("image.dds", Os_FM_WriteBinary);
	Image_SaveDDS(image, fhi);

	auto dst = Image_CompressAMDBC1(image, nullptr);

	REQUIRE(dst);
	REQUIRE(dst->width == 260);
	REQUIRE(dst->height == 260);
	REQUIRE(dst->depth == image->depth);
	REQUIRE(dst->slices == image->slices);
	REQUIRE(dst->format == ImageFormat_BC1_RGB_UNORM_BLOCK);

	VFile::ScopedFile fh = VFile::File::FromFile("test.dds", Os_FM_WriteBinary);
	Image_SaveDDS(dst, fh);

	Image_Destroy(image);
	Image_Destroy(dst);
}