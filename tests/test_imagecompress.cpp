#include <al2o3_vfile/vfile.hpp>
#include "al2o3_platform/platform.h"
#include "al2o3_catch2/catch2.hpp"
#include "gfx_image/create.h"
#include "gfx_imagecompress/imagecompress.h"
#include "gfx_imageio/io.h"


#define LOAD_DDS(filename, img) { VFile::ScopedFile fh = VFile::File::FromFile("artifacts/" filename, Os_FM_ReadBinary); \
	REQUIRE(fh); img = Image_LoadDDS(fh); }
#define SAVE_DDS(img, filename) { VFile::ScopedFile fh = VFile::File::FromFile("artifacts/" filename, Os_FM_WriteBinary); \
	REQUIRE(fh); Image_SaveDDS(img, fh); }

static void GenerateTestPatternRGB(Image_ImageHeader const* image) {
	Image_PixelD const red = { 1, 0, 0, 1 };
	Image_PixelD const green = { 0, 1, 0, 1 };
	Image_PixelD const blue = { 0, 0, 1, 1 };

	for(auto y = 0u;y < image->height;++y) {
		for(auto x = 0u;x < image->width;++x) {
			size_t index = Image_CalculateIndex(image, x, y, 0, 0);
			Image_PixelD col = green;

			if(((x / 2) & 2) || ((y / 2) & 2)) col = red;
			if(((x / 3) % 3) && ((y / 3) % 3)) col = blue;

			Image_SetPixelAt(image, &col, index);
		}
	}
}

static void GenerateTestPatternRGB_Punchthrough(Image_ImageHeader const* image) {
	Image_PixelD const red = { 1, 0, 0, 1 };
	Image_PixelD const green = { 0, 1, 0, 1 };
	Image_PixelD const blue = { 0, 0, 1, 1 };
	Image_PixelD const nada = { 0, 0, 0, 0 };

	for(auto y = 0u;y < image->height;++y) {
		for(auto x = 0u;x < image->width;++x) {
			size_t index = Image_CalculateIndex(image, x, y, 0, 0);
			Image_PixelD col = green;

			if(((x / 2) & 2) || ((y / 2) & 2)) col = red;
			if(((x / 3) % 3) && ((y / 3) % 3)) col = blue;
			if(x > image->width/2 && y > image->height/2) {col = nada; }

			Image_SetPixelAt(image, &col, index);
		}
	}
}

static void GenerateTestPatternRGBA(Image_ImageHeader const* image) {
	Image_PixelD const red = { 1, 0, 0, 1 };
	Image_PixelD const green = { 0, 1, 0, 1 };
	Image_PixelD const blue = { 0, 0, 1, 1 };

	for(auto y = 0u;y < image->height;++y) {
		for(auto x = 0u;x < image->width;++x) {
			size_t index = Image_CalculateIndex(image, x, y, 0, 0);
			Image_PixelD col = green;

			if(((x / 2) & 2) || ((y / 2) & 2)) col = red;
			if(((x / 3) % 3) && ((y / 3) % 3)) col = blue;

			col.a	= (float)x / (float)image->width;

			Image_SetPixelAt(image, &col, index);
		}
	}
}

static void GenerateTestPatternFloatRGBA(Image_ImageHeader const* image) {
	Image_PixelD const red = { 1, 0, 0, 1 };
	Image_PixelD const green = { 0, 1, 0, 1 };
	Image_PixelD const blue = { 0, 0, 1, 1 };

	for(auto y = 0u;y < image->height;++y) {
		for(auto x = 0u;x < image->width;++x) {
			size_t index = Image_CalculateIndex(image, x, y, 0, 0);
			Image_PixelD col = green;

			if(((x / 2) & 2) || ((y / 2) & 2)) col = red;
			if(((x / 3) % 3) && ((y / 3) % 3)) col = blue;

			col.a	= (float)x / (float)image->width;

			Image_SetPixelAt(image, &col, index);
		}
	}
}

TEST_CASE("AMD BC1 Direct RGB", "[ImageCompress]") {
	auto image = Image_Create2D(256, 256, ImageFormat_R8G8B8_UNORM);
	REQUIRE(image);
	GenerateTestPatternRGB(image);
	SAVE_DDS(image, "compress_ref_R8G8B8_UNORM_256x256.dds")

	auto dst = Image_CompressAMDBC1(image, nullptr, nullptr, nullptr, nullptr);

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
	GenerateTestPatternRGB(image);
	SAVE_DDS(image, "compress_ref_R8G8B8_UNORM_257x257.dds")

	auto dst = Image_CompressAMDBC1(image, nullptr, nullptr, nullptr, nullptr);

	REQUIRE(dst);
	REQUIRE(dst->width == 260);
	REQUIRE(dst->height == 260);
	REQUIRE(dst->depth == image->depth);
	REQUIRE(dst->slices == image->slices);
	REQUIRE(dst->format == ImageFormat_BC1_RGB_UNORM_BLOCK);

	SAVE_DDS(dst, "compress_BC1_RGB_UNORM_260x260.dds")

	Image_Destroy(image);
	Image_Destroy(dst);
}


TEST_CASE("AMD BC1 Direct RGB Alpha options src no alpha", "[ImageCompress]") {
	Image_ImageHeader const* image = nullptr;
	LOAD_DDS("compress_ref_R8G8B8_UNORM_256x256.dds", image)
	REQUIRE(image);

	Image_CompressBC1Options options {};
	options.UseAlpha = true;
	options.AlphaThreshold = 128;

	auto dst = Image_CompressAMDBC1(image, nullptr, &options, nullptr, nullptr);

	REQUIRE(dst);
	REQUIRE(dst->width == image->width);
	REQUIRE(dst->height == image->height);
	REQUIRE(dst->depth == image->depth);
	REQUIRE(dst->slices == image->slices);
	REQUIRE(dst->format == ImageFormat_BC1_RGBA_UNORM_BLOCK);

	SAVE_DDS(dst, "compress_BC1_RGBA_UNORM_256x256.dds")

	Image_Destroy(image);
	Image_Destroy(dst);
}

TEST_CASE("AMD BC1 Direct RGB Alpha options src with alpha", "[ImageCompress]") {
	auto image = Image_Create2D(256, 256, ImageFormat_R8G8B8A8_UNORM);
	REQUIRE(image);
	GenerateTestPatternRGB_Punchthrough(image);
	SAVE_DDS(image, "compress_ref_punchthrough_R8G8B8A8_UNORM_256x256.dds")


	Image_CompressBC1Options options {};
	options.UseAlpha = true;
	options.AlphaThreshold = 128;

	auto dst = Image_CompressAMDBC1(image, nullptr, &options, nullptr, nullptr);

	REQUIRE(dst);
	REQUIRE(dst->width == image->width);
	REQUIRE(dst->height == image->height);
	REQUIRE(dst->depth == image->depth);
	REQUIRE(dst->slices == image->slices);
	REQUIRE(dst->format == ImageFormat_BC1_RGBA_UNORM_BLOCK);

	SAVE_DDS(dst, "compress_punchthrough_BC1_RGBA_UNORM_256x256.dds")

	Image_Destroy(image);
	Image_Destroy(dst);
}
TEST_CASE("AMD BC2 Direct RGBA", "[ImageCompress]") {
	auto image = Image_Create2D(256, 256, ImageFormat_R8G8B8A8_UNORM);
	REQUIRE(image);
	GenerateTestPatternRGBA(image);
	SAVE_DDS(image, "compress_ref_R8G8B8A8_UNORM_256x256.dds")


	auto dst = Image_CompressAMDBC2(image, nullptr, nullptr, nullptr);

	REQUIRE(dst);
	REQUIRE(dst->width == image->width);
	REQUIRE(dst->height == image->height);
	REQUIRE(dst->depth == image->depth);
	REQUIRE(dst->slices == image->slices);
	REQUIRE(dst->format == ImageFormat_BC2_UNORM_BLOCK);

	SAVE_DDS(dst, "compress_BC2_UNORM_256x256.dds")

	Image_Destroy(image);
	Image_Destroy(dst);
}

TEST_CASE("AMD BC3 Direct RGBA", "[ImageCompress]") {
	auto image = Image_Create2D(256, 256, ImageFormat_R8G8B8A8_UNORM);
	REQUIRE(image);
	GenerateTestPatternRGBA(image);
	SAVE_DDS(image, "compress_ref_R8G8B8A8_UNORM_256x256.dds")

	auto dst = Image_CompressAMDBC3(image, nullptr, nullptr, nullptr);

	REQUIRE(dst);
	REQUIRE(dst->width == image->width);
	REQUIRE(dst->height == image->height);
	REQUIRE(dst->depth == image->depth);
	REQUIRE(dst->slices == image->slices);
	REQUIRE(dst->format == ImageFormat_BC3_UNORM_BLOCK);

	SAVE_DDS(dst, "compress_BC3_UNORM_256x256.dds")

	Image_Destroy(image);
	Image_Destroy(dst);
}

TEST_CASE("AMD BC4 Direct R", "[ImageCompress]") {
	auto image = Image_Create2D(256, 256, ImageFormat_R8G8B8A8_UNORM);
	REQUIRE(image);
	GenerateTestPatternRGBA(image);
	SAVE_DDS(image, "compress_ref_R8G8B8A8_UNORM_256x256.dds")

	auto dst = Image_CompressAMDBC4(image, nullptr, nullptr);

	REQUIRE(dst);
	REQUIRE(dst->width == image->width);
	REQUIRE(dst->height == image->height);
	REQUIRE(dst->depth == image->depth);
	REQUIRE(dst->slices == image->slices);
	REQUIRE(dst->format == ImageFormat_BC4_UNORM_BLOCK);

	SAVE_DDS(dst, "compress_BC4_UNORM_256x256.dds")

	Image_Destroy(image);
	Image_Destroy(dst);
}
TEST_CASE("AMD BC5 Direct RG", "[ImageCompress]") {
	auto image = Image_Create2D(256, 256, ImageFormat_R8G8B8A8_UNORM);
	REQUIRE(image);
	GenerateTestPatternRGBA(image);
	SAVE_DDS(image, "compress_ref_R8G8B8A8_UNORM_256x256.dds")

	auto dst = Image_CompressAMDBC5(image, nullptr, nullptr);

	REQUIRE(dst);
	REQUIRE(dst->width == image->width);
	REQUIRE(dst->height == image->height);
	REQUIRE(dst->depth == image->depth);
	REQUIRE(dst->slices == image->slices);
	REQUIRE(dst->format == ImageFormat_BC5_UNORM_BLOCK);

	SAVE_DDS(dst, "compress_BC5_UNORM_256x256.dds")

	Image_Destroy(image);
	Image_Destroy(dst);
}
/*
TEST_CASE("AMD BC6H Direct RGBA", "[ImageCompress]") {
	auto image = Image_Create2D(32, 32, ImageFormat_R32G32B32A32_SFLOAT);
	REQUIRE(image);
	GenerateTestPatternFloatRGBA(image);
	SAVE_DDS(image, "compress_ref_R8G8B8A8_SFLOAT_32x32.dds")

	auto dst = Image_CompressAMDBC6H(image, nullptr, nullptr, nullptr);

	REQUIRE(dst);
	REQUIRE(dst->width == image->width);
	REQUIRE(dst->height == image->height);
	REQUIRE(dst->depth == image->depth);
	REQUIRE(dst->slices == image->slices);
	REQUIRE(dst->format == ImageFormat_BC6H_SFLOAT_BLOCK);

	SAVE_DDS(dst, "compress_BC6H_SFLOAT_32x32.dds")

	Image_Destroy(image);
	Image_Destroy(dst);
}
*/
TEST_CASE("AMD BC7 Direct RGB", "[ImageCompress]") {
	auto image = Image_Create2D(32, 32, ImageFormat_R8G8B8_UNORM);
	REQUIRE(image);
	GenerateTestPatternRGB(image);
	SAVE_DDS(image, "compress_ref_R8G8B8_UNORM_32x32.dds")

	auto dst = Image_CompressAMDBC7(image, nullptr, nullptr, nullptr);

	REQUIRE(dst);
	REQUIRE(dst->width == image->width);
	REQUIRE(dst->height == image->height);
	REQUIRE(dst->depth == image->depth);
	REQUIRE(dst->slices == image->slices);
	REQUIRE(dst->format == ImageFormat_BC7_UNORM_BLOCK);

	SAVE_DDS(dst, "compress_BC7_UNORM_32x32.dds")

	Image_Destroy(image);
	Image_Destroy(dst);
}