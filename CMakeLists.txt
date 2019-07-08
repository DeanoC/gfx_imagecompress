#if building standalone and not specified build tests
option(BUILD_TESTS "Whether to build tests" ON)
include(FetchContent)
FetchContent_Declare( al2o3 GIT_REPOSITORY https://github.com/DeanoC/al2o3 GIT_TAG master )
FetchContent_GetProperties(al2o3)
if(NOT al2o3_POPULATED)
	FetchContent_Populate(al2o3)
	add_subdirectory(${al2o3_SOURCE_DIR} ${al2o3_BINARY_DIR})
endif()

SET_MIN_VERSIONS()
set(LibName gfx_imagecompress)
project(${LibName})

set(Interface
		imagecompress.h
		)

set(Src
		imagecompress.cpp
		amd_bcx_compressor.cpp
		amd_bcx_compressor.hpp
		amd_bc1_normalspeed.cpp
		block_utils.hpp
		block_utils.cpp
		)

set(Deps
		al2o3_platform
		al2o3_cmath
		al2o3_stb
		gfx_imageformat
		gfx_image_interface
		utils_misccpp
		)
ADD_LIB(${LibName} "${Interface}" "${Src}" "${Deps}")

set( Tests
		runner.cpp
		test_imagecompress.cpp
		)
set( TestDeps
		utils_simple_logmanager
		gfx_image_impl_basic
		gfx_imageio
		al2o3_catch2 )
ADD_LIB_TESTS(${LibName} "${Interface}" "${Tests}" "${TestDeps}")

