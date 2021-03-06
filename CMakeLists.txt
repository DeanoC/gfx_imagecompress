cmake_minimum_required(VERSION 3.12.4)
get_directory_property(hasParent PARENT_DIRECTORY)
if(NOT hasParent)
	option(unittests "unittests" OFF)
	get_filename_component(_PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
	set_property(GLOBAL PROPERTY GLOBAL_FETCHDEPS_BASE ${_PARENT_DIR}/al2o3 )
	include(FetchContent)
	FetchContent_Declare( al2o3 GIT_REPOSITORY https://github.com/DeanoC/al2o3 GIT_TAG master )
	FetchContent_GetProperties(al2o3)
	if(NOT al2o3_POPULATED)
		FetchContent_Populate(al2o3)
		add_subdirectory(${al2o3_SOURCE_DIR} ${al2o3_BINARY_DIR})
	endif()
	INIT_AL2O3(${CMAKE_CURRENT_SOURCE_DIR})
endif ()

SET_MIN_VERSIONS()
set(LibName gfx_imagecompress)
project(${LibName})

set(Interface
		imagecompress.h
		)

set(Src
		imagecompress.cpp

		amd_bc1_compressor.cpp
		amd_bc2_compressor.cpp
		amd_bc3_compressor.cpp
		amd_bc4_compressor.cpp
		amd_bc5_compressor.cpp
		amd_bc6h_compressor.cpp
		amd_bc7_compressor.cpp

		amd_bcx_helpers.cpp
		amd_bcx_helpers.hpp
		amd_bcx_body.cpp
		amd_bcx_body.hpp
		amd_bc6h_body.cpp
		amd_bc6h_body.hpp
		amd_bc7_body.cpp
		amd_bc7_body.hpp
		amd_hdr_encode.cpp
		amd_hdr_encode.hpp
		amd_bc7_partitions.cpp
		amd_bc7_partitions.hpp
		amd_bc7_3dquant_vpc.cpp
		amd_bc7_3dquant_vpc.hpp
		amd_shake.cpp
		amd_shake.hpp

		block_utils.hpp
		block_utils.cpp

		richgel999_bc7enc16.cpp
		richgel999_bc7enc16.h
		)

set(Deps
		al2o3_platform
		al2o3_cmath
		al2o3_stb
		tiny_imageformat
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

