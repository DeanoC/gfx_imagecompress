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
//  BC7_Encode.cpp : A reference encoder for BC7
//

#include "al2o3_platform/platform.h"
#include "al2o3_cmath/scalar.h"
#include "amd_bc7_body.hpp"
#include "amd_hdr_encode.hpp"
#include "amd_bc7_partitions.hpp"
#include "amd_bc7_3dquant_vpc.hpp"
#include "amd_shake.hpp"
#include <float.h>

#define TRUE 1
#define FALSE 0

// Largest possible size for an individual subset
#define MAX_SUBSET_SIZE         16

// Maximum number of possible subsets
#define MAX_SUBSETS             3

// Maximum number of index bits
#define MAX_INDEX_BITS          4

// Maximum number of partition types
#define MAX_PARTITIONS          64

// Number of block types in the format
#define NUM_BLOCK_TYPES         8

// Size of a compressed block in bytes
#define COMPRESSED_BLOCK_SIZE   16

// If this define is set then 6-bit weights will be used for the ramp.
// Otherwise the ramp will use a pure linear interpolation
#define USE_FINAL_BC7_WEIGHTS   1

// If this is defined, ramp calculation is done via math floor and division.
// Otherwise, ramp calculation is done by bit shifting
#define USE_HIGH_PRECISION_INTERPOLATION_BC7

#define MAX_PARTITIONS_TABLE (1+64+64)

#define MAX_ENTRIES_QUANT_TRACE     16
#define MAX_CLUSTERS_QUANT_TRACE    8

typedef enum _COMPONENT {
	COMP_RED = 0,
	COMP_GREEN = 1,
	COMP_BLUE = 2,
	COMP_ALPHA = 3
} COMPONENT;

//
// Block encoding information for all block types
// {Component Encoding, PartitionBits, RotationBits, indexSwapBits,
//  scalarBits, vectorBits, pBitType, subsetCount, {index0Bits, index1Bits}}
//
CMP_BTI bti[NUM_BLOCK_TYPES] =
		{
				{NO_ALPHA, 4, 0, 0, 0, 12, TWO_PBIT, 3, {3, 0}},  // Format Mode 0
				{NO_ALPHA, 6, 0, 0, 0, 18, ONE_PBIT, 2, {3, 0}},  // Format Mode 1
				{NO_ALPHA, 6, 0, 0, 0, 15, NO_PBIT, 3, {2, 0}},  // Format Mode 2
				{NO_ALPHA, 6, 0, 0, 0, 21, TWO_PBIT, 2, {2, 0}},  // Format Mode 3
				{SEPARATE_ALPHA, 0, 2, 1, 6, 15, NO_PBIT, 1, {2, 3}},  // Format Mode 4
				{SEPARATE_ALPHA, 0, 2, 0, 8, 21, NO_PBIT, 1, {2, 2}},  // Format Mode 5
				{COMBINED_ALPHA, 0, 0, 0, 0, 28, TWO_PBIT, 1, {4, 0}},  // Format Mode 6
				{COMBINED_ALPHA, 6, 0, 0, 0, 20, TWO_PBIT, 2, {2, 0}}   // Format Mode 7
		};

// Write one bit to a buffer at a specified offset
//
// Used by BC7_Encode

void WriteBit(uint8_t *base,
							int offset,
							uint8_t bitVal) {
	int byteLocation;
	int remainder;
	uint8_t val;
	uint8_t mask;

	byteLocation = offset / 8;
	remainder = offset % 8;
	mask = ~(1 << remainder);
	bitVal &= 1;

	val = base[byteLocation];
	val = val & mask;
	val |= bitVal << remainder;

	base[byteLocation] = val & 0xff;
}

// Used by BC7_Decode
const float rampLerpWeights[5][1 << MAX_INDEX_BITS] =
		{
#if USE_FINAL_BC7_WEIGHTS
				{0.0},  // 0 bit index
				{0.0, 1.0}, // 1 bit index
				{0.0, 21.0 / 64.0, 43.0 / 64.0, 1.0}, // 2 bit index
				{0.0, 9.0 / 64.0, 18.0 / 64.0, 27.0 / 64.0, 37.0 / 64.0, 46.0 / 64.0, 55.0 / 64.0, 1.0}, // 3 bit index
				{0.0, 4.0 / 64.0, 9.0 / 64.0, 13.0 / 64.0, 17.0 / 64.0, 21.0 / 64.0, 26.0 / 64.0, 30.0 / 64.0,
				 34.0 / 64.0, 38.0 / 64.0, 43.0 / 64.0, 47.0 / 64.0, 51.0 / 64.0, 55.0 / 64.0, 60.0 / 64.0, 1.0} // 4 bit index
#else
		// Pure linear weights
		{0.0},  // 0 bit index
		{0.0, 1.0}, // 1 bit index
		{0.0, 1.0/3.0, 2.0/3.0, 1.0}, // 2 bit index
		{0.0, 1.0/7.0, 2.0/7.0, 3.0/7.0, 4.0/7.0, 5.0/7.0, 6.0/7.0, 1.0}, // 3 bit index
		{0.0, 1.0/15.0, 2.0/15.0, 3.0/15.0, 4.0/15.0, 5.0/15.0, 6.0/15.0, 7.0/15.0,
		 8.0/15.0, 9.0/15.0, 10.0/15.0, 11.0/15.0, 12.0/15.0, 13.0/15.0, 14.0/15.0, 1.0} // 4 bit index
#endif
		};


//
// This routine generates the ramp colours with correct rounding
//
//
//
// Used by BC7_Decode

#ifndef USE_HIGH_PRECISION_INTERPOLATION_BC7
uint16_t aWeight2[] = { 0, 21, 43, 64 };
uint16_t aWeight3[] = { 0, 9, 18, 27, 37, 46, 55, 64 };
uint16_t aWeight4[] = { 0, 4, 9, 13, 17, 21, 26, 30, 34, 38, 43, 47, 51, 55, 60, 64 };

uint8_t interpolate(uint8_t e0, uint8_t e1, uint8_t index, uint8_t indexprecision)
{
		if (indexprecision == 2)
				return (uint8_t)(((64 - aWeight2[index])*uint16_t(e0) + aWeight2[index] * uint16_t(e1) + 32) >> 6);
		else if (indexprecision == 3)
				return (uint8_t)(((64 - aWeight3[index])*uint16_t(e0) + aWeight3[index] * uint16_t(e1) + 32) >> 6);
		else // indexprecision == 4
				return (uint8_t)(((64 - aWeight4[index])*uint16_t(e0) + aWeight4[index] * uint16_t(e1) + 32) >> 6);
}
#endif

void GetRamp(uint32_t endpoint[][MAX_DIMENSION_BIG],
						 float ramp[MAX_DIMENSION_BIG][(1 << MAX_INDEX_BITS)],
						 uint32_t clusters[2],
						 uint32_t componentBits[MAX_DIMENSION_BIG]) {
	float ep[2][MAX_DIMENSION_BIG];
	uint32_t i;

	// Expand each endpoint component to 8 bits by shifting the MSB to bit 7
	// and then replicating the high bits to the low bits revealed by
	// the shift
	for (i = 0; i < MAX_DIMENSION_BIG; i++) {
		ep[0][i] = 0.;
		ep[1][i] = 0.;
		if (componentBits[i]) {
			ep[0][i] = (float) (endpoint[0][i] << (8 - componentBits[i]));
			ep[1][i] = (float) (endpoint[1][i] << (8 - componentBits[i]));
			ep[0][i] += (float) ((uint32_t) ep[0][i] >> componentBits[i]);
			ep[1][i] += (float) ((uint32_t) ep[1][i] >> componentBits[i]);

			ep[0][i] = Math_MinF(255.0f, Math_MaxF(0.0f, ep[0][i]));
			ep[1][i] = Math_MinF(255.0f, Math_MaxF(0.0f, ep[1][i]));
		}
	}

	// If this block type has no explicit alpha channel
	// then make sure alpha is 1.0 for all points on the ramp
	if (!componentBits[COMP_ALPHA]) {
		ep[0][COMP_ALPHA] = ep[1][COMP_ALPHA] = 255.0f;
	}

	uint32_t rampIndex = clusters[0];

	rampIndex = (uint32_t) (log((float) rampIndex) / log(2.0));

	// Generate colours for the RGB ramp
	for (i = 0; i < clusters[0]; i++) {
#ifdef USE_HIGH_PRECISION_INTERPOLATION_BC7
		ramp[COMP_RED][i] = floorf((ep[0][COMP_RED] * (1.0f - rampLerpWeights[rampIndex][i])) +
				(ep[1][COMP_RED] * rampLerpWeights[rampIndex][i]) + 0.5f);
		ramp[COMP_RED][i] = Math_MinF(255.0f, Math_MaxF(0.0f, ramp[COMP_RED][i]));
		ramp[COMP_GREEN][i] = floorf((ep[0][COMP_GREEN] * (1.0f - rampLerpWeights[rampIndex][i])) +
				(ep[1][COMP_GREEN] * rampLerpWeights[rampIndex][i]) + 0.5f);
		ramp[COMP_GREEN][i] = Math_MinF(255.0f, Math_MaxF(0.0f, ramp[COMP_GREEN][i]));
		ramp[COMP_BLUE][i] = floorf((ep[0][COMP_BLUE] * (1.0f - rampLerpWeights[rampIndex][i])) +
				(ep[1][COMP_BLUE] * rampLerpWeights[rampIndex][i]) + 0.5f);
		ramp[COMP_BLUE][i] = Math_MinF(255.0f, Math_MaxF(0.0f, ramp[COMP_BLUE][i]));
#else
		ramp[COMP_RED][i] = interpolate(ep[0][COMP_RED], ep[1][COMP_RED], i, rampIndex);
				ramp[COMP_GREEN][i] = interpolate(ep[0][COMP_GREEN], ep[1][COMP_GREEN], i, rampIndex);
				ramp[COMP_BLUE][i] = interpolate(ep[0][COMP_BLUE], ep[1][COMP_BLUE], i, rampIndex);
#endif
	}

	rampIndex = clusters[1];
	rampIndex = (uint32_t) (log((float) rampIndex) / log(2.0));

	if (!componentBits[COMP_ALPHA]) {
		for (i = 0; i < clusters[1]; i++) {
			ramp[COMP_ALPHA][i] = 255.;
		}
	} else {

		// Generate alphas
		for (i = 0; i < clusters[1]; i++) {
#ifdef USE_HIGH_PRECISION_INTERPOLATION_BC7
			ramp[COMP_ALPHA][i] = floorf((ep[0][COMP_ALPHA] * (1.0f - rampLerpWeights[rampIndex][i])) +
					(ep[1][COMP_ALPHA] * rampLerpWeights[rampIndex][i]) + 0.5f);
			ramp[COMP_ALPHA][i] = Math_MinF(255.0f, Math_MaxF(0.f, ramp[COMP_ALPHA][i]));
#else
			ramp[COMP_ALPHA][i] = interpolate(ep[0][COMP_ALPHA], ep[1][COMP_ALPHA], i, rampIndex);
#endif
		}

	}
}
// Threshold quality below which we will always run fast quality and shaking
// Selfnote: User should be able to set this?
// Default FQuality is at 0.1 < g_qFAST_THRESHOLD which will cause the SingleIndex compression to start skipping shape blocks
// during compression
// if user sets a value above this then all shapes will be used for compression scan for quality
float g_qFAST_THRESHOLD = 0.5f;

// This limit is used for DualIndex Block and if fQuality is above this limit then Quantization shaking will always be performed
// on all indexs
float  g_HIGHQULITY_THRESHOLD = 0.7f;

//
// For a given block mode this sets up the data needed by the compressor
//
// Note that BC7 only uses NO_PBIT, ONE_PBIT and TWO_PBIT encodings
// for endpoints
//
void BC7BlockEncoder::BlockSetup(uint32_t blockMode) {
	switch (bti[blockMode].pBitType) {
	case NO_PBIT: m_parityBits = CART;
		break;
	case ONE_PBIT: m_parityBits = SAME_PAR;
		break;
	case TWO_PBIT: m_parityBits = BCC;
		break;
	case THREE_PBIT: m_parityBits = SAME_FCC;
		break;
	case FOUR_PBIT: m_parityBits = FCC;
		break;
	case FIVE_PBIT: m_parityBits = FCC_SAME_BCC;
		break;
	}

	if (bti[blockMode].encodingType == NO_ALPHA) {
		m_componentBits[COMP_RED] = bti[blockMode].vectorBits / 3;
		m_componentBits[COMP_GREEN] = bti[blockMode].vectorBits / 3;
		m_componentBits[COMP_BLUE] = bti[blockMode].vectorBits / 3;
		m_componentBits[COMP_ALPHA] = 0;

		m_clusters[0] = 1 << bti[blockMode].indexBits[0];
		m_clusters[1] = 0;
	} else if (bti[blockMode].encodingType == COMBINED_ALPHA) {
		m_componentBits[COMP_RED] = bti[blockMode].vectorBits / 4;
		m_componentBits[COMP_GREEN] = bti[blockMode].vectorBits / 4;
		m_componentBits[COMP_BLUE] = bti[blockMode].vectorBits / 4;
		m_componentBits[COMP_ALPHA] = bti[blockMode].vectorBits / 4;

		m_clusters[0] = 1 << bti[blockMode].indexBits[0];
		m_clusters[1] = 0;
	} else if (bti[blockMode].encodingType == SEPARATE_ALPHA) {
		m_componentBits[COMP_RED] = bti[blockMode].vectorBits / 3;
		m_componentBits[COMP_GREEN] = bti[blockMode].vectorBits / 3;
		m_componentBits[COMP_BLUE] = bti[blockMode].vectorBits / 3;
		m_componentBits[COMP_ALPHA] = bti[blockMode].scalarBits;

		m_clusters[0] = 1 << bti[blockMode].indexBits[0];
		m_clusters[1] = 1 << bti[blockMode].indexBits[1];
	}
}

//
// This function sorts out the bit encoding for the BC7 block and packs everything
// in the right order for the hardware decoder
//
//
//

void BC7BlockEncoder::EncodeSingleIndexBlock(uint32_t blockMode,
																						 uint32_t partition,
																						 uint32_t colour[MAX_SUBSETS][2],
																						 int indices[MAX_SUBSETS][MAX_SUBSET_SIZE],
																						 uint8_t block[COMPRESSED_BLOCK_SIZE]) {
	uint32_t i, j, k;
	uint32_t *partitionTable;
	int bitPosition = 0;    // Position the pointer at the LSB
	uint8_t *basePtr = (uint8_t *) block;
	uint32_t blockIndices[MAX_SUBSET_SIZE];

	// Generate Unary header
	for (i = 0; i < (int) blockMode; i++) {
		WriteBit(basePtr, bitPosition++, 0);
	}
	WriteBit(basePtr, bitPosition++, 1);

	// Write partition bits
	for (i = 0; i < bti[blockMode].partitionBits; i++) {
		WriteBit(basePtr, bitPosition++, (uint8_t) (partition >> i) & 0x1);
	}

	// Extract the index bits from the partitions
	partitionTable = (uint32_t *) BC7_PARTITIONS[bti[blockMode].subsetCount - 1][partition];

	uint32_t idxCount[3] = {0, 0, 0};
	bool flipColours[3] = {false, false, false};

	// Sort out the index set and tag whether we need to flip the
	// endpoints to get the correct state in the implicit index bits
	// The implicitly encoded MSB of the fixup index must be 0
	uint32_t fixup[3] = {0, 0, 0};
	switch (bti[blockMode].subsetCount) {
	case 3: fixup[1] = BC7_FIXUPINDICES[2][partition][1];
		fixup[2] = BC7_FIXUPINDICES[2][partition][2];
		break;
	case 2: fixup[1] = BC7_FIXUPINDICES[1][partition][1];
		break;
	default: break;
	}

	// Extract indices and mark subsets that need to have their colours flipped to get the
	// right state for the implicit MSB of the fixup index
	for (i = 0; i < MAX_SUBSET_SIZE; i++) {
		uint32_t p = partitionTable[i];
		blockIndices[i] = indices[p][idxCount[p]++];

		for (j = 0; j < (int) bti[blockMode].subsetCount; j++) {
			if (i == fixup[j]) {
				if (blockIndices[i] & (1 << (bti[blockMode].indexBits[0] - 1))) {
					flipColours[j] = true;
				}
			}
		}
	}

	// Now we must flip the endpoints where necessary so that the implicitly encoded
	// index bits have the correct state
	for (i = 0; i < (int) bti[blockMode].subsetCount; i++) {
		if (flipColours[i]) {
			uint32_t temp;
			temp = colour[i][0];
			colour[i][0] = colour[i][1];
			colour[i][1] = temp;
		}
	}

	// ...next flip the indices where necessary
	for (i = 0; i < MAX_SUBSET_SIZE; i++) {
		uint32_t p = partitionTable[i];
		if (flipColours[p]) {
			blockIndices[i] = ((1 << bti[blockMode].indexBits[0]) - 1) - blockIndices[i];
		}
	}

	uint32_t subset, ep, component;

	// Endpoints are stored in the following order RRRR GGGG BBBB (AAAA) (PPPP)
	// i.e. components are packed together
	uint32_t unpackedColours[MAX_SUBSETS][2][MAX_DIMENSION_BIG];
	uint32_t parityBits[MAX_SUBSETS][2];

	// Unpack the colour values for the subsets
	for (i = 0; i < bti[blockMode].subsetCount; i++) {
		uint32_t packedColours[2] = {colour[i][0],
																 colour[i][1]};

		if (bti[blockMode].pBitType == TWO_PBIT) {
			parityBits[i][0] = packedColours[0] & 1;
			parityBits[i][1] = packedColours[1] & 1;
			packedColours[0] >>= 1;
			packedColours[1] >>= 1;
		} else if (bti[blockMode].pBitType == ONE_PBIT) {
			parityBits[i][0] = packedColours[1] & 1;
			parityBits[i][1] = packedColours[1] & 1;
			packedColours[0] >>= 1;
			packedColours[1] >>= 1;
		} else {
			parityBits[i][0] = 0;
			parityBits[i][1] = 0;
		}

		uint32_t component1;
		for (component1 = 0; component1 < MAX_DIMENSION_BIG; component1++) {
			if (m_componentBits[component1]) {
				unpackedColours[i][0][component1] = packedColours[0] & ((1 << m_componentBits[component1]) - 1);
				unpackedColours[i][1][component1] = packedColours[1] & ((1 << m_componentBits[component1]) - 1);
				packedColours[0] >>= m_componentBits[component1];
				packedColours[1] >>= m_componentBits[component1];
			}
		}
	}

	// Loop over components
	for (component = 0; component < MAX_DIMENSION_BIG; component++) {
		// loop over subsets
		for (subset = 0; subset < (int) bti[blockMode].subsetCount; subset++) {
			// Loop over endpoints and write colour bits
			for (ep = 0; ep < 2; ep++) {
				// Write this component
				for (k = 0; k < m_componentBits[component]; k++) {
					WriteBit(basePtr,
									 bitPosition++,
									 (uint8_t) (unpackedColours[subset][ep][component] >> k) & 0x1);
				}
			}
		}
	}

	// Now write parity bits if present
	if (bti[blockMode].pBitType != NO_PBIT) {
		for (subset = 0; subset < (int) bti[blockMode].subsetCount; subset++) {
			if (bti[blockMode].pBitType == ONE_PBIT) {
				WriteBit(basePtr,
								 bitPosition++,
								 parityBits[subset][0] & 1);
			} else if (bti[blockMode].pBitType == TWO_PBIT) {
				WriteBit(basePtr,
								 bitPosition++,
								 parityBits[subset][0] & 1);
				WriteBit(basePtr,
								 bitPosition++,
								 parityBits[subset][1] & 1);
			}
		}
	}

	// Now encode the index bits
	for (i = 0; i < MAX_SUBSET_SIZE; i++) {
		uint32_t p = partitionTable[i];
		// If this is a fixup index then drop the MSB which is implicitly 0
		if (i == fixup[p]) {
			for (j = 0; j < (bti[blockMode].indexBits[0] - 1); j++) {
				WriteBit(basePtr, bitPosition++, (uint8_t) (blockIndices[i] >> j));
			}
		} else {
			for (j = 0; j < bti[blockMode].indexBits[0]; j++) {
				WriteBit(basePtr, bitPosition++, (uint8_t) (blockIndices[i] >> j));
			}
		}
	}

	// Check that we encoded exactly the right number of bits
	if (bitPosition != (COMPRESSED_BLOCK_SIZE * 8)) {
		return;
	}
}

//
// This routine can be used to compress a block to any of the modes with a shared index set
//
// It will encode the best result for this mode into a BC7 block
//
//
//
float BC7BlockEncoder::CompressSingleIndexBlock(float const in[MAX_SUBSET_SIZE * MAX_DIMENSION_BIG],
																								 uint8_t out[COMPRESSED_BLOCK_SIZE],
																								 uint32_t blockMode) {
	uint32_t i, k, n;
	uint32_t dimension;

	// Figure out the effective dimension of this block mode
	if (bti[blockMode].encodingType == NO_ALPHA) {
		dimension = 3;
	} else {
		dimension = 4;
	}

	uint32_t numPartitionModes = 1 << bti[blockMode].partitionBits;
	uint32_t partitionsToTry = numPartitionModes;

	// Linearly reduce the number of partitions to try as the quality falls below a threshold
	if (m_quality < g_qFAST_THRESHOLD) {
		partitionsToTry = (uint32_t) floor((float) (partitionsToTry * m_partitionSearchSize) + 0.5);
		partitionsToTry = Math_MinU32(numPartitionModes, Math_MaxU32(1, partitionsToTry));
	}

	uint32_t blockPartition;
	float partition[MAX_SUBSETS][MAX_SUBSET_SIZE][MAX_DIMENSION_BIG];
	uint32_t entryCount[MAX_SUBSETS];
	uint32_t subset;

	// Loop over the available partitions for the block mode and quantize them
	// to figure out the best candidates for further refinement
	for (blockPartition = 0;
			 blockPartition < partitionsToTry;
			 blockPartition++) {
		Partition(blockPartition,
							in,
							partition,
							entryCount,
							blockMode,
							dimension);

		float error = 0.;
		float outB[MAX_SUBSET_SIZE][MAX_DIMENSION_BIG];
		float direction[MAX_DIMENSION_BIG];
		float step;

		for (subset = 0; subset < bti[blockMode].subsetCount; subset++) {
			int indices[MAX_SUBSETS][MAX_SUBSET_SIZE];

			if (entryCount[subset]) {

				if ((m_clusters[0] > 8) ||
						(m_blockMaxRange <= m_quantizerRangeThreshold)) {
					error += optQuantAnD_d(partition[subset],
																 entryCount[subset],
																 m_clusters[0],
																 indices[subset],
																 outB,
																 direction,
																 &step,
																 dimension);

				} else {
					error += optQuantTrace_d(partition[subset],
																	 entryCount[subset],
																	 m_clusters[0],
																	 indices[subset],
																	 outB,
																	 direction,
																	 &step,
																	 dimension);

				}

				// Store off the indices for later
				for (uint32_t idx = 0; idx < entryCount[subset]; idx++) {
					m_storedIndices[blockPartition][subset][idx] = indices[subset][idx];
				}
			}
		}

		m_storedError[blockPartition] = error;
	}

	// Sort the results
	sortProjection(m_storedError,
								 m_sortedModes,
								 partitionsToTry);


	// Run shaking (endpoint refinement) pass for partitions that gave the
	// best set of errors from quantization

	// ep_shaker will take its endpoint information from bits[0-2]
	// ep_shaker_2_d will take its information from bits[3]
	int bits[4] = {0, 0, 0, 0};

	// ep_shaker_d needs bits specified individually per channel including parity
	bits[0] = m_componentBits[COMP_RED] + (m_parityBits ? 1 : 0);
	bits[1] = m_componentBits[COMP_GREEN] + (m_parityBits ? 1 : 0);
	bits[2] = m_componentBits[COMP_BLUE] + (m_parityBits ? 1 : 0);

	// ep_shaker_2_d needs bits specified as total bits for both endpoints including parity
	for (i = 0; i < dimension; i++) {
		bits[3] += m_componentBits[i];
	}
	bits[3] *= 2;
	if (m_parityBits == BCC) {
		bits[3] += 2;
	} else if (m_parityBits == SAME_PAR) {
		bits[3] += 1;
	}

	int epo_code[MAX_SUBSETS][2][MAX_DIMENSION_BIG];
	float epo[2][MAX_DIMENSION_BIG];
	float outB[MAX_SUBSET_SIZE][MAX_DIMENSION_BIG];

	int bestEndpoints[MAX_SUBSETS][2][MAX_DIMENSION_BIG];
	int bestIndices[MAX_SUBSETS][MAX_SUBSET_SIZE];
	uint32_t bestEntryCount[MAX_SUBSETS];
	uint32_t bestPartition = 0;
	float bestError = FLT_MAX;

	// Extensive shaking is most important when the ramp is short, and
	// when we have less indices. On a long ramp the quality of the
	// initial quantizing is relatively more important
	// We modulate the shake size according to the number of ramp indices
	// - the more indices we have the less shaking should be required to find a near
	// optimal match

	// shakeSize gives the size of the shake cube (for ep_shaker_2_d)
	// ep_shaker always runs on a 1x1x1 cube on both endpoints
	uint32_t shakeSize = 8 - (uint32_t) floorf(1.5f * bti[blockMode].indexBits[0]);
	shakeSize = Math_MaxU32(2, Math_MinU32((uint32_t) floorf(shakeSize * (float) m_quality + 0.5f), 6));

	// Shake attempts indicates how many partitions to try to shake
	uint32_t numShakeAttempts =
			Math_MaxU32(1, Math_MinU32((uint32_t) floorf(8.0f * (float) m_quality + 0.5f), partitionsToTry));

	// Set up all the parameters for the shakers
	// Must increase shake size if these block endpoints use parity
	if ((m_parityBits == SAME_PAR) ||
			(m_parityBits == BCC)) {
		shakeSize += 2;
	}

	// Now do the endpoint shaking
	for (i = 0; i < numShakeAttempts; i++) {
		float error = 0;

		blockPartition = m_sortedModes[i];

		Partition(blockPartition,
							in,
							partition,
							entryCount,
							blockMode,
							dimension);

		for (subset = 0; subset < bti[blockMode].subsetCount; subset++) {
			if (entryCount[subset]) {
				// If quality is set low or the dimension is not compatible with
				// shaker_d then just run shaker_2_d
				if ((m_blockMaxRange > m_shakerRangeThreshold) ||
						(dimension != 3)) {
					error += ep_shaker_2(partition[subset],
																 entryCount[subset],
																 m_storedIndices[blockPartition][subset],
																 outB,
																 epo_code[subset],
																 shakeSize,
																 m_clusters[0] - 1,
																 bits[3],
																 dimension,
																 epo);
				} else {
					float tempError[2];
					int tempIndices[MAX_SUBSET_SIZE];
					int temp_epo_code[2][MAX_DIMENSION_BIG];

					// Step one - run ep_shaker and ep_shaker_2 in parallel, and get the error from each

					for (k = 0; k < entryCount[subset]; k++) {
						tempIndices[k] = m_storedIndices[blockPartition][subset][k];
					}

					tempError[0] = ep_shaker(partition[subset],
																		 entryCount[subset],
																		 tempIndices,
																		 outB,
																		 temp_epo_code,
																		 m_clusters[0] - 1,
																		 bits,
																		 (CMP_qt) m_parityBits,
																		 dimension);

					tempError[1] = ep_shaker_2(partition[subset],
																			 entryCount[subset],
																			 m_storedIndices[blockPartition][subset],
																			 outB,
																			 epo_code[subset],
																			 shakeSize,
																			 m_clusters[0] - 1,
																			 bits[3],
																			 dimension,
																			 epo);

					if (tempError[0] < tempError[1]) {
						// If ep_shaker did better than ep_shaker_2 then we need to reshake
						// the output from ep_shaker using ep_shaker_2 for further refinement

						tempError[1] = ep_shaker_2(partition[subset],
																				 entryCount[subset],
																				 tempIndices,
																				 outB,
																				 temp_epo_code,
																				 shakeSize,
																				 m_clusters[0] - 1,
																				 bits[3],
																				 dimension,
																				 epo);

						// Copy the results into the expected location
						for (k = 0; k < entryCount[subset]; k++) {
							m_storedIndices[blockPartition][subset][k] = tempIndices[k];
						}

						for (k = 0; k < MAX_DIMENSION_BIG; k++) {
							epo_code[subset][0][k] = temp_epo_code[0][k];
							epo_code[subset][1][k] = temp_epo_code[1][k];
						}
					}

					error += tempError[1];
				}
			}
		}

		if (error < bestError) {
			bestPartition = blockPartition;

			for (subset = 0; subset < bti[blockMode].subsetCount; subset++) {
				bestEntryCount[subset] = entryCount[subset];

				if (entryCount[subset]) {
					for (k = 0; k < dimension; k++) {
						bestEndpoints[subset][0][k] = epo_code[subset][0][k];
						bestEndpoints[subset][1][k] = epo_code[subset][1][k];
					}

					for (n = 0; n < entryCount[subset]; n++) {
						bestIndices[subset][n] = m_storedIndices[blockPartition][subset][n];
					}
				}
			}

			bestError = error;
		}

		// Early out if we  found we can compress with error below the quality threshold
		if (m_errorThreshold > 0) {
			if (bestError <= m_errorThreshold) {
				break;
			}
		}
	}

	// Now we have all the data needed to encode the block
	// We need to pack the endpoints prior to encoding
	uint32_t packedEndpoints[3][2];
	for (subset = 0; subset < bti[blockMode].subsetCount; subset++) {
		if (bestEntryCount[subset]) {
			uint32_t rightAlignment = 0;
			packedEndpoints[subset][0] = 0;
			packedEndpoints[subset][1] = 0;

			// Sort out parity bits
			if (m_parityBits != CART) {
				packedEndpoints[subset][0] = bestEndpoints[subset][0][0] & 1;
				packedEndpoints[subset][1] = bestEndpoints[subset][1][0] & 1;
				for (k = 0; k < MAX_DIMENSION_BIG; k++) {
					bestEndpoints[subset][0][k] >>= 1;
					bestEndpoints[subset][1][k] >>= 1;
				}
				rightAlignment++;
			}

			// Fixup endpoints
			for (k = 0; k < dimension; k++) {
				if (m_componentBits[k]) {
					packedEndpoints[subset][0] |= bestEndpoints[subset][0][k] << rightAlignment;
					packedEndpoints[subset][1] |= bestEndpoints[subset][1][k] << rightAlignment;
					rightAlignment += m_componentBits[k];
				}
			}
		}
	}

	// Save the data to output
	EncodeSingleIndexBlock(blockMode,
												 bestPartition,
												 packedEndpoints,
												 bestIndices,
												 out);
	return bestError;
}

static uint32_t componentRotations[4][4] =
		{
				{COMP_ALPHA, COMP_RED, COMP_GREEN, COMP_BLUE},
				{COMP_RED, COMP_ALPHA, COMP_GREEN, COMP_BLUE},
				{COMP_GREEN, COMP_RED, COMP_ALPHA, COMP_BLUE},
				{COMP_BLUE, COMP_RED, COMP_GREEN, COMP_ALPHA}
		};

void BC7BlockEncoder::EncodeDualIndexBlock(uint32_t blockMode,
																					 uint32_t indexSelection,
																					 uint32_t componentRotation,
																					 int endpoint[2][2][MAX_DIMENSION_BIG],
																					 int indices[2][MAX_SUBSET_SIZE],
																					 uint8_t out[COMPRESSED_BLOCK_SIZE]) {
#ifdef USE_DBGTRACE
	DbgTrace(());
#endif
	uint32_t i, j, k;
	int bitPosition = 0;    // Position the pointer at the LSB
	uint8_t *basePtr = out;
	uint32_t idxBits[2];
	bool swapIndices;

	// Generate Unary header for this mode
	for (i = 0; i < blockMode; i++) {
		WriteBit(basePtr, bitPosition++, 0);
	}
	WriteBit(basePtr, bitPosition++, 1);

	// Write rotation bits
	for (i = 0; i < bti[blockMode].rotationBits; i++) {
		WriteBit(basePtr, bitPosition++, (uint8_t) ((componentRotation >> i) & 0xff));
	}

	// Write index selector bits
	for (i = 0; i < bti[blockMode].indexModeBits; i++) {
		WriteBit(basePtr, bitPosition++, (uint8_t) (indexSelection ? 1 : 0));
	}

	if (indexSelection) {
		swapIndices = TRUE;
		idxBits[0] = bti[blockMode].indexBits[1];
		idxBits[1] = bti[blockMode].indexBits[0];
	} else {
		swapIndices = FALSE;
		idxBits[0] = bti[blockMode].indexBits[0];
		idxBits[1] = bti[blockMode].indexBits[1];
	}

	bool flipColours[2] = {false, false};

	// Indicate if we need to fixup the indices
	if (indices[0][0] & (1 << (idxBits[0] - 1))) {
		flipColours[0] = true;
	}
	if (indices[1][0] & (1 << (idxBits[1] - 1))) {
		flipColours[1] = true;
	}

	// Fixup the indices
	for (i = 0; i < 2; i++) {
		if (flipColours[i]) {
			for (j = 0; j < MAX_SUBSET_SIZE; j++) {
				indices[i][j] = ((1 << idxBits[i]) - 1) - indices[i][j];
			}
		}
	}

	// Now fixup the endpoints so that the implicitly encoded
	// index bits have the correct state
	for (i = 0; i < 2; i++) {
		if (flipColours[i]) {
			for (k = 0; k < 4; k++) {
				uint32_t temp;
				temp = endpoint[i][0][k];
				endpoint[i][0][k] = endpoint[i][1][k];
				endpoint[i][1][k] = temp;
			}
		}
	}

	uint32_t ep, component;
	// Encode the colour and alpha information
	uint32_t vectorComponentBits = bti[blockMode].vectorBits / 3;

	// Loop over components
	for (component = 0; component < MAX_DIMENSION_BIG; component++) {
		if (component != COMP_ALPHA) {
			for (ep = 0; ep < 2; ep++) {
				for (k = 0; k < vectorComponentBits; k++) {
					WriteBit(basePtr,
									 bitPosition++,
									 (uint8_t) ((endpoint[0][ep][component] >> k) & 0x1));
				}
			}
		} else {
			for (ep = 0; ep < 2; ep++) {
				for (j = 0; j < bti[blockMode].scalarBits; j++) {
					WriteBit(basePtr,
									 bitPosition++,
									 (uint8_t) ((endpoint[1][ep][0] >> j) & 0x1));
				}
			}
		}
	}

	// Now encode the index bits
	for (i = 0; i < 2; i++) {
		uint32_t idxSelect = i;

		if (swapIndices) {
			idxSelect = i ^ 1;
		}
		for (j = 0; j < MAX_SUBSET_SIZE; j++) {
			if (j == 0) {
				for (k = 0; k < (idxBits[idxSelect] - 1); k++) {
					WriteBit(basePtr, bitPosition++, (uint8_t) (indices[idxSelect][j] >> k));
				}
			} else {
				for (k = 0; k < idxBits[idxSelect]; k++) {
					WriteBit(basePtr, bitPosition++, (uint8_t) (indices[idxSelect][j] >> k));
				}
			}
		}
	}

	// Check that we encoded exactly the right number of bits
	if (bitPosition != (COMPRESSED_BLOCK_SIZE * 8)) {
		return;
	}
}

float BC7BlockEncoder::CompressDualIndexBlock(float const in[MAX_SUBSET_SIZE * MAX_DIMENSION_BIG],
																							 uint8_t out[COMPRESSED_BLOCK_SIZE],
																							 uint32_t blockMode) {
	uint32_t i;
	float cBlock[MAX_SUBSET_SIZE][MAX_DIMENSION_BIG];
	float aBlock[MAX_SUBSET_SIZE][MAX_DIMENSION_BIG];

	uint32_t maxRotation = 1 << bti[blockMode].rotationBits;
	uint32_t rotation;

	uint32_t maxIndexSelection = 1 << bti[blockMode].indexModeBits;
	uint32_t indexSelection;

	int indices[2][MAX_SUBSET_SIZE];
	float outQ[2][MAX_SUBSET_SIZE][MAX_DIMENSION_BIG];
	float direction[MAX_DIMENSION_BIG];
	float step;

	float quantizerError;
	float bestQuantizerError = FLT_MAX;
	float overallError;
	float bestOverallError = FLT_MAX;

	// Go through each possible rotation and selection of indices
	for (rotation = 0; rotation < maxRotation; rotation++) { // A


		for (i = 0; i < MAX_SUBSET_SIZE; i++) {
			cBlock[i][COMP_RED] = in[i * 4  + componentRotations[rotation][1]];
			cBlock[i][COMP_GREEN] = in[i * 4  + componentRotations[rotation][2]];
			cBlock[i][COMP_BLUE] = in[i * 4  + componentRotations[rotation][3]];

			aBlock[i][COMP_RED] = in[i * 4  + componentRotations[rotation][0]];
			aBlock[i][COMP_GREEN] = in[i * 4  + componentRotations[rotation][0]];
			aBlock[i][COMP_BLUE] = in[i * 4  + componentRotations[rotation][0]];
		}

		for (indexSelection = 0; indexSelection < maxIndexSelection; indexSelection++) { // B
			quantizerError = 0.;
			// Quantize the vector block
			if (m_blockMaxRange <= m_quantizerRangeThreshold) {

				quantizerError = optQuantAnD_d(cBlock,
																			 MAX_SUBSET_SIZE,
																			 (1 << bti[blockMode].indexBits[0 ^ indexSelection]),
																			 indices[0],
																			 outQ[0],
																			 direction,
																			 &step,
																			 3);

			} else {
				quantizerError = optQuantTrace_d(cBlock,
																				 MAX_SUBSET_SIZE,
																				 (1 << bti[blockMode].indexBits[0 ^ indexSelection]),
																				 indices[0],
																				 outQ[0],
																				 direction,
																				 &step,
																				 3);

			}

			// Quantize the scalar block
			if (m_blockMaxRange <= m_quantizerRangeThreshold) {
				quantizerError += optQuantAnD_d(aBlock,
																				MAX_SUBSET_SIZE,
																				(1 << bti[blockMode].indexBits[1 ^ indexSelection]),
																				indices[1],
																				outQ[1],
																				direction,
																				&step,
																				3) / 3.0f;
			} else {
				quantizerError += optQuantTrace_d(aBlock,
																					MAX_SUBSET_SIZE,
																					(1 << bti[blockMode].indexBits[1 ^ indexSelection]),
																					indices[1],
																					outQ[1],
																					direction,
																					&step,
																					3) / 3.0f;

			}

			// If quality is high then run the full shaking for this config and
			// store the result if it beats the best overall error
			// Otherwise only run the shaking if the error is better than the best
			// quantizer error
			if ((m_quality > g_HIGHQULITY_THRESHOLD) || (quantizerError <= bestQuantizerError)) {
				// Shake size gives the size of the shake cube
				uint32_t shakeSize;

				shakeSize = Math_MaxU32(2, Math_MinU32((uint32_t) (6 * m_quality), 6));

				int bits[2][4];

				// Specify number of bits for vector block
				bits[0][COMP_RED] = m_componentBits[COMP_RED];
				bits[0][COMP_GREEN] = m_componentBits[COMP_GREEN];
				bits[0][COMP_BLUE] = m_componentBits[COMP_BLUE];
				bits[0][3] = 2 * (m_componentBits[COMP_RED] + m_componentBits[COMP_GREEN] + m_componentBits[COMP_BLUE]);

				// Specify number of bits for scalar block
				bits[1][0] = m_componentBits[COMP_ALPHA];
				bits[1][1] = m_componentBits[COMP_ALPHA];
				bits[1][2] = m_componentBits[COMP_ALPHA];
				bits[1][3] = 6 * m_componentBits[COMP_ALPHA];

				overallError = 0;
				int epo_code[2][2][MAX_DIMENSION_BIG];
				float epo[2][MAX_DIMENSION_BIG];

				if (m_blockMaxRange > m_shakerRangeThreshold) {
					overallError += ep_shaker_2(cBlock,
																				MAX_SUBSET_SIZE,
																				indices[0],
																				outQ[0],
																				epo_code[0],
																				shakeSize,
																				(1 << bti[blockMode].indexBits[0 ^ indexSelection]) - 1,
																				bits[0][3],
																				3,
																				epo);
				} else {
					ep_shaker(cBlock,
											MAX_SUBSET_SIZE,
											indices[0],
											outQ[0],
											epo_code[0],
											(1 << bti[blockMode].indexBits[0 ^ indexSelection]) - 1,
											bits[0],
											(CMP_qt) 0,
											3);

					overallError += ep_shaker_2(cBlock,
																				MAX_SUBSET_SIZE,
																				indices[0],
																				outQ[0],
																				epo_code[0],
																				shakeSize,
																				(1 << bti[blockMode].indexBits[0 ^ indexSelection]) - 1,
																				bits[0][3],
																				3,
																				epo);
				}

				if (m_blockMaxRange > m_shakerRangeThreshold) {
					overallError += ep_shaker_2(aBlock,
																				MAX_SUBSET_SIZE,
																				indices[1],
																				outQ[1],
																				epo_code[1],
																				shakeSize,
																				(1 << bti[blockMode].indexBits[1 ^ indexSelection]) - 1,
																				bits[1][3],
																				3,
																				epo) / 3.0f;
				} else {
					ep_shaker(aBlock,
											MAX_SUBSET_SIZE,
											indices[1],
											outQ[1],
											epo_code[1],
											(1 << bti[blockMode].indexBits[1 ^ indexSelection]) - 1,
											bits[1],
											(CMP_qt) 0,
											3);

					overallError += ep_shaker_2(aBlock,
																				MAX_SUBSET_SIZE,
																				indices[1],
																				outQ[1],
																				epo_code[1],
																				shakeSize,
																				(1 << bti[blockMode].indexBits[1 ^ indexSelection]) - 1,
																				bits[1][3],
																				3,
																				epo) / 3.0f;
				}

				// If we beat the previous best then encode the block
				if (overallError < bestOverallError) {
					EncodeDualIndexBlock(blockMode,
															 indexSelection,
															 rotation,
															 epo_code,
															 indices,
															 out);

					bestOverallError = overallError;
				}

				if (quantizerError < bestQuantizerError) {
					bestQuantizerError = quantizerError;
				}

			}
		} // B
	} // A
	return bestOverallError;
}



//
// This routine compresses a block and returns the RMS error
//
//
//
//

float BC7BlockEncoder::CompressBlock(float const inN[MAX_SUBSET_SIZE * MAX_DIMENSION_BIG],
																			uint8_t out[COMPRESSED_BLOCK_SIZE]) {
	uint32_t i, j;
	bool blockNeedsAlpha = FALSE;
	bool blockAlphaZeroOne = FALSE;
	uint32_t validModeMask = m_validModeMask;
	bool encodedBlock = FALSE;

	for (i = 0; i < MAX_DIMENSION_BIG; i++) {
		m_blockMin[i] = FLT_MAX;
		m_blockMax[i] = 0.0;
		m_blockRange[i] = 0.0;
	}
	
	float in[MAX_SUBSET_SIZE * MAX_DIMENSION_BIG];

	// Check if the input block has any alpha values that are not 1
	// We assume 8-bit input here, so 1 is mapped to 255.
	// Also check if the block encodes an explicit zero or one in the
	// alpha channel. If so then we might need also need special as the
	// block may have a thresholded or punch-through alpha

	for (i = 0; i < MAX_SUBSET_SIZE; i++) {
		if (inN[i * 4  + COMP_ALPHA] < 1.0) {
			blockNeedsAlpha = TRUE;
		} else if ((in[i * 4  + COMP_ALPHA] == 1.0) ||
				(inN[i * 4  + COMP_ALPHA] == 0.0)) {
			blockAlphaZeroOne = TRUE;
		}
		for (j = 0; j < MAX_DIMENSION_BIG; j++) {
			in[i * 4 + j] = inN[i * 4 + j] * 255.0f;
			m_blockMin[j] = (in[i * 4  + j] < m_blockMin[j]) ? in[i * 4  + j] : m_blockMin[j];
			m_blockMax[j] = (in[i * 4  + j] > m_blockMax[j]) ? in[i * 4  + j] : m_blockMax[j];
		}

	}

	m_blockRange[0] = m_blockMax[0] - m_blockMin[0];
	m_blockRange[1] = m_blockMax[1] - m_blockMin[1];
	m_blockRange[2] = m_blockMax[2] - m_blockMin[2];
	m_blockRange[3] = m_blockMax[3] - m_blockMin[3];
	m_blockMaxRange = Math_MaxF(m_blockRange[0], m_blockRange[1]);
	m_blockMaxRange = Math_MaxF(m_blockMaxRange, m_blockRange[2]);
	m_blockMaxRange = Math_MaxF(m_blockMaxRange, m_blockRange[3]);

	bool solidColour = false;
	if (m_blockMaxRange < 1e-10)
		solidColour = true;

	// Initial loop - go through the block modes and get the ones that are valid
	for (uint32_t blockMode = 0; blockMode < NUM_BLOCK_TYPES; blockMode++) {
		// Check if this mode is allowed based on the global settings
		if (!(validModeMask & (1 << blockMode))) {
			continue;
		}

		// If the block needs Alpha and this mode doesn't support alpha then
		// indicate that this is not a valid mode and continue
		if ((blockNeedsAlpha == TRUE) &&
				(bti[blockMode].encodingType == NO_ALPHA)) {
			validModeMask &= ~(1 << blockMode);
		}

		// Optional restriction for colour-only blocks so that they
		// don't use modes that have combined colour+alpha - this
		// avoids the possibility that the encoder might choose an
		// alpha other than 1.0 (due to parity) and cause something to
		// become accidentally slightly transparent (it's possible that
		// when encoding 3-component texture applications will assume that
		// the 4th component can safely be assumed to be 1.0 all the time)
		if (!solidColour &&
				(blockNeedsAlpha == FALSE) &&
				(m_colourRestrict == TRUE) &&
				(bti[blockMode].encodingType == COMBINED_ALPHA)) {
			validModeMask &= ~(1 << blockMode);
		}

		// Optional restriction for blocks with alpha to avoid issues with
		// punch-through or thresholded alpha encoding
		if ((blockNeedsAlpha == TRUE) &&
				(m_alphaRestrict == TRUE) &&
				(blockAlphaZeroOne == TRUE) &&
				(bti[blockMode].encodingType == COMBINED_ALPHA)) {
			validModeMask &= ~(1 << blockMode);
		}
	}

	ASSERT(validModeMask != 0);

	// Try all the legal block modes that we flagged

	uint8_t temporaryOutputBlock[COMPRESSED_BLOCK_SIZE];
	float bestError = FLT_MAX;
	float thisError;
	uint32_t bestblockMode = 99;

	// We change the order in which we visit the block modes to try to maximize the chance
	// that we manage to early out as quickly as possible.
	// This is a significant performance optimization for the lower quality modes where the
	// exit threshold is higher, and also tends to improve quality (as the generally higher quality
	// modes are now enumerated earlier, so the first encoding that passes the threshold will
	// tend to pass by a greater margin than if we used a dumb ordering, and thus overall error will
	// be improved)
	// Deano pushed 6 to the front, as solid colour block have error 0 and so are pick the first
	// block. 6 is the best for pure RGB colour blocks
	// mode 0 is technically illegal so only try it last of all (6 should do its job well)
	uint32_t blockModeOrder[NUM_BLOCK_TYPES] = {6, 4, 3, 1, 2, 0, 7, 5};

	for (uint32_t j1 = 0; j1 < NUM_BLOCK_TYPES; j1++) {
		uint32_t blockMode = blockModeOrder[j1];
		uint32_t Mode = 0x0001 << blockMode;

		if (!(validModeMask & Mode)) {
			continue;
		}

		// Setup mode parameters for this block
		BlockSetup(blockMode);

		if (bti[blockMode].encodingType != SEPARATE_ALPHA) {

			thisError = CompressSingleIndexBlock(in, temporaryOutputBlock, blockMode);
		} else {
			thisError = CompressDualIndexBlock(in, temporaryOutputBlock, blockMode);
		}

		// If this compression did better than all previous attempts then copy the result
		// to the output block
		if (thisError < bestError) {
			for (i = 0; i < COMPRESSED_BLOCK_SIZE; i++) {
				out[i] = temporaryOutputBlock[i];
			}
			bestError = thisError;
			encodedBlock = TRUE;
			bestblockMode = blockMode;
		}

		// If we have achieved an error lower than the requirement threshold then just exit now
		// Early out if we  found we can compress with error below the quality threshold
		if (m_errorThreshold > 0) {
			if (bestError <= m_errorThreshold) {
				break;
			}
		}
	}

	if (bestError < m_smallestError) {
		m_smallestError = bestError;
	}
	if (bestError > m_largestError) {
		m_largestError = bestError;
	}

	if (!encodedBlock) {
		// return some sort of error and abort sequence!
		encodedBlock = FALSE;
	}

	return bestError;
}
