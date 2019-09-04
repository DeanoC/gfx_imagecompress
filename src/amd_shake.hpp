void init_ramps ();


float ep_shaker_2_(
		float data[MAX_ENTRIES][DIMENSION],
		int numEntries,
		int index_[MAX_ENTRIES],
		float out[MAX_ENTRIES][DIMENSION],
		int epo_code[2][DIMENSION],
		int size,
		int Mi_,                // last cluster
		int bits[3],            // including parity
		CMP_qt type
);


float ep_shaker_2(
		float data[MAX_ENTRIES][MAX_DIMENSION_BIG],
		int numEntries,
		int index_[MAX_ENTRIES],
		float out[MAX_ENTRIES][MAX_DIMENSION_BIG],
		int epo_code[2][MAX_DIMENSION_BIG],
		int size,
		int Mi_,                // last cluster
		int bits, // [3],            // including parity
		//   CMP_qt type,
		int dimension,
		float epo[2][MAX_DIMENSION_BIG]
);

float ep_shaker_(
		float data[MAX_ENTRIES][DIMENSION],
		int numEntries,
		int index_[MAX_ENTRIES],
		float out[MAX_ENTRIES][DIMENSION],
		int epo_code[2][DIMENSION],
		int size,
		int Mi_,                // last cluster
		int bits[3],            // including parity
		CMP_qt type
);


float ep_shaker(
		float data[MAX_ENTRIES][MAX_DIMENSION_BIG],
		int numEntries,
		int index_[MAX_ENTRIES],
		float out[MAX_ENTRIES][MAX_DIMENSION_BIG],
		int epo_code[2][MAX_DIMENSION_BIG],
		// Note:    int size,
		int Mi_,                // last cluster
		int bits[3],            // including parity
		CMP_qt type ,
		int dimension
);

#define MAX_PARITY_CASES     8

#ifdef USE_BC7
extern int npv_nd[][2*MAX_DIMENSION_BIG];
extern int par_vectors_nd[][2*MAX_DIMENSION_BIG][(1<<(2*MAX_DIMENSION_BIG-1))][2][MAX_DIMENSION_BIG];
#endif
