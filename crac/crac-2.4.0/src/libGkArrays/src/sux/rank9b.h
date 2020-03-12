/*		 
 * Sux: Succinct data structures
 *
 * Copyright (C) 2007-2010 Sebastiano Vigna 
 *
 *  This library is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as published by the Free
 *  Software Foundation; either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  This library is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
 *  for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef rank9b_h
#define rank9b_h
#include <stdint.h>
#include <assert.h>
#include "popcount.h"
#include "macros.h"

class rank9b {
private:
	const uint64_t *bits;
	uint64_t *counts, *inventory;
	uint64_t num_words, num_counts, inventory_size, ones_per_inventory, log2_ones_per_inventory, num_ones;

	/** Counts the number of bits in x. */
	__inline static int count( const uint64_t x ) {
		register uint64_t byte_sums = x - ( ( x & 0xa * ONES_STEP_4 ) >> 1 );
		byte_sums = ( byte_sums & 3 * ONES_STEP_4 ) + ( ( byte_sums >> 2 ) & 3 * ONES_STEP_4 );
		byte_sums = ( byte_sums + ( byte_sums >> 4 ) ) & 0x0f * ONES_STEP_8;
		return byte_sums * ONES_STEP_8 >> 56;
	}

	/* Selects the k-th (k>=0) bit in l. Returns 72 if no such bit exists. */

	__inline static int select_in_word( const uint64_t x, const int k ) {
		assert( k < count( x ) );

#ifdef SELPOPCOUNT
		for( int i = 0, c = k; i < 64; i+=8 )
			if ( ( c -= popcount[ x >> i & 0xFF ] ) < 0 ) {
				c += popcount[ x >> i & 0xFF ];
				for( int j = 0; j < 8; j++ ) if ( ( x & 1ULL << ( i + j ) ) && c-- == 0 ) return i + j;
			}
		return -1;
#endif

		// Phase 1: sums by byte
		register uint64_t byte_sums = x - ( ( x & 0xa * ONES_STEP_4 ) >> 1 );
		byte_sums = ( byte_sums & 3 * ONES_STEP_4 ) + ( ( byte_sums >> 2 ) & 3 * ONES_STEP_4 );
		byte_sums = ( byte_sums + ( byte_sums >> 4 ) ) & 0x0f * ONES_STEP_8;
		byte_sums *= ONES_STEP_8;

		// Phase 2: compare each byte sum with k
		const uint64_t k_step_8 = k * ONES_STEP_8;
		const uint64_t place = ( LEQ_STEP_8( byte_sums, k_step_8 ) * ONES_STEP_8 >> 53 ) & ~0x7;

		// Phase 3: Locate the relevant byte and make 8 copies with incrental masks
		const int byte_rank = k - ( ( ( byte_sums << 8 ) >> place ) & 0xFF );

		const uint64_t spread_bits = ( x >> place & 0xFF ) * ONES_STEP_8 & INCR_STEP_8;
		const uint64_t bit_sums = ZCOMPARE_STEP_8( spread_bits ) * ONES_STEP_8;

		// Compute the inside-byte location and return the sum
		const uint64_t byte_rank_step_8 = byte_rank * ONES_STEP_8;

		return place + ( LEQ_STEP_8( bit_sums, byte_rank_step_8 ) * ONES_STEP_8 >> 56 );
	}


public:
	rank9b( const uint64_t * const bits, const uint64_t num_bits );
	~rank9b();
	uint64_t rank( const uint64_t pos );
	uint64_t select( const uint64_t rank );
	// Just for analysis purposes
	void print_counts();
	uint64_t bit_count();
};

#endif
