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

#include <cstdio>
#include <ctime>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <sys/time.h>
#include <sys/resource.h>
#include "rank9b.h"

uint64_t shifts, dichotomous;

rank9b::rank9b( const uint64_t * const bits, const uint64_t num_bits ) {
	this->bits = bits;
	num_words = ( num_bits + 63 ) / 64;
	num_counts = ( ( num_bits + 64 * 8 - 1 ) / ( 64 * 8 ) ) * 2;
	
	// Init rank/select structure
	counts = new uint64_t[ num_counts + 1 ];
	memset( counts, 0, ( num_counts + 1 ) * sizeof *counts );

	uint64_t c = 0;
	uint64_t pos = 0;
	for( uint64_t i = 0; i < num_words; i += 8, pos += 2 ) {
		counts[ pos ] = c;
		c += count( bits[ i ] );
		for( int j = 1;  j < 8; j++ ) {
			counts[ pos + 1 ] |= ( c - counts[ pos ] ) << 9 * ( j - 1 );
			if ( i + j < num_words ) c += count( bits[ i + j ] );
		}
	}

	counts[ num_counts ] = c;

	assert( c <= num_bits );

#ifndef NOINVENTORY

	ones_per_inventory = ( c * 16 * 64 + num_bits - 1 ) / num_bits;
	// Make ones_per_inventory into a power of 2
	for( log2_ones_per_inventory = 64; log2_ones_per_inventory-- != 1; ) if ( ones_per_inventory & 1ULL << log2_ones_per_inventory ) break;
	ones_per_inventory = 1ULL << log2_ones_per_inventory;
	inventory_size = ( c + ones_per_inventory - 1 ) / ones_per_inventory;

	assert( ones_per_inventory <= 16 * 64 );

	inventory = new uint64_t[ inventory_size + 1 ];
	uint64_t d = 0;
	const uint64_t mask = ones_per_inventory - 1;
	for( uint64_t i = 0; i < num_words; i++ )
		for( int j = 0; j < 64; j++ )
			if ( bits[ i ] & 1ULL << j ) {
				if ( ( d & mask ) == 0 ) {
					inventory[ d >> log2_ones_per_inventory ] = ( i / 8 ) * 2;
					assert( counts[ ( i / 8 ) * 2 ] <= d ); 
					assert( counts[ ( i / 8 ) * 2 + 2 ] > d ); 
				}

				d++;
			}

	assert( c == d );
	inventory[ inventory_size ] = ( num_words / 8 ) * 2;


#ifdef DEBUG
	printf("First inventories: %lld %lld %lld %lld\n", inventory[ 0 ], inventory[ 1 ], inventory[ 2 ], inventory[ 3 ] );
#endif

#endif

#ifndef NDEBUG
	uint64_t r, t;
	for( uint64_t i = 0; i < c; i++ ) {
		t = select( i );
		r = rank( t );
		if ( r != i ) {
			printf( "i: %lld s: %lld r: %lld\n", i, t, r );
			assert( r == i );
		}
	}

	for( uint64_t i = 0; i < num_bits; i++ ) {
		r = rank( i );
		if ( r < c ) {
			t = select( r );
			if ( t < i ) {
				printf( "i: %lld r: %lld s: %lld\n", i, r, t );
				assert( t >= i );
			}
		}
	}
#endif
}

rank9b::~rank9b() {
	delete [] counts;
#ifndef NOINVENTORY
	delete [] inventory;
#endif
}


uint64_t rank9b::rank( const uint64_t k ) {
	const uint64_t word = k / 64;
	const uint64_t block = word / 4 & ~1;
	const int offset = word % 8 - 1;

#ifdef RANKPOPCOUNT
	uint64_t b = bits[ word ] & ( ( 1ULL << k % 64 ) - 1 );
	int c = 0;
	for( int i = ( k % 64 + 7 ) / 8; i-- != 0; ) {
			c += popcount[ b & 0xFF ];
			b >>= 8;
	}
	return counts[ block ] + ( counts[ block + 1 ] >> ( offset + ( offset >> sizeof offset * 8 - 4 & 0x8 ) ) * 9 & 0x1FF ) + c;
#else
#ifdef RANKPOPCOUNT2
	const uint64_t b = bits[ word ] & ( ( 1ULL << k % 64 ) - 1 );
	return counts[ block ] + ( counts[ block + 1 ] >> ( offset + ( offset >> sizeof offset * 8 - 4 & 0x8 ) ) * 9 & 0x1FF ) + 
		popcount[ b & 0xFF ] + popcount[ b >> 8 & 0xFF ] + popcount[ b >> 16 & 0xFF ] + popcount[ b >> 24 & 0xFF ] +
		popcount[ b >> 32 & 0xFF ] + popcount[ b >> 40 & 0xFF ] + popcount[ b >> 48 & 0xFF ] + popcount[ b >> 56 & 0xFF ];
#else
	return counts[ block ] + ( counts[ block + 1 ] >> ( offset + ( offset >> sizeof offset * 8 - 4 & 0x8 ) ) * 9 & 0x1FF ) + count( bits[ word ] & ( ( 1ULL << k % 64 ) - 1 ) );
#endif
#endif
}



uint64_t rank9b::select( const uint64_t rank ) {
#ifndef NOINVENTORY
	const uint64_t inventory_index_left = rank >> log2_ones_per_inventory;
	assert( inventory_index_left < inventory_size );

	uint64_t block_left = inventory[ inventory_index_left ];
	uint64_t block_right = inventory[ inventory_index_left + 1 ];
#else
	uint64_t block_left = 0, block_right = num_words / 8 * 2;
#endif
#ifdef DEBUG
	printf( "Initially,    rank: %lld block_left: %lld block_right: %lld counts_left: %lld counts_right: %lld\n", rank, block_left, block_right, counts[ block_left ], counts[ block_right ] );
#endif

	if ( rank >= counts[ block_right ] ) {
#ifdef COUNTS
		shifts++;
#endif
		block_right = ( block_left = block_right ) + 2;
#ifdef DEBUG
printf( "After shift,  rank: %lld block_left: %lld block_right: %lld counts_left: %lld counts_right: %lld\n", rank, block_left, block_right, counts[ block_left ], counts[ block_right ] );
#endif
	}
	else {
		uint64_t block_middle;

		while( block_right - block_left > 2 ) {
			block_middle = ( block_right + block_left ) / 2 & ~1;
			if ( rank >= counts[ block_middle ] ) block_left = block_middle;
			else block_right = block_middle;
#ifdef COUNTS
			dichotomous++;
#endif
		}
	}

	assert( counts[ block_left ] <= rank );
	assert( counts[ block_right ] > rank );

	const uint64_t rank_in_block = rank - counts[ block_left ];	

#ifdef DEBUG
	printf( "After search, rank: %lld block_left: %lld block_right: %lld counts_left: %lld counts_right: %lld rank_in_block: %lld\n", rank, block_left, block_right, counts[ block_left ], counts[ block_right ], rank_in_block );
#endif

	assert( counts[ block_left ] <= rank );
	assert( counts[ block_right ] > rank );
	assert( block_right - block_left == 2 );
	assert( rank_in_block < counts[ block_right ] - counts[ block_left ] );
	assert( rank_in_block < 512 );
	
	const uint64_t rank_in_block_step_9 = rank_in_block * ONES_STEP_9;
	const uint64_t subcounts = counts[ block_left + 1 ];
	const uint64_t offset_in_block = ( ULEQ_STEP_9( subcounts, rank_in_block_step_9 ) * ONES_STEP_9 >> 54 & 0x7 );

	const uint64_t word = block_left * 4 + offset_in_block;
	const uint64_t rank_in_word = rank_in_block - ( subcounts >> ( offset_in_block - 1 & 7 ) * 9 & 0x1FF );
#ifdef DEBUG
	printf( "rank_in_block: %lld offset_in_block: %lld rank_in_word: %lld compare: %016llx shift: %lld\n", rank_in_block, offset_in_block, rank_in_word, UCOMPARE_STEP_9( rank_in_block_step_9, subcounts ), subcounts >> ( offset_in_block - 1 ) * 9 & 0x1FF );
#endif
	assert( offset_in_block >= 0 );
	assert( offset_in_block <= 7 );

#ifdef DEBUG
	printf( "rank_in_block: %lld offset_in_block: %lld rank_in_word: %lld\n", rank_in_block, offset_in_block, rank_in_word );
	printf( "subcounts: " );
	for( int i = 0; i < 7; i++ ) printf( "%lld ", subcounts >> i * 9 & 0x1FF );
	printf( "\n" );
	fflush( stdout );
#endif

	assert( rank_in_word < 64 );
	assert( rank_in_word >= 0 );

	return word * 64 + select_in_word( bits[ word ], rank_in_word );
}

uint64_t rank9b::bit_count() {
#ifdef NOINVENTORY
	return num_counts * 64;
#else
	return ( num_counts + inventory_size ) * 64;
#endif
}

void rank9b::print_counts() {}
