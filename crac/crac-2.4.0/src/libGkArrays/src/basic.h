/* basics.h

   Code mostly inspired from original file written by Rodrigo Gonzalez.

   Copyright (C) 2005, Rodrigo Gonzalez, all rights reserved.

   Some preliminary stuff

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/


#ifndef BASICSINCLUDED
#define BASICSINCLUDED

#include <limits>
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

/*number of bits of a machine word */
const size_t W = std::numeric_limits<u_int64_t>::digits;

/* reads bit p from e */
template <class T>
bool bitget(const T * const e, const u_int64_t p) {
  return ((e[p/W] >> (p % W)) & u_int64_t(1));
}

/* sets bit p in e */
template <class T>
void bitset(T * const e, const u_int64_t p) {
  e[p/W] |= u_int64_t(1) << (p % W);
}

/* cleans bit p in e */
template <class T>
void bitclean(T * const e, const u_int64_t p) {
  e[p/W] &= ~(u_int64_t(1) << (p % W));
}

/* number of machine words requiered to represents e elements of length n */
inline u_int64_t nb_digits(const u_int64_t &e, const u_int64_t &n) {
  return (e * n) / W + ((e * n) % W > 0);
}

#endif

// Local Variables:
// mode:c++
// End:
