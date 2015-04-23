/*
   libit - Library for basic source and channel coding functions
   Copyright (C) 2005-2005 Vivien Chappelier, Herve Jegou

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public
   License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU Library General Public
   License along with this library; if not, write to the Free
   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/ 

/*
  Encoding and Decoding algorithms for Variables Length Codes
  Copyright (C) 2005 Herve Jegou
*/


#ifndef __it_vlc_coding_h
#define __it_vlc_coding_h

#include <it/vlc.h>

#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------*/

/* Return the number of bits required to encode a given source          */
int vlc_nb_bits_required( const vlc_t * vlc, ivec S );

/*----------------------------------------------------------------------*/
/* Encode a symbol flow using concatenation of codewords                */

/* Encode a source with a VLC. The codeword and transmitted concatenated*/
bvec vlc_encode_concat( const vlc_t * vlc, ivec S );

/* Decode a VLC encoded bitstream                                       */
ivec vlc_decode_concat( const vlc_t * vlc, bvec E );

/* The same, but with an expected number of symbols equal to N.
   Hence, the bitstream may be truncated or filled with -1;             */
ivec vlc_decode_concat_N( const vlc_t * vlc, bvec E, size_t N );

/*----------------------------------------------------------------------*/
/* Encode a symbol flow using algorithm CMA                             */

/* Encode a source with a VLC. The codeword and transmitted concatenated*/
bvec vlc_encode_CMA( const vlc_t * vlc, ivec S );


/* Decode a bitstream encoded with vlc_encode_CMA                       */
ivec vlc_decode_CMA( const vlc_t * vlc, bvec E, size_t K );

/*----------------------------------------------------------------------
  Layered encoding and decoding of a symbol flow with VLC               */
 
bvec vlc_encode_layer( const vlc_t * vlc, ivec S );
ivec vlc_decode_layer( const vlc_t * vlc, bvec E, size_t K );

/*----------------------------------------------------------------------
  Encoding and decoding of a symbol flow with algorithm SMA-stack       */

bvec vlc_encode_sma_stack( const vlc_t * vlc, ivec S );
ivec vlc_decode_sma_stack( const vlc_t * vlc, bvec E, size_t K );


#ifdef __cplusplus
}
#endif /* extern "C" */
#endif
