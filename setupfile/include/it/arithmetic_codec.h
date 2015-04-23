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
  Arithmetic coder/decoder
  Copyright (C) 2005 Vivien Chappelier
*/

/*---------------------------------------------------------------------------*/
/*! \defgroup entropycoding                                                  */
/* @{                                                                        */
/*---------------------------------------------------------------------------*/

#ifndef __it_arithmetic_codec_h
#define __it_arithmetic_codec_h

typedef unsigned char arithmetic_codec_bit_t;
typedef unsigned long arithmetic_codec_register_t;

typedef struct _arithmetic_coder_t_ {
  int  precision;                       /* # bits of precision */
  arithmetic_codec_register_t lower;	/* lower bound of the interval */
  arithmetic_codec_register_t range;	/* range of the interval */
  int  pending;                         /* # pending bits to be written */
  bvec buffer;	                        /* buffer for the sequence of bits */
  arithmetic_codec_bit_t *sequence;	/* pointer to the current bit */
} arithmetic_coder_t;

typedef struct _arithmetic_decoder_t_ {
  int  precision;                       /* # bits of precision */
  arithmetic_codec_register_t lower;	/* lower bound of the interval */
  arithmetic_codec_register_t range;	/* range of the interval */
  arithmetic_codec_register_t value;	/* encoded value */
  bvec buffer;	                        /* sequence of bits */
  arithmetic_codec_bit_t *sequence;	/* pointer to the current bit */
} arithmetic_decoder_t;


/* encoder constructor/destructor */
arithmetic_coder_t *arithmetic_coder_new (int precision);
void arithmetic_coder_delete (arithmetic_coder_t * arithmetic_coder);

/* start encoding */
void arithmetic_coder_start (arithmetic_coder_t * arithmetic_coder,
			     bvec buffer);
/* stop encoding and return the number of bits written to the buffer */
int  arithmetic_coder_stop (arithmetic_coder_t * arithmetic_coder);

/* encode one bit with probability prob_0 of being 0 */
void arithmetic_coder_encode_bit (arithmetic_coder_t * arithmetic_coder,
				  double prob_0, arithmetic_codec_bit_t bit);
/* encode a symbol with the binary arithmetic coder using unary binarization */
void arithmetic_coder_encode_symbol (arithmetic_coder_t * arithmetic_coder,
				     vec pdf, int symbol);

/* decoder constructor/destructor */
arithmetic_decoder_t *arithmetic_decoder_new (int precision);
void arithmetic_decoder_delete (arithmetic_decoder_t * arithmetic_decoder);

/* start decoding */
void arithmetic_decoder_start (arithmetic_decoder_t * arithmetic_decoder,
			       bvec buffer);
/* stop decoding and return the number of bits read from the buffer */
int  arithmetic_decoder_stop (arithmetic_decoder_t * arithmetic_decoder);

/* decode one bit with probability prob_zero of being 0 */
arithmetic_codec_bit_t arithmetic_decoder_decode_bit (arithmetic_decoder_t *
						      arithmetic_decoder,
						      double prob_zero);

/* decode a symbol with the binary arithmetic coder using unary binarization */
int  arithmetic_decoder_decode_symbol (arithmetic_decoder_t *
				       arithmetic_decoder, vec pdf);

/*@}*/

#endif
