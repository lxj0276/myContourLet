/* libezbc - EZBC subband coding/decoding
      Copyright (C) 2004 Vivien Chappelier

      This library is inspired from code originaly developped by Yongjun Wu Jan in MC-EZBC. This rewritten version may still contain code similar to his original implementation.

      Reference paper: "Embedded Image Coding Using Zeroblocks Of Subband/Wavelet Coefficient And Context Modeling", Shih-Ta Hsiang and John W. Woods, Data Compression Conference (DCC '01), March 27 - 29, 2001, Snowbird, Utah 
*/

#ifndef __ARITHMETIC_CODEC_H
#define __ARITHMETIC_CODEC_H

#include "bitbuffer.h"

typedef unsigned char arithmetic_codec_bit_t;
typedef unsigned long arithmetic_codec_register_t;

#define ARITHMETIC_CODEC_BITS    (8*sizeof(arithmetic_codec_register_t))
#define ARITHMETIC_CODEC_HALF    (1UL << (ARITHMETIC_CODEC_BITS-1))
#define ARITHMETIC_CODEC_QUARTER (1UL << (ARITHMETIC_CODEC_BITS-2))

#define ARITHMETIC_CODEC_EOS 0xff /* end of sequence */

typedef struct _arithmetic_codec_t_ {
  arithmetic_codec_register_t lower;
  arithmetic_codec_register_t range;
  arithmetic_codec_register_t value;
  int bits_to_follow;
  int first_bit;
  int length;
  arithmetic_codec_bit_t *sequence;
  bitbuffer_t *bitbuffer;
} arithmetic_codec_t;

/* create a new arithmetic codec reading from/writing to buffer */
arithmetic_codec_t *arithmetic_codec_new(unsigned char *buffer, int length);

/* delete the arithmetic codec object */
void arithmetic_codec_delete(arithmetic_codec_t *arithmetic_codec);

/* start arithmetic coding */
void arithmetic_encoder_start(arithmetic_codec_t *arithmetic_codec);

/* start arithmetic decoding */
void arithmetic_decoder_start(arithmetic_codec_t *arithmetic_codec);

/* stop arithmetic coding, returning the length of the coded sequence */
int arithmetic_encoder_stop(arithmetic_codec_t *arithmetic_codec);

/* stop arithmetic decoding */
void arithmetic_decoder_stop(arithmetic_codec_t *arithmetic_codec);

/* encode a binary symbol knowing the probability of 0 */
void arithmetic_codec_encode(arithmetic_codec_t *arithmetic_codec,
			     float prob_0,
			     arithmetic_codec_bit_t bit);

/* decode a binary symbol knowing the probability of 0 */
int arithmetic_codec_decode(arithmetic_codec_t *arithmetic_codec,
			    float prob_zero);


#endif
