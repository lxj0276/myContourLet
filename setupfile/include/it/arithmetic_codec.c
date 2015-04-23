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

#include <it/types.h>
#include <it/vec.h>
#include <it/arithmetic_codec.h>

/* output bits from the register to the bitstream */
static void write_bit (arithmetic_coder_t * arithmetic_coder,
		       unsigned char bit)
{
  *arithmetic_coder->sequence++ = bit;

  while (arithmetic_coder->pending > 0) {
    *arithmetic_coder->sequence++ = !bit;
    arithmetic_coder->pending--;
  }
}

/* input bits from the bitstream to the register */
static void read_bit (arithmetic_decoder_t * arithmetic_decoder)
{
  int  i;

  arithmetic_decoder->sequence++;

  i = arithmetic_decoder->sequence[arithmetic_decoder->precision - 1] & 1;
  arithmetic_decoder->value += arithmetic_decoder->value + i;
}

/* renormalize the encoder */
static void renormalize_enc (arithmetic_coder_t * arithmetic_coder)
{
  arithmetic_codec_register_t half = 1UL << (arithmetic_coder->precision - 1);
  arithmetic_codec_register_t quarter =
    1UL << (arithmetic_coder->precision - 2);

  while (arithmetic_coder->range < quarter) {
    if (arithmetic_coder->lower >= half) {
      /* range entirely contained in upper half */
      write_bit (arithmetic_coder, 1);
      arithmetic_coder->lower -= half;
    }
    else if (arithmetic_coder->lower + arithmetic_coder->range <= half) {
      /* range entirely contained in lower half */
      write_bit (arithmetic_coder, 0);
    }
    else {
      /* bit cannot be determined yet */
      arithmetic_coder->pending++;
      arithmetic_coder->lower -= quarter;
    }

    arithmetic_coder->lower <<= 1;
    arithmetic_coder->range <<= 1;
  }
}

/* renormalize the decoder */
static void renormalize_dec (arithmetic_decoder_t * arithmetic_decoder)
{
  arithmetic_codec_register_t half =
    1UL << (arithmetic_decoder->precision - 1);
  arithmetic_codec_register_t quarter =
    1UL << (arithmetic_decoder->precision - 2);

  while (arithmetic_decoder->range < quarter) {
    if (arithmetic_decoder->lower >= half) {
      /* range entirely contained in upper half */
      arithmetic_decoder->value -= half;
      arithmetic_decoder->lower -= half;
    }
    else if (arithmetic_decoder->lower + arithmetic_decoder->range <= half) {
      /* range entirely contained in lower half */
      /* nothing to do */
    }
    else {
      arithmetic_decoder->value -= quarter;
      arithmetic_decoder->lower -= quarter;
    }

    arithmetic_decoder->lower <<= 1;
    arithmetic_decoder->range <<= 1;

    read_bit (arithmetic_decoder);
  }
}

void arithmetic_coder_start (arithmetic_coder_t * arithmetic_coder,
			     bvec buffer)
{
  arithmetic_coder->buffer = buffer;
  arithmetic_coder->lower = 0;
  arithmetic_coder->range = (1UL << arithmetic_coder->precision) - 1;
  arithmetic_coder->pending = 0;
  arithmetic_coder->sequence = arithmetic_coder->buffer;
}

int arithmetic_coder_stop (arithmetic_coder_t * arithmetic_coder)
{
  int  i;
  arithmetic_codec_register_t bits;

  /* flush the register */
  bits = arithmetic_coder->lower + 1;
  for (i = 0; i < arithmetic_coder->precision; i++)
    write_bit (arithmetic_coder,
	       (unsigned
		char) ((bits >> (arithmetic_coder->precision - 1 - i))
		       & 1));

  return (arithmetic_coder->sequence - arithmetic_coder->buffer);
}

arithmetic_coder_t *arithmetic_coder_new (int precision)
{
  arithmetic_coder_t *arithmetic_coder;

  assert (precision <= 8 * sizeof (arithmetic_codec_register_t));
  arithmetic_coder =
    (arithmetic_coder_t *) malloc (sizeof (arithmetic_coder_t));
  memset (arithmetic_coder, 0, sizeof (arithmetic_coder_t));
  arithmetic_coder->precision = precision;
  arithmetic_coder->buffer = 0;

  return (arithmetic_coder);
}

void arithmetic_coder_delete (arithmetic_coder_t * arithmetic_coder)
{
  free (arithmetic_coder);
}

void arithmetic_coder_encode_bit (arithmetic_coder_t * arithmetic_coder,
				  double prob_0, arithmetic_codec_bit_t bit)
{
  double prob_1 = 1.0 - prob_0;
  int  LPS = prob_0 > prob_1;
  double prob_LPS = LPS ? prob_1 : prob_0;
  arithmetic_codec_register_t range_LPS;

  range_LPS =
    (arithmetic_codec_register_t) (arithmetic_coder->range * prob_LPS);
  if (!range_LPS)
    range_LPS = 1;		/* avoid ranges below the coder precision */

  if (bit == LPS) {
    arithmetic_coder->lower += arithmetic_coder->range - range_LPS;
    arithmetic_coder->range = range_LPS;
  }
  else
    arithmetic_coder->range -= range_LPS;

  renormalize_enc (arithmetic_coder);
}


void arithmetic_decoder_start (arithmetic_decoder_t * arithmetic_decoder,
			       bvec buffer)
{
  int  i, j;

  arithmetic_decoder->buffer = buffer;
  arithmetic_decoder->sequence = buffer;
  arithmetic_decoder->value = 0;

  for (i = 0; i < arithmetic_decoder->precision; i++) {
    j = arithmetic_decoder->sequence[i] & 1;
    arithmetic_decoder->value += arithmetic_decoder->value + j;
  }

  arithmetic_decoder->lower = 0;
  arithmetic_decoder->range = (1UL << arithmetic_decoder->precision) - 1;
}

int arithmetic_decoder_stop (arithmetic_decoder_t * arithmetic_decoder)
{
  int  i;

  /* flush the register */
  for (i = 0; i < arithmetic_decoder->precision; i++)
    read_bit (arithmetic_decoder);

  return (arithmetic_decoder->sequence - arithmetic_decoder->buffer);
}

arithmetic_decoder_t *arithmetic_decoder_new (int precision)
{
  arithmetic_decoder_t *arithmetic_decoder;

  assert (precision <= 8 * sizeof (arithmetic_codec_register_t));
  arithmetic_decoder =
    (arithmetic_decoder_t *) malloc (sizeof (arithmetic_decoder_t));
  memset (arithmetic_decoder, 0, sizeof (arithmetic_decoder_t));
  arithmetic_decoder->precision = precision;
  return (arithmetic_decoder);
}

void arithmetic_decoder_delete (arithmetic_decoder_t * arithmetic_decoder)
{
  free (arithmetic_decoder);
}

arithmetic_codec_bit_t arithmetic_decoder_decode_bit (arithmetic_decoder_t *
						      arithmetic_decoder,
						      double prob_0)
{
  int  bit;
  double prob_1 = 1.0 - prob_0;
  int  LPS = prob_0 > prob_1;
  double prob_LPS = LPS ? prob_1 : prob_0;
  arithmetic_codec_register_t range_LPS;

  range_LPS =
    (arithmetic_codec_register_t) (arithmetic_decoder->range * prob_LPS);
  if (!range_LPS)
    range_LPS = 1;		/* avoid ranges below the coder precision */

  if (arithmetic_decoder->value - arithmetic_decoder->lower >=
      arithmetic_decoder->range - range_LPS) {
    bit = LPS;
    arithmetic_decoder->lower += arithmetic_decoder->range - range_LPS;
    arithmetic_decoder->range = range_LPS;
  }
  else {
    bit = !LPS;
    arithmetic_decoder->range -= range_LPS;
  }

  renormalize_dec (arithmetic_decoder);
  return (bit);
}

/* encode a symbol with the binary arithmetic coder using unary binarization */
/* BUG: this may fail if the probabilities are too close to zero */
/* FIXME: force the probabilities to be > 1 / (2^precision) */
void arithmetic_coder_encode_symbol (arithmetic_coder_t * arithmetic_coder,
				     vec pdf, int symbol)
{
  int  i, l;
  double p0, p1;

  l = vec_length (pdf);
  assert (l > 1);
  assert (symbol < l);

  p0 = 0;
  p1 = 1;

  for (i = 0; i < l; i++) {

    p0 = pdf[i];
    p1 -= pdf[i];

    if (symbol == i) {
      arithmetic_coder_encode_bit (arithmetic_coder, p0 / (p1 + p0), 0);
      break;
    }
    else {
      arithmetic_coder_encode_bit (arithmetic_coder, p0 / (p1 + p0), 1);
    }
  }
}

/* decode a symbol with the binary arithmetic coder using unary binarization */
int arithmetic_decoder_decode_symbol (arithmetic_decoder_t *
				      arithmetic_decoder, vec pdf)
{
  int  i, l;
  double p0, p1;

  l = vec_length (pdf);
  assert (l > 1);

  p0 = 0;
  p1 = 1;

  for (i = 0; i < l; i++) {

    p0 = pdf[i];
    p1 -= pdf[i];

    if (arithmetic_decoder_decode_bit (arithmetic_decoder, p0 / (p1 + p0))
	== 0)
      return (i);
  }
  return (-1);
}
