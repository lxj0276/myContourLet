/* libezbc - EZBC subband coding/decoding
      Copyright (C) 2004 Vivien Chappelier

      This library is inspired from code originaly developped by Yongjun Wu Jan in MC-EZBC. This rewritten version may still contain code similar to his original implementation.

      Reference paper: "Embedded Image Coding Using Zeroblocks Of Subband/Wavelet Coefficient And Context Modeling", Shih-Ta Hsiang and John W. Woods, Data Compression Conference (DCC '01), March 27 - 29, 2001, Snowbird, Utah 
*/

/* Binary Arithmetic Coder/Decoder 
      Copyright (C) 2004 Vivien Chappelier

      Originally inspired from the MPEG-4 VM binary arithmetic coder.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "arithmetic_codec.h"

static void renormalize_encoder(arithmetic_codec_t *arithmetic_codec);
static void renormalize_decoder(arithmetic_codec_t *arithmetic_codec);

arithmetic_codec_t *arithmetic_codec_new(unsigned char *buffer,
					 int length);
void arithmetic_decoder_delete(arithmetic_codec_t *arithmetic_codec);
void arithmetic_codec_start(arithmetic_codec_t *arithmetic_codec);
void arithmetic_decoder_stop(arithmetic_codec_t *arithmetic_codec);
int arithmetic_encoder_stop(arithmetic_codec_t *arithmetic_codec);
int arithmetic_codec_decode(arithmetic_codec_t *arithmetic_codec,
			    float prob_zero);
void arithmetic_codec_encode(arithmetic_codec_t *arithmetic_codec,
			     float prob_zero,
			     arithmetic_codec_bit_t bit);

int ffs(int n)
{
  int i;
  for(i = 0; n; n >>= 1, i++);
  return(i - 1);
}

static void follow(arithmetic_codec_t *arithmetic_codec, unsigned char bit)
{
  if (!arithmetic_codec->first_bit)
    bitbuffer_write_bit(arithmetic_codec->bitbuffer, bit);
  else
    arithmetic_codec->first_bit = 0;

  while(arithmetic_codec->bits_to_follow > 0) {
    bitbuffer_write_bit(arithmetic_codec->bitbuffer, !bit);
    arithmetic_codec->bits_to_follow--;
  }
}

static void renormalize_encoder(arithmetic_codec_t *arithmetic_codec)
{
  while (arithmetic_codec->range < ARITHMETIC_CODEC_QUARTER) {
    if (arithmetic_codec->lower >= ARITHMETIC_CODEC_HALF) { 
      follow(arithmetic_codec, 1);
      arithmetic_codec->lower -= ARITHMETIC_CODEC_HALF;
    } else {
      if (arithmetic_codec->lower + arithmetic_codec->range <=
	  ARITHMETIC_CODEC_HALF)
	follow(arithmetic_codec, 0);
      else {
	arithmetic_codec->bits_to_follow++;
	arithmetic_codec->lower -= ARITHMETIC_CODEC_QUARTER;
      }
    }
    arithmetic_codec->lower += arithmetic_codec->lower;
    arithmetic_codec->range += arithmetic_codec->range;
  }
}

static void renormalize_decoder(arithmetic_codec_t *arithmetic_codec)
{
  while (arithmetic_codec->range < ARITHMETIC_CODEC_QUARTER) {
    if (arithmetic_codec->lower >= ARITHMETIC_CODEC_HALF) {
      arithmetic_codec->value -= ARITHMETIC_CODEC_HALF;
      arithmetic_codec->lower -= ARITHMETIC_CODEC_HALF;
      arithmetic_codec->bits_to_follow = 0;
    } else
      if (arithmetic_codec->lower + arithmetic_codec->range <=
	  ARITHMETIC_CODEC_HALF)
	arithmetic_codec->bits_to_follow = 0;
      else{
	arithmetic_codec->value -= ARITHMETIC_CODEC_QUARTER;
	arithmetic_codec->lower -= ARITHMETIC_CODEC_QUARTER;
	arithmetic_codec->bits_to_follow++;
      }

    arithmetic_codec->lower += arithmetic_codec->lower;
    arithmetic_codec->range += arithmetic_codec->range;
             
    arithmetic_codec->value <<= 1;
    arithmetic_codec->value |= bitbuffer_read_bit(arithmetic_codec->bitbuffer);
  }
}

void arithmetic_encoder_start(arithmetic_codec_t *arithmetic_codec)
{
  arithmetic_codec->value = 0;
  arithmetic_codec->lower = 0;
  arithmetic_codec->range = ARITHMETIC_CODEC_HALF - 1;
  arithmetic_codec->bits_to_follow = 0;
  arithmetic_codec->first_bit = 1;
  bitbuffer_clear(arithmetic_codec->bitbuffer);
}

void arithmetic_decoder_start(arithmetic_codec_t *arithmetic_codec)
{
  int i, j;

  arithmetic_codec->value = 0;

  for (i = 0; i < ARITHMETIC_CODEC_BITS - 1; i++) {
    j = bitbuffer_read_bit(arithmetic_codec->bitbuffer);
    arithmetic_codec->value += arithmetic_codec->value + j;
  }

  arithmetic_codec->lower = 0;
  arithmetic_codec->range = ARITHMETIC_CODEC_HALF - 1;
  arithmetic_codec->bits_to_follow = 0;
  arithmetic_codec->first_bit = 1;
}

void arithmetic_decoder_stop(arithmetic_codec_t *arithmetic_codec)
{
  //  int nbits = 32 - ffs(arithmetic_codec->range);
  int nbits = 32;

  /* put back unused bits */
  bitbuffer_seek(arithmetic_codec->bitbuffer,
		 bitbuffer_pos(arithmetic_codec->bitbuffer) - ARITHMETIC_CODEC_BITS + nbits);

  /* re-align on a byte boundary */
  bitbuffer_align(arithmetic_codec->bitbuffer);
}

int arithmetic_encoder_stop(arithmetic_codec_t *arithmetic_codec)
{
  int nbits, bits, i;

  bits = arithmetic_codec->lower + arithmetic_codec->range / 2;
  /* make sure last symbols are decodable */
  //TEMP  nbits = 32 - ffs(arithmetic_codec->range);
  nbits = 32;

  for (i = 0; i < nbits-1; i++)
    follow(arithmetic_codec, (bits >> (nbits - 1 - i)) & 1);

  /* re-align on a byte boundary */
  return(bitbuffer_flush(arithmetic_codec->bitbuffer));
}


arithmetic_codec_t *arithmetic_codec_new(unsigned char *buffer,
					 int length)
{
  arithmetic_codec_t *arithmetic_codec;

  arithmetic_codec = (arithmetic_codec_t *) malloc(sizeof(arithmetic_codec_t));
  memset(arithmetic_codec, 0, sizeof(arithmetic_codec_t));
  arithmetic_codec->sequence = buffer;
  arithmetic_codec->length = length;
  arithmetic_codec->bitbuffer = bitbuffer_new(buffer, length);

  return(arithmetic_codec);
}

void arithmetic_codec_delete(arithmetic_codec_t *arithmetic_codec)
{
  bitbuffer_delete(arithmetic_codec->bitbuffer);
  free(arithmetic_codec);
}

int arithmetic_codec_decode(arithmetic_codec_t *arithmetic_codec,
			      float prob_0) {
  int bit;
  float prob_1 = 1.0 - prob_0;
  int LPS = prob_0 > prob_1;
  float pLPS = LPS ? prob_1 : prob_0;
  arithmetic_codec_register_t rLPS;

  rLPS = (arithmetic_codec_register_t) (arithmetic_codec->range * pLPS);

  if (arithmetic_codec->value - arithmetic_codec->lower >=
      arithmetic_codec->range - rLPS) {
    bit = LPS;
    arithmetic_codec->lower += arithmetic_codec->range - rLPS;
    arithmetic_codec->range = rLPS;
  } else {
    bit = !LPS;
    arithmetic_codec->range -= rLPS;
  }

  renormalize_decoder(arithmetic_codec);

  return(bit);
}


void arithmetic_codec_encode(arithmetic_codec_t *arithmetic_codec,
			     float prob_0,
			     arithmetic_codec_bit_t bit)
{
  float prob_1 = 1.0 - prob_0;
  int LPS = prob_0 > prob_1;
  float pLPS = LPS ? prob_1 : prob_0;
  arithmetic_codec_register_t rLPS;

  rLPS = (arithmetic_codec_register_t) (arithmetic_codec->range * pLPS);

  if (bit == LPS) {
    arithmetic_codec->lower += arithmetic_codec->range - rLPS;
    arithmetic_codec->range = rLPS;
  } else
    arithmetic_codec->range -= rLPS;

  renormalize_encoder(arithmetic_codec);
}


