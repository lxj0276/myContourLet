/* libezbc - EZBC subband coding/decoding
      Copyright (C) 2004 Vivien Chappelier

      This library is inspired from code originaly developped by Yongjun Wu Jan in MC-EZBC. This rewritten version may still contain code similar to his original implementation.

      Reference paper: "Embedded Image Coding Using Zeroblocks Of Subband/Wavelet Coefficient And Context Modeling", Shih-Ta Hsiang and John W. Woods, Data Compression Conference (DCC '01), March 27 - 29, 2001, Snowbird, Utah 
*/

/* Binary Input/Output
      Copyright (C) 2004 Vivien Chappelier
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "bitbuffer.h"

bitbuffer_t *bitbuffer_new(unsigned char *buffer, int size)
{
  bitbuffer_t *bb;

  bb = (bitbuffer_t *) malloc(sizeof(bitbuffer_t));

  bb->data = buffer;
  bb->pos = 0;
  bb->len = 8 * size;

  return(bb);
}

void bitbuffer_delete(bitbuffer_t *bb)
{
  free(bb);
}

int bitbuffer_write_bit(bitbuffer_t *bb, int bit)
{
  if(bb->pos >= bb->len) return(-1);
  bb->data[bb->pos >> 3] |= bit << (7 - (bb->pos & 7));
  bb->pos++;
  return(1);
}

int bitbuffer_read_bit(bitbuffer_t *bb)
{
  int bit;

  if(bb->pos >= bb->len) return(-1);
  bit = (bb->data[bb->pos >> 3] >> (7 - (bb->pos & 7))) & 1;
  bb->pos++;
  return(bit);
}

void bitbuffer_write(bitbuffer_t *bb, int data, int length)
{
  int i;

  for(i = 0; i < length; i++)
    bitbuffer_write_bit(bb, (data >> (length - 1 - i)) & 1);
}

int bitbuffer_read(bitbuffer_t *bb, int length)
{
  int i;
  int data;

  data = 0;
  for(i = 0; i < length; i++) {
    data <<= 1;
    data |= bitbuffer_read_bit(bb);
  }

  return(data);
}

int bitbuffer_flush(bitbuffer_t *bb)
{
  bb->pos = (bb->pos + 7) & (~7);
  return(bb->pos >> 3);
}

int bitbuffer_align(bitbuffer_t *bb)
{
  bb->pos = (bb->pos + 7) & (~7);
  return(bb->pos >> 3);
}

