/* libezbc - EZBC subband coding/decoding
      Copyright (C) 2004 Vivien Chappelier

      This library is inspired from code originaly developped by Yongjun Wu Jan in MC-EZBC. This rewritten version may still contain code similar to his original implementation.

      Reference paper: "Embedded Image Coding Using Zeroblocks Of Subband/Wavelet Coefficient And Context Modeling", Shih-Ta Hsiang and John W. Woods, Data Compression Conference (DCC '01), March 27 - 29, 2001, Snowbird, Utah 
*/

#ifndef __BITBUFFER_H
#define __BITBUFFER_H

/* bit storage */
typedef struct _bitbuffer_ {
  unsigned char *data;
  int pos;
  int len;
} bitbuffer_t;

/* create a new bit buffer of 'size' bytes */
bitbuffer_t *bitbuffer_new(unsigned char *buffer, int size);

/* free buffer resources */
void bitbuffer_delete(bitbuffer_t *bb);

/* write a bit to the buffer */
int bitbuffer_write_bit(bitbuffer_t *bb, int bit);

/* read a bit from the buffer */
int bitbuffer_read_bit(bitbuffer_t *bb);

/* write some bits to the buffer */
void bitbuffer_write(bitbuffer_t *bb, int data, int length);

/* read some bits from the buffer */
int bitbuffer_read(bitbuffer_t *bb, int length);

/* byte-align the buffer */
int bitbuffer_flush(bitbuffer_t *bb); /* write */
int bitbuffer_align(bitbuffer_t *bb); /* read */

#define bitbuffer_clear(bb) do { memset(bb->data, 0, bb->len / 8); bitbuffer_rewind(bb); } while(0)

#define bitbuffer_pos(bb) ((bb)->pos)
#define bitbuffer_seek(bb, p) ((bb)->pos = (p))
#define bitbuffer_rewind(bb) bitbuffer_seek(bb, 0)
#define bitbuffer_length(bb) ((bb->pos + 7) >> 3)

#endif
