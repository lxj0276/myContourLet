/*
   contourlet - Implementation of the contourlet transform for image coding
   Copyright (C) 2005 Vivien Chappelier - IRISA/University of Rennes 1

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public
   License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU General Public
   License along with this program; if not, write to the Free
   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#include "stdafx.h"
#include "it/types.h"
#include "it/mat.h"
#include <math.h>

#include "libezbc/arithmetic_codec.h"
#include "libezbc/ezbc_types.h"
#include "libezbc/ezbc_tables.h"
#include "libezbc/ezbc_codec.h"
#include "libezbc/ezbc_encoder.h"
#include "libezbc/ezbc_decoder.h"

#include "contourlet.h"

typedef struct _ezbc_subband_ {
  int x, y; /* subband position */
  int w, h; /* subband size */
  ezbc_coeff_t *coeffs; /* subband coefficients (binarized) */
  struct _ezbc_subband_ *parent; /* parent subband */
  struct _ezbc_subband_ *child; /* child subband */
  ezbc_subband_codec_t *codec; /* the associated codec */
  int orientation; /* its orientation */
} ezbc_subband_t;

/* get the maximum amplitude in the subband */
static double find_max(mat image)
{
  size_t w, h;
  size_t x, y;
  double a, max = 0;

  w = mat_width(image);
  h = mat_height(image);

  for(y = 0; y < h; y++) {
    for(x = 0; x < w; x++) {
      a = fabs(image[y][x]);
      if(a > max) max = a;
    }
  }
  return(max);
}

/* extract and binarize coefficients at location +X,Y-WxH */
double binarize_subbands(contourlet_t *ct, ezbc_subband_t *bands, int msb, int bpp)
{
  int i, j;
  int ct_levels = ct->ct_levels;
  int wt_levels = ct->wt_levels;
  ivec dfb_levels = ct->dfb_levels;
  ezbc_coeff_t *coeffs;
  double delta;
  ezbc_coeff_t c;
  double v;
  int w, h;
  int x, y;
  int b;
  double mean = 0;

  w = bands[0].w;
  h = bands[0].h;
  coeffs = bands[0].coeffs;

  /* binarization: use fixed length uniform dead-zone quantization */
  /* warning! the bits are left aligned */
  delta = pow(2.0, msb + 1 - bpp);

  /* compute the mean of the low-frequency band */
  for(y = 0; y < h; y++)
    for(x = 0; x < w; x++)
      mean += ct->dwt[0][y][x];
  mean /= w*h;
  
  /* low-frequency band */
  for(y = 0; y < h; y++)
    for(x = 0; x < w; x++) {
      v = ct->dwt[0][y][x] - mean;

      c = ((int) (fabs(v) / delta)) << (EZBC_MSB + 1 - bpp);

      if(v < 0)	c |= EZBC_SIGN;
      *coeffs++ = c;
    }
  
  b = 1;

  /* wavelet subbands */
  for(i = 0; i < wt_levels; i++)
    for(j = 0; j < 3; j++) {
      w = bands[b].w;
      h = bands[b].h;
      coeffs = bands[b].coeffs;

      for(y = 0; y < h; y++)
	for(x = 0; x < w; x++) {
	  /* transpose the vertical subband so that we can reuse the */
	  /* context model of the horizontal subband */
	  if(bands[b].orientation == EZBC_ORIENT_VERTICAL)
	    v = ct->dwt[1+3*i+j][x][y];
	  else
	    v = ct->dwt[1+3*i+j][y][x];

	  c = ((int) (fabs(v) / delta)) << (EZBC_MSB + 1 - bpp);
	  
	  if(v < 0)	c |= EZBC_SIGN;
	  *coeffs++ = c;
	}

      b++;
    }


  /* contourlet subbands */
  for(i = ct_levels - 1; i >= 0; i--)
    for(j = 0; j < (1 << dfb_levels[i]); j++) {
      w = bands[b].w;
      h = bands[b].h;
      coeffs = bands[b].coeffs;

      for(y = 0; y < h; y++)
	for(x = 0; x < w; x++) {
	  /* transpose the vertical subband so that we can reuse the */
	  /* context model of the horizontal subband */
	  if(bands[b].orientation == EZBC_ORIENT_VERTICAL)
	    v = ct->high[i][j][x][y];
	  else
	    v = ct->high[i][j][y][x];

	  c = ((int) (fabs(v) / delta)) << (EZBC_MSB + 1 - bpp);
	  
	  if(v < 0)	c |= EZBC_SIGN;
	  *coeffs++ = c;
	}

      b++;
    }

  return(mean);
}

/* unbinarize coefficients at location +X,Y-WxH and insert them in the image */
void unbinarize_subband(contourlet_t *ct, ezbc_subband_t *bands, int msb, int bpp, double mean)
{
  int i, j, b;
  int ct_levels = ct->ct_levels;
  int wt_levels = ct->wt_levels;
  ivec dfb_levels = ct->dfb_levels;
  ezbc_coeff_t *coeffs;
  double delta;
  ezbc_coeff_t c;
  double v;
  int w, h;
  int x, y;
  double deltas[EZBC_MSB];

  for(i = 0; i < EZBC_MSB; i++)
    deltas[i] = pow(2.0, msb + 1 - i);

  w = bands[0].w;
  h = bands[0].h;
  coeffs = bands[0].coeffs;

  /* binarization: use fixed length uniform dead-zone quantization */
  /* warning! the bits are left aligned */

  /* low-frequency band */
  for(y = 0; y < h; y++)
    for(x = 0; x < w; x++) {
      c = *coeffs++;

      if(c) {
	bpp = (c & EZBC_LENGTH) + 1;
	delta = deltas[bpp];
	v = ((c & EZBC_ABS) >> (EZBC_MSB + 1 - bpp)) * delta + delta / 2;
      }
      else
	v = 0;

      if(c & EZBC_SIGN) v = -v;

      ct->dwt[0][y][x] = v + mean;
    }

  b = 1;

  for(i = 0; i < wt_levels; i++)
    for(j = 0; j < 3; j++) {
      w = bands[b].w;
      h = bands[b].h;
      coeffs = bands[b].coeffs;

      /* binarization: use fixed length uniform dead-zone quantization */
      /* warning! the bits are left aligned */
      for(y = 0; y < h; y++)
	for(x = 0; x < w; x++) {
	  c = *coeffs++;

	  if(c) {
	    bpp = (c & EZBC_LENGTH) + 1;
	    delta = deltas[bpp];
	    v = ((c & EZBC_ABS) >> (EZBC_MSB + 1 - bpp)) * delta + delta / 2;
	  }
	  else
	    v = 0;

	  if(c & EZBC_SIGN) v = -v;

	  /* transpose the vertical subband so that we can reuse the */
	  /* context model of the horizontal subband */
	  if(bands[b].orientation == EZBC_ORIENT_VERTICAL)
	    ct->dwt[1+3*i+j][x][y] = v;
	  else
	    ct->dwt[1+3*i+j][y][x] = v;
	}
      b++;
    }


  /* contourlet subbands */
  for(i = ct_levels - 1; i >= 0; i--)
    for(j = 0; j < (1 << dfb_levels[i]); j++) {
      w = bands[b].w;
      h = bands[b].h;
      coeffs = bands[b].coeffs;

      /* binarization: use fixed length uniform dead-zone quantization */
      /* warning! the bits are left aligned */
      for(y = 0; y < h; y++)
	for(x = 0; x < w; x++) {
	  c = *coeffs++;

	  if(c) {
	    bpp = (c & EZBC_LENGTH) + 1;
	    delta = deltas[bpp];
	    v = ((c & EZBC_ABS) >> (EZBC_MSB + 1 - bpp)) * delta + delta / 2;
	  }
	  else
	    v = 0;

	  if(c & EZBC_SIGN) v = -v;

	  /* transpose the vertical subband so that we can reuse the */
	  /* context model of the horizontal subband */
	  if(bands[b].orientation == EZBC_ORIENT_VERTICAL)
	    ct->high[i][j][x][y] = v;
	  else
	    ct->high[i][j][y][x] = v;
	}
      b++;
    }
}

/* create the contourlet hierarchy */
void create_contourlet_hierarchy(contourlet_t *ct, ezbc_subband_t *bands)
{
  int i, j;
  int b = 0;
  int ct_levels = ct->ct_levels;
  int wt_levels = ct->wt_levels;
  ivec dfb_levels = ct->dfb_levels;

  /* low-frequency band */
  bands[b].x = 0;
  bands[b].y = 0;
#if 0
  bands[b].w = mat_width(ct->low);
  bands[b].h = mat_height(ct->low);
#else
  bands[b].w = mat_width(ct->dwt[0]);
  bands[b].h = mat_height(ct->dwt[0]);
#endif
  bands[b].child = NULL; /* LL does not really have children */
  bands[b].parent = NULL; /* LL has no parent */
  bands[b].orientation = EZBC_ORIENT_DIAGONAL; // TEMP
  bands[b].coeffs = (ezbc_coeff_t *) calloc(bands[b].w * bands[b].h, sizeof(ezbc_coeff_t));
  b++;

  /* wavelet subbands */
  for(i = 0; i < wt_levels; i++) {
    for(j = 0; j < 3; j++) {
      bands[b].x = bands[b].y = 0; // unused
      switch(j) {
      case 0:
	/* special case as the subband will be transposed */
	bands[b].w = mat_height(ct->dwt[1+3*i+j]);
	bands[b].h = mat_width(ct->dwt[1+3*i+j]);
	bands[b].orientation = EZBC_ORIENT_VERTICAL;
	break;
      case 1:
	bands[b].w = mat_width(ct->dwt[1+3*i+j]);
	bands[b].h = mat_height(ct->dwt[1+3*i+j]);
	bands[b].orientation = EZBC_ORIENT_HORIZONTAL;
	break;
      default:
	bands[b].w = mat_width(ct->dwt[1+3*i+j]);
	bands[b].h = mat_height(ct->dwt[1+3*i+j]);
	bands[b].orientation = EZBC_ORIENT_DIAGONAL;
	break;
      }

      /* don't use inter subband dependency */
      bands[b].parent = bands[b].child = NULL;
      bands[b].coeffs = (ezbc_coeff_t *) calloc(bands[b].w * bands[b].h,
						sizeof(ezbc_coeff_t));
      
      b++;
    }
  }

  /* contourlet subbands */
  for(i = ct_levels - 1; i >= 0; i--) {
    for(j = 0; j < (1 << dfb_levels[i]); j++) {
      bands[b].x = bands[b].y = 0; // unused
      if(j < (1 << (dfb_levels[i] - 1))) {
	/* special case as the subband will be transposed */
	bands[b].w = mat_height(ct->high[i][j]);
	bands[b].h = mat_width(ct->high[i][j]);
	bands[b].orientation = EZBC_ORIENT_VERTICAL;
      } else {
	bands[b].w = mat_width(ct->high[i][j]);
	bands[b].h = mat_height(ct->high[i][j]);
	bands[b].orientation = EZBC_ORIENT_HORIZONTAL;
      }

      /* don't use inter subband dependency */
      bands[b].parent = bands[b].child = NULL;
      bands[b].coeffs = (ezbc_coeff_t *) calloc(bands[b].w * bands[b].h,
						sizeof(ezbc_coeff_t));

      b++;
    }
  }
}

int ezbc_encode(contourlet_t *ct,
		unsigned char *buffer,
		int length,
		int rate)
{
  int i, j, b, count;
  ezbc_subband_t *bands;
  ezbc_subband_codec_t *codec, *child_codec, *parent_codec;
  arithmetic_codec_t *arithmetic_codec;
  double mean;
  int bpp = 26;
  int msb;
  int ct_levels = ct->ct_levels;
  int wt_levels = ct->wt_levels;
  ivec dfb_levels = ct->dfb_levels;

  /* create the contourlet hierarchy */
  count = 1;
  for(i = 0; i < wt_levels; i++)
    count += 3;
  for(i = 0; i < ct_levels; i++)
    count += (1 << dfb_levels[i]);

  bands = (ezbc_subband_t *) malloc(count * sizeof(ezbc_subband_t));
  create_contourlet_hierarchy(ct, bands);

  /* find the most significant bit */
  msb = intlog2((int) ceil(find_max(ct->dwt[0])));
  for(i = 0; i < wt_levels; i++)
    for(j = 0; j < 3; j++)
      if(intlog2((int) ceil(find_max(ct->dwt[1+3*i+j]))) > msb)
	msb = intlog2((int) ceil(find_max(ct->dwt[1+3*i+j])));

  for(i = 0; i < ct_levels; i++)
    for(j = 0; j < (1 << dfb_levels[i]); j++)
      if(intlog2((int) ceil(find_max(ct->high[i][j]))) > msb)
	msb = intlog2((int) ceil(find_max(ct->high[i][j])));

  //TEMP
  msb --;

  *buffer++ = msb; /* code it */
  length--;

  /* binarization: use fixed length uniform dead-zone quantization */
  mean = binarize_subbands(ct, bands, msb, bpp);

  *((double *) buffer) = mean;
  buffer += sizeof(double);
  length -= sizeof(double);

  /* create a new arithmetic encoder */
  arithmetic_codec = arithmetic_codec_new(buffer, length);
  arithmetic_encoder_start(arithmetic_codec);

  /* initialize contexts look up tables */
  ezbc_init();

  /* allocate the subband codecs */
  for(i = 0; i < count; i++)
    bands[i].codec = (ezbc_subband_codec_t *) malloc(sizeof(ezbc_subband_codec_t));

  /* initialize the subband codecs */
  for(i = 0; i < count; i++) {
    codec = bands[i].codec;
    child_codec = (bands[i].child) ? bands[i].child->codec : NULL;
    parent_codec = (bands[i].parent) ? bands[i].parent->codec : NULL;

    ezbc_subband_init(codec, parent_codec, child_codec, arithmetic_codec,
		      bands[i].coeffs,
		      bands[i].w, bands[i].h, bands[i].orientation,
		      bpp, (int) rate);
  }

  /* encode the subbands */
  for(b = 0; b < bpp; b++)
    for(i = 0; i < count; i++)
      ezbc_subband_encode(bands[i].codec, EZBC_MSB - b, NULL);

  //  dump_ezbc_cost();

  /* finish encoding */
  length = arithmetic_encoder_stop(arithmetic_codec);
  arithmetic_codec_delete(arithmetic_codec);

  length += sizeof(double); /* low pass mean */
  length++;                      /* msb */

  /* free resources */
  for(i = 0; i < count; i++) {
    ezbc_subband_cleanup(bands[i].codec);
    free(bands[i].codec);
  }

  free(bands);

  return(length);
}

int ezbc_decode(contourlet_t *ct,
		unsigned char *buffer,
		int length,
		int rate)
{
  int i, b, count;
  ezbc_subband_t *bands;
  ezbc_subband_codec_t *codec, *child_codec, *parent_codec;
  arithmetic_codec_t *arithmetic_codec;
  double mean;
  int bpp = 26;
  int msb;
  int ct_levels = ct->ct_levels;
  int wt_levels = ct->wt_levels;
  ivec dfb_levels = ct->dfb_levels;

  /* create the contourlet hierarchy */
  count = 1;
  for(i = 0; i < wt_levels; i++)
    count += 3;
  for(i = 0; i < ct_levels; i++)
    count += (1 << dfb_levels[i]);
  bands = (ezbc_subband_t *) malloc(count * sizeof(ezbc_subband_t));
  create_contourlet_hierarchy(ct, bands);

  /* decode the index of the most significant bit */
  msb = *buffer++;
  length--;

  /* decode the mean of the low pass subband (float) */
  mean = *((double *) buffer);
  buffer += sizeof(double);
  length -= sizeof(double);

  /* initialize the ezbc codec */
  ezbc_init();

  /* create a new arithmetic decoder */
  arithmetic_codec = arithmetic_codec_new(buffer, length);
  arithmetic_decoder_start(arithmetic_codec);

  /* allocate the subband codecs */
  for(i = 0; i < count; i++)
    bands[i].codec = (ezbc_subband_codec_t *) malloc(sizeof(ezbc_subband_codec_t));

  /* initialize the subband codecs */
  for(i = 0; i < count; i++) {
    codec = bands[i].codec;
    child_codec = (bands[i].child) ? bands[i].child->codec : NULL;
    parent_codec = (bands[i].parent) ? bands[i].parent->codec : NULL;

    ezbc_subband_init(codec, parent_codec, child_codec, arithmetic_codec,
		      bands[i].coeffs,
		      bands[i].w, bands[i].h, bands[i].orientation,
		      bpp, (int) rate);
  }

  /* decode the subbands */
  for(b = 0; b < bpp; b++)
    for(i = 0; i < count; i++)
      ezbc_subband_decode(bands[i].codec, EZBC_MSB - b);

  /* binarization: use fixed length uniform dead-zone quantization */
  unbinarize_subband(ct, bands, msb, bpp, mean);

  /* stop decoding */
  arithmetic_decoder_stop(arithmetic_codec);
  arithmetic_codec_delete(arithmetic_codec);

  /* free resources */
  for(i = 0; i < count; i++) {
    ezbc_subband_cleanup(bands[i].codec);
    free(bands[i].codec);
  }

  free(bands);

  return(length);
}
