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
  Separable 2D wavelet transform using lifting.
  Copyright (C) 2005 Vivien Chappelier, Herve Jegou
*/

#include <it/types.h>
#include <it/wavelet.h>
#include <it/wavelet2D.h>
#include <it/io.h>

/* This computes the wavelet decomposition using the lifting implementation.
   It uses a buffer of exactly twice the size of the image. As the memory 
   layout may look a bit complex, these comments address
   how one stage of the decomposition is done. 

   1. The image (or low-pass band of the previous level) to decompose is
   initially in the upper half of the buffer, in other terms, at the end.

   2. First the horizontal (H) decomposition is done which results in a low and
   high band. The low band is stored at the very beginning of the buffer and
   has half the number of coefficients of the initial image. The high band is
   stored just behind and has the same number of coefficients. The initial 
   image is not needed anymore and can now be discarded (overwritten).

   3. The second decomposition is a vertical (V) decomposition on the horizontal
   high band. It computes the horizontal high vertical low (high_low) and
   horizontal high vertical high (high_high) bands which both have one quarter
   the number of coefficients of the original image and are stored at the 
   beginning of the second half of the buffer (overwritting partially the
   original image). The horizontal high band is discarded.

   4. The final decomposition is a vertical decomposition too, but on the
   horizontal low band. The low_low and low_high bands resulting from the
   decomposition are stored in the second quarter of the buffer, overwritting
   the old horizontal high band. The second and third eighth of the buffer can
   now be seen as a new buffer of quarter the size of the initial one,
   containing the low pass image to process in its upper half.

   The decomposition can therefore be reiterated. After all levels are
   processed, all subbands appear in the initial buffer at a small offset from
   the beggining. They are read back, put in a nice matrix and given to the
   user.

   The following drawing explains how one stage of the decomposition is done.
   For clarity the buffer has been split in 8 slots. What is called a 'page'
   in the code corresponds to the length of 4 of those slots.

buffer:
     step 1      step 2      step 3         step 4

  0: empty       low         low            empty       
  1: empty       low         low            empty       } seen as the buffer
  2: empty       high        empty          low_low     } for the next stage
  3: empty   H   high    V   empty      V   low_high      
  4: image    => empty    => high_low    => high_low     
  5: image       empty       high_high      high_high     
  6: image       empty       empty          empty       
  7: image       empty       empty          empty       
*/

/*
  X[0] = X[0] + alpha * (X[-1] + X[1]); odd
  X[0] = X[0] + beta  * (X[-1] + X[1]);  even
  X[0] = X[0] + gamma * (X[-1] + X[1]); odd
  X[0] = X[0] + delta * (X[-1] + X[1]);  even
  X[0] = -scale * X[0];                 odd
  X[0] = X[0] / scale;                   even
*/
#define reset_pointers()                        \
do {						\
  pixels = buffer + page;			\
  low = buffer;					\
  high = buffer + (page >> 1);			\
  low_low =   buffer + (page >> 1) + 0 * (page >> 2); \
  low_high =  buffer + (page >> 1) + 1 * (page >> 2); \
  high_low =  buffer + (page >> 1) + 2 * (page >> 2); \
  high_high = buffer + (page >> 1) + 3 * (page >> 2); \
} while(0)

#define increment_pointers()			\
do {						\
  pixels+=2;					\
  low++;					\
  high++;					\
  high_low++;					\
  high_high++;					\
  low_high++;					\
  low_low++;					\
} while(0)

#define adjust_pointers(odd)			\
do {						\
  pixels+=odd;					\
  if(odd == 1) low++;				\
  if(odd == -1) high--;				\
} while(0)


#define hlifting(dest, source, constant, alt_source_pre, alt_source, alt_source_post, w, h, odd)	\
do {									\
  reset_pointers();							\
  for(y = 0; y < (h); y++) {						\
    (dest) = (source) + (constant) * (alt_source_pre);			\
    increment_pointers();						\
									\
    for(x = 1; x < (w)-1; x++) {					\
      (dest) = (source) + (constant) * (alt_source);			\
      increment_pointers();						\
    }									\
    									\
    (dest) = (source) + (constant) * (alt_source_post);			\
    increment_pointers();						\
    adjust_pointers(odd);						\
  }									\
} while(0)

#define hlifting(dest, source, constant, alt_source_pre, alt_source, alt_source_post, w, h, odd)	\
do {									\
  reset_pointers();							\
  for(y = 0; y < (h); y++) {						\
    (dest) = (source) + (constant) * (alt_source_pre);			\
    increment_pointers();						\
									\
    for(x = 1; x < (w)-1; x++) {					\
      (dest) = (source) + (constant) * (alt_source);			\
      increment_pointers();						\
    }									\
    									\
    (dest) = (source) + (constant) * (alt_source_post);			\
    increment_pointers();						\
    adjust_pointers(odd);						\
  }									\
} while(0)

#define vlifting(dest, source, constant, alt_source_pre, alt_source, alt_source_post, line, w, h)	\
do {									\
  reset_pointers();							\
  for(x = 0; x < w; x ++) {						\
    (dest) = (source) + (constant) * (alt_source_pre);			\
    increment_pointers();						\
  }									\
  (line) += p;								\
  for(y = 1; y < h-1; y ++) {						\
    for(x = 0; x < w; x ++) {						\
      (dest) = (source) + (constant) * (alt_source);			\
      increment_pointers();						\
    }									\
    (line) += p;							\
  }									\
  for(x = 0; x < w; x ++) {						\
    (dest) = (source) + (constant) * (alt_source_post);			\
    increment_pointers();						\
  }									\
} while(0)

#define scale(dest, scale, count)		\
do {						\
  reset_pointers();				\
  for(x = 0; x < (count); x++) {		\
    (dest)[0] *= (scale);			\
    increment_pointers();			\
  }						\
} while(0)

#define shift_up(x, l) ((x + ~(-1 << l)) >> l)
#define round_up(x, l) (shift_up(x, l) << l)

/* compute the next level decomposition using the lifting method. */
static int __wavelet2D_split (it_wavelet2D_t * wavelet)
{
  int  i, x, y, w, h, p;
  int  width;
  int  height;
  int  level;
  int  levels;
  double *buffer;
  double *pixels;
  double *low, *high, *low_high, *low_low, *high_low, *high_high;
  int  page;
  int  count;
  double const *step;
  double scale;

  assert (wavelet);

  if (wavelet->level == wavelet->levels)
    return (-IT_EINVAL);

  width = wavelet->width;
  height = wavelet->height;
  level = wavelet->level;
  levels = wavelet->levels;
  buffer = wavelet->buffer;
  count = wavelet->lifting->count / 2;
  step = wavelet->lifting->step;
  scale = wavelet->lifting->scale;

  w = shift_up (width, level);
  h = shift_up (height, level);
  page = (round_up (width, levels) * round_up (height, levels)) >> level;

  /* horizontal filtering */

  /* stage 1 : odd samples lifting */
  if (w & 1) {
    hlifting (high[0], pixels[1], step[0],
	      pixels[0] + pixels[2],
	      pixels[0] + pixels[2], pixels[0] + pixels[2], w / 2, h, 1);
  }
  else {
    hlifting (high[0], pixels[1], step[0],
	      pixels[0] + pixels[2],
	      pixels[0] + pixels[2], pixels[0] + pixels[0], w / 2, h, 0);
  }

  /* stage 2 : even samples lifting */
  if (w & 1) {
    hlifting (low[0], pixels[0], step[1],
	      high[0] + high[0],
	      high[-1] + high[0], high[-1] + high[-1], w / 2 + 1, h, -1);
  }
  else {
    hlifting (low[0], pixels[0], step[1],
	      high[0] + high[0],
	      high[-1] + high[0], high[-1] + high[0], w / 2, h, 0);
  }

  reset_pointers ();

  for (i = 1; i < count; i++) {
    /* stage 3 : odd samples lifting */
    if (w & 1) {
      hlifting (high[0], high[0], step[2 * i],
		low[0] + low[1],
		low[0] + low[1], low[0] + low[1], w / 2, h, 1);
    }
    else {
      hlifting (high[0], high[0], step[2 * i],
		low[0] + low[1],
		low[0] + low[1], low[0] + low[0], w / 2, h, 0);
    }

    /* stage 4 : even samples lifting */
    if (w & 1) {
      hlifting (low[0], low[0], step[2 * i + 1],
		high[0] + high[0],
		high[-1] + high[0], high[-1] + high[-1], w / 2 + 1, h, -1);
    }
    else {
      hlifting (low[0], low[0], step[2 * i + 1],
		high[0] + high[0],
		high[-1] + high[0], high[-1] + high[0], w / 2, h, 0);
    }
  }

  scale (low, scale, ((w + 1) / 2) * h);
  scale (high, 1.0 / scale, (w / 2) * h);

  /* vertical filtering (on horizontal high band) */
  p = w / 2;
  /* stage 1 : odd samples lifting */
  if (h & 1) {
    vlifting (high_high[0], high[p], step[0],
	      high[0] + high[2 * p],
	      high[0] + high[2 * p],
	      high[0] + high[2 * p], high, w / 2, h / 2);
  }
  else {
    vlifting (high_high[0], high[p], step[0],
	      high[0] + high[2 * p],
	      high[0] + high[2 * p], high[0] + high[0], high, w / 2, h / 2);
  }

  /* stage 2 : even samples lifting */
  if (h & 1) {
    vlifting (high_low[0], high[0], step[1],
	      high_high[0] + high_high[0],
	      high_high[-p] + high_high[0],
	      high_high[-p] + high_high[-p], high, w / 2, h / 2 + 1);
  }
  else {
    vlifting (high_low[0], high[0], step[1],
	      high_high[0] + high_high[0],
	      high_high[-p] + high_high[0],
	      high_high[-p] + high_high[0], high, w / 2, h / 2);
  }

  for (i = 1; i < count; i++) {
    /* stage 3 : odd samples lifting */
    if (h & 1) {
      vlifting (high_high[0], high_high[0], step[2 * i],
		high_low[0] + high_low[p],
		high_low[0] + high_low[p],
		high_low[0] + high_low[p], high, w / 2, h / 2);
    }
    else {
      vlifting (high_high[0], high_high[0], step[2 * i],
		high_low[0] + high_low[p],
		high_low[0] + high_low[p],
		high_low[0] + high_low[0], high, w / 2, h / 2);
    }

    /* stage 4 : even samples lifting */
    if (h & 1) {
      vlifting (high_low[0], high_low[0], step[2 * i + 1],
		high_high[0] + high_high[0],
		high_high[-p] + high_high[0],
		high_high[-p] + high_high[-p], high, w / 2, h / 2 + 1);
    }
    else {
      vlifting (high_low[0], high_low[0], step[2 * i + 1],
		high_high[0] + high_high[0],
		high_high[-p] + high_high[0],
		high_high[-p] + high_high[0], high, w / 2, h / 2);
    }
  }
  scale (high_low, scale, p * (h + 1) / 2);
  scale (high_high, 1.0 / scale, p * (h + 1) / 2);

  /* vertical filtering (on horizontal low band) */
  p = (w + 1) / 2;
  /* stage 1 : odd samples lifting */
  if (h & 1) {
    vlifting (low_high[0], low[p], step[0],
	      low[0] + low[2 * p],
	      low[0] + low[2 * p],
	      low[0] + low[2 * p], low, (w + 1) / 2, h / 2);
  }
  else {
    vlifting (low_high[0], low[p], step[0],
	      low[0] + low[2 * p],
	      low[0] + low[2 * p], low[0] + low[0], low, (w + 1) / 2, h / 2);
  }

  /* stage 2 : even samples lifting */
  if (h & 1) {
    vlifting (low_low[0], low[0], step[1],
	      low_high[0] + low_high[0],
	      low_high[-p] + low_high[0],
	      low_high[-p] + low_high[-p], low, (w + 1) / 2, h / 2 + 1);
  }
  else {
    vlifting (low_low[0], low[0], step[1],
	      low_high[0] + low_high[0],
	      low_high[-p] + low_high[0],
	      low_high[-p] + low_high[0], low, (w + 1) / 2, h / 2);
  }

  for (i = 1; i < count; i++) {
    /* stage 3 : odd samples lifting */
    if (h & 1) {
      vlifting (low_high[0], low_high[0], step[2 * i],
		low_low[0] + low_low[p],
		low_low[0] + low_low[p],
		low_low[0] + low_low[p], low, (w + 1) / 2, h / 2);
    }
    else {
      vlifting (low_high[0], low_high[0], step[2 * i],
		low_low[0] + low_low[p],
		low_low[0] + low_low[p],
		low_low[0] + low_low[0], low, (w + 1) / 2, h / 2);
    }

    /* stage 4 : even samples lifting */
    if (h & 1) {
      vlifting (low_low[0], low_low[0], step[2 * i + 1],
		low_high[0] + low_high[0],
		low_high[-p] + low_high[0],
		low_high[-p] + low_high[-p], low, (w + 1) / 2, h / 2 + 1);
    }
    else {
      vlifting (low_low[0], low_low[0], step[2 * i + 1],
		low_high[0] + low_high[0],
		low_high[-p] + low_high[0],
		low_high[-p] + low_high[0], low, (w + 1) / 2, h / 2);
    }
  }

  scale (low_low, scale, p * (h + 1) / 2);
  scale (low_high, 1.0 / scale, p * (h + 1) / 2);

  wavelet->level++;

  return (0);
}

/* reconstruct the previous level of decomposition
   by inverting the lifting steps. */
static int __wavelet2D_merge (it_wavelet2D_t * wavelet)
{
  int  i, x, y, w, h, p;
  int  width;
  int  height;
  int  level;
  int  levels;
  double *pixels;
  double *buffer;
  double *low, *high, *low_high, *low_low, *high_low, *high_high;
  int  page;
  int  count;
  double const *step;
  double scale;

  assert (wavelet);

  if (wavelet->level == 0)
    return (-IT_EINVAL);

  wavelet->level--;

  width = wavelet->width;
  height = wavelet->height;
  level = wavelet->level;
  levels = wavelet->levels;
  buffer = wavelet->buffer;
  count = wavelet->lifting->count / 2;
  step = wavelet->lifting->step;
  scale = wavelet->lifting->scale;

  w = shift_up (width, level);
  h = shift_up (height, level);
  page = (round_up (width, levels) * round_up (height, levels)) >> level;

  /* vertical filtering (on horizontal low band) */
  p = (w + 1) / 2;

  scale (low_low, 1.0 / scale, p * (h + 1) / 2);
  scale (low_high, scale, p * (h + 1) / 2);

  for (i = 1; i < count; i++) {
    /* stage 4 : even samples lifting */
    if (h & 1) {
      vlifting (low_low[0], low_low[0], -step[2 * i + 1],
		low_high[0] + low_high[0],
		low_high[-p] + low_high[0],
		low_high[-p] + low_high[-p], low, (w + 1) / 2, h / 2 + 1);
    }
    else {
      vlifting (low_low[0], low_low[0], -step[2 * i + 1],
		low_high[0] + low_high[0],
		low_high[-p] + low_high[0],
		low_high[-p] + low_high[0], low, (w + 1) / 2, h / 2);
    }

    /* stage 3 : odd samples lifting */
    if (h & 1) {
      vlifting (low_high[0], low_high[0], -step[2 * i],
		low_low[0] + low_low[p],
		low_low[0] + low_low[p],
		low_low[0] + low_low[p], low, (w + 1) / 2, h / 2);
    }
    else {
      vlifting (low_high[0], low_high[0], -step[2 * i],
		low_low[0] + low_low[p],
		low_low[0] + low_low[p],
		low_low[0] + low_low[0], low, (w + 1) / 2, h / 2);
    }
  }

  /* stage 2 : even samples lifting */
  if (h & 1) {
    vlifting (low[0], low_low[0], -step[1],
	      low_high[0] + low_high[0],
	      low_high[-p] + low_high[0],
	      low_high[-p] + low_high[-p], low, (w + 1) / 2, h / 2 + 1);
  }
  else {
    vlifting (low[0], low_low[0], -step[1],
	      low_high[0] + low_high[0],
	      low_high[-p] + low_high[0],
	      low_high[-p] + low_high[0], low, (w + 1) / 2, h / 2);
  }

  /* stage 1 : odd samples lifting */
  if (h & 1) {
    vlifting (low[p], low_high[0], -step[0],
	      low[0] + low[2 * p],
	      low[0] + low[2 * p],
	      low[0] + low[2 * p], low, (w + 1) / 2, h / 2);
  }
  else {
    vlifting (low[p], low_high[0], -step[0],
	      low[0] + low[2 * p],
	      low[0] + low[2 * p], low[0] + low[0], low, (w + 1) / 2, h / 2);
  }

  /* vertical filtering (on horizontal high band) */
  p = w / 2;

  scale (high_low, 1.0 / scale, p * (h + 1) / 2);
  scale (high_high, scale, p * (h + 1) / 2);

  for (i = 1; i < count; i++) {
    /* stage 4 : even samples lifting */
    if (h & 1) {
      vlifting (high_low[0], high_low[0], -step[2 * i + 1],
		high_high[0] + high_high[0],
		high_high[-p] + high_high[0],
		high_high[-p] + high_high[-p], high, w / 2, h / 2 + 1);
    }
    else {
      vlifting (high_low[0], high_low[0], -step[2 * i + 1],
		high_high[0] + high_high[0],
		high_high[-p] + high_high[0],
		high_high[-p] + high_high[0], high, w / 2, h / 2);
    }

    /* stage 3 : odd samples lifting */
    if (h & 1) {
      vlifting (high_high[0], high_high[0], -step[2 * i],
		high_low[0] + high_low[p],
		high_low[0] + high_low[p],
		high_low[0] + high_low[p], high, w / 2, h / 2);
    }
    else {
      vlifting (high_high[0], high_high[0], -step[2 * i],
		high_low[0] + high_low[p],
		high_low[0] + high_low[p],
		high_low[0] + high_low[0], high, w / 2, h / 2);
    }
  }

  /* stage 2 : even samples lifting */
  if (h & 1) {
    vlifting (high[0], high_low[0], -step[1],
	      high_high[0] + high_high[0],
	      high_high[-p] + high_high[0],
	      high_high[-p] + high_high[-p], high, w / 2, h / 2 + 1);
  }
  else {
    vlifting (high[0], high_low[0], -step[1],
	      high_high[0] + high_high[0],
	      high_high[-p] + high_high[0],
	      high_high[-p] + high_high[0], high, w / 2, h / 2);
  }

  /* stage 1 : odd samples lifting */
  if (h & 1) {
    vlifting (high[p], high_high[0], -step[0],
	      high[0] + high[2 * p],
	      high[0] + high[2 * p],
	      high[0] + high[2 * p], high, w / 2, h / 2);
  }
  else {
    vlifting (high[p], high_high[0], -step[0],
	      high[0] + high[2 * p],
	      high[0] + high[2 * p], high[0] + high[0], high, w / 2, h / 2);
  }

  /* horizontal filtering */
  scale (low, 1.0 / scale, ((w + 1) / 2) * h);
  scale (high, scale, (w / 2) * h);

  for (i = 1; i < count; i++) {
    /* stage 4 : even samples lifting */
    if (w & 1) {
      hlifting (low[0], low[0], -step[2 * i + 1],
		high[0] + high[0],
		high[-1] + high[0], high[-1] + high[-1], w / 2 + 1, h, -1);
    }
    else {
      hlifting (low[0], low[0], -step[2 * i + 1],
		high[0] + high[0],
		high[-1] + high[0], high[-1] + high[0], w / 2, h, 0);
    }

    /* stage 3 : odd samples lifting */
    if (w & 1) {
      hlifting (high[0], high[0], -step[2 * i],
		low[0] + low[1],
		low[0] + low[1], low[0] + low[1], w / 2, h, 1);
    }
    else {
      hlifting (high[0], high[0], -step[2 * i],
		low[0] + low[1],
		low[0] + low[1], low[0] + low[0], w / 2, h, 0);
    }
  }

  /* stage 2 : even samples lifting */
  if (w & 1) {
    hlifting (pixels[0], low[0], -step[1],
	      high[0] + high[0],
	      high[-1] + high[0], high[-1] + high[-1], w / 2 + 1, h, -1);
  }
  else {
    hlifting (pixels[0], low[0], -step[1],
	      high[0] + high[0],
	      high[-1] + high[0], high[-1] + high[0], w / 2, h, 0);
  }

  /* stage 1 : odd samples lifting */
  if (w & 1) {
    hlifting (pixels[1], high[0], -step[0],
	      pixels[0] + pixels[2],
	      pixels[0] + pixels[2], pixels[0] + pixels[2], w / 2, h, 1);
  }
  else {
    hlifting (pixels[1], high[0], -step[0],
	      pixels[0] + pixels[2],
	      pixels[0] + pixels[2], pixels[0] + pixels[0], w / 2, h, 0);
  }

  return (0);
}

/* arrange the transformed coefficients in an image */
/* where the top left corner contains the lowest band.  */
static int __wavelet_flatten (it_wavelet2D_t * it_this, mat image)
{
  int  x, y, wf, hf, wc, hc, ps;
  int  width, height, level, page;
  int  levels;
  double *band;
  double *buffer = it_this->buffer;

  width = it_this->width;
  height = it_this->height;
  levels = it_this->levels;
  page = round_up (width, levels) * round_up (height, levels);
  buffer = it_this->buffer + page;
  ps = width;

  if (width != mat_width (image))
    return (-IT_EINVAL);
  if (height != mat_height (image))
    return (-IT_EINVAL);

  hc = height;
  wc = width;

  for (level = 0; level < levels; level++) {
    page = (round_up (width, levels) * round_up (height, levels)) >> level;
    buffer = (it_this->buffer) + (page >> 1);

    ps = (ps + 1) / 2;		/* upper rounded */
    hf = hc / 2;
    hc = (hc + 1) / 2;
    wf = wc / 2;
    wc = (wc + 1) / 2;

    band = buffer + 2 * (page >> 2);	/* HL */
    for (y = 0; y < hc; y++)
      for (x = 0; x < wf; x++)
	image[y][x + wc] = band[y * wf + x];

    band = buffer + 1 * (page >> 2);	/* LH */
    for (y = 0; y < hf; y++)
      for (x = 0; x < wc; x++)
	image[y + hc][x] = band[y * wc + x];

    band = buffer + 3 * (page >> 2);	/* HH */
    for (y = 0; y < hf; y++)
      for (x = 0; x < wf; x++)
	image[y + hc][x + wc] = band[y * wf + x];
  }

  band = buffer;		/* LL */
  for (y = 0; y < hc; y++)
    for (x = 0; x < wc; x++)
      image[y][x] = band[y * wc + x];

  return (0);
}

/* read the transformed coefficients from an image */
/* where the top left corner contains the lowest band.  */
static int __wavelet_unflatten (it_wavelet2D_t * it_this, mat image)
{
  int  x, y, wf, hf, wc, hc, ps;
  int  width, height, level, page;
  int  levels;
  double *band;
  double *buffer = it_this->buffer;

  width = it_this->width;
  height = it_this->height;
  levels = it_this->levels;
  page = round_up (width, levels) * round_up (height, levels);
  buffer = it_this->buffer + page;
  ps = width;

  if (width != mat_width (image))
    return (-IT_EINVAL);
  if (height != mat_height (image))
    return (-IT_EINVAL);

  /* we are given a fully decomposed image */
  it_this->level = levels;
  hc = height;
  wc = width;

  for (level = 0; level < levels; level++) {
    page = (round_up (width, levels) * round_up (height, levels)) >> level;
    buffer = (it_this->buffer) + (page >> 1);

    ps = (ps + 1) / 2;		/* upper rounded */
    hf = hc / 2;
    hc = (hc + 1) / 2;
    wf = wc / 2;
    wc = (wc + 1) / 2;

    band = buffer + 2 * (page >> 2);	/* HL */
    for (y = 0; y < hc; y++)
      for (x = 0; x < wf; x++)
	band[y * wf + x] = image[y][x + wc];

    band = buffer + 1 * (page >> 2);	/* LH */
    for (y = 0; y < hf; y++)
      for (x = 0; x < wc; x++)
	band[y * wc + x] = image[y + hc][x];

    band = buffer + 3 * (page >> 2);	/* HH */
    for (y = 0; y < hf; y++)
      for (x = 0; x < wf; x++)
	band[y * wf + x] = image[y + hc][x + wc];
  }

  band = buffer;		/* LL */
  for (y = 0; y < hc; y++)
    for (x = 0; x < wc; x++)
      band[y * wc + x] = image[y][x];

  return (0);
}

/* compute the wavelet transform of the image */
static Mat __it_wavelet2D_transform (it_transform2D_t * transform,
				     Mat __image)
{
  it_wavelet2D_t *wavelet = IT_WAVELET2D (transform);
  int  x, y;
  int  width, height;
  int  page;
  int  levels;
  double *buffer;
  mat  flat;
  mat  image = (mat) __image;
  int  free_on_exit = 0;

  assert (image);

  /* check the input is actually a mat */
  assert (Vec_header (image[0]).element_size == sizeof (double));

  width = mat_width (image);
  height = mat_height (image);

  if (!wavelet->width || !wavelet->height) {
    it_transform2D_set_size (wavelet, width, height);
    free_on_exit = 1;
  }

  flat = mat_new (height, width);
  levels = wavelet->levels;
  page = round_up (width, levels) * round_up (height, levels);
  buffer = wavelet->buffer + page;

  for (y = 0; y < height; y++)
    for (x = 0; x < width; x++)
      buffer[y * width + x] = image[y][x];

  while (wavelet->level < levels)
    __wavelet2D_split (wavelet);

  __wavelet_flatten (IT_WAVELET2D (transform), flat);

  if (free_on_exit)
    it_transform2D_clear_size (wavelet);

  return ((Mat) flat);
}

/* compute the inverse wavelet transform of the coefficients */
static Mat __it_wavelet2D_itransform (it_transform2D_t * transform,
				      Mat __flat)
{
  it_wavelet2D_t *wavelet = IT_WAVELET2D (transform);
  int  width, height;
  int  x, y;
  int  page;
  int  levels;
  double *buffer;
  mat  image;
  mat  flat = (mat) __flat;
  int  free_on_exit = 0;

  assert (flat);

  /* check the input is actually a vec */
  assert (Vec_header (flat[0]).element_size == sizeof (double));

  width = mat_width (flat);
  height = mat_height (flat);

  if (!wavelet->width || !wavelet->height) {
    it_transform2D_set_size (wavelet, width, height);
    free_on_exit = 1;
  }

  image = mat_new (height, width);
  levels = wavelet->levels;
  page = round_up (width, levels) * round_up (height, levels);
  buffer = wavelet->buffer + page;

  __wavelet_unflatten (IT_WAVELET2D (transform), flat);

  while (wavelet->level)
    __wavelet2D_merge (wavelet);

  for (y = 0; y < height; y++)
    for (x = 0; x < width; x++)
      image[y][x] = buffer[y * width + x];

  if (free_on_exit)
    it_transform2D_clear_size (wavelet);

  return ((Mat) image);
}

static int __it_wavelet2D_copy (it_wavelet2D_t * it_this,
				it_wavelet2D_t * source)
{
  assert (it_this->levels == source->levels);

  if (source->width && source->height) {
    it_transform2D_set_size (it_this, source->width, source->height);
    vec_copy (it_this->buffer, source->buffer);
  }

  return (0);
}

static void __wavelet2D_get_output_size (it_transform2D_t * transform,
					 idx_t * width, idx_t * height)
{
  /* this transform is critically sampled */
}

static void __wavelet2D_set_size (it_transform2D_t * transform,
				  idx_t width, idx_t height)
{
  it_wavelet2D_t *wavelet = IT_WAVELET2D (transform);
  int  levels = wavelet->levels;

  if (wavelet->width && wavelet->height)
    vec_delete (wavelet->buffer);

  if (width && height)
    wavelet->buffer =
      vec_new (2 * round_up (width, levels) * round_up (height, levels));

  wavelet->width = width;
  wavelet->height = height;
}

static void __wavelet2D_get_size (it_transform2D_t * transform,
				  idx_t * width, idx_t * height)
{
  it_wavelet2D_t *wavelet = IT_WAVELET2D (transform);

  *width = wavelet->width;
  *height = wavelet->height;
}

static void __wavelet2D_destructor (it_object_t * it_this)
{
  it_wavelet2D_t *wavelet = IT_WAVELET2D (it_this);

  if (wavelet->width && wavelet->height)
    vec_delete (wavelet->buffer);

  /* call the parent destructor */
  wavelet->it_overloaded (destructor) (it_this);
}

it_instanciate (it_wavelet2D_t)
{
  it_new_args_start ();
  it_construct (it_transform2D_t);
  it_set_magic (it_this, it_wavelet2D_t);

  /* overload the virtual destructor */
  it_overload (it_this, it_object_t, destructor, __wavelet2D_destructor);

  IT_TRANSFORM2D (it_this)->transform = __it_wavelet2D_transform;
  IT_TRANSFORM2D (it_this)->itransform = __it_wavelet2D_itransform;
  IT_TRANSFORM2D (it_this)->get_output_size = __wavelet2D_get_output_size;
  IT_TRANSFORM2D (it_this)->set_size = __wavelet2D_set_size;
  IT_TRANSFORM2D (it_this)->get_size = __wavelet2D_get_size;

  it_this->copy = __it_wavelet2D_copy;

  it_this->level = 0;
  it_this->lifting = it_new_args_next (it_wavelet_lifting_t *);
  it_this->levels = it_new_args_next (int);
  it_this->width = 0;
  it_this->height = 0;

  it_new_args_stop ();

  return (it_this);
}

/*--------------------------------------------------------------------*/
mat *it_wavelet2D_split (mat wav, int nb_levels)
{
  int  mid_row, mid_col, nb_row, nb_col, l;
  mat *subbands;

  assert (wav);

  nb_row = mat_height (wav);
  nb_col = mat_width (wav);
  mid_row = (nb_row + 1) / 2;
  mid_col = (nb_col + 1) / 2;

  subbands = (mat *) malloc (sizeof (mat) * (3 * nb_levels + 1));

  for (l = nb_levels; l > 0; l--) {
    subbands[l * 3 - 2] =
      mat_get_submatrix (wav, 0, mid_col, mid_row - 1, nb_col - 1);
    subbands[l * 3 - 1] =
      mat_get_submatrix (wav, mid_row, 0, nb_row - 1, mid_col - 1);
    subbands[l * 3] =
      mat_get_submatrix (wav, mid_row, mid_col, nb_row - 1, nb_col - 1);

    nb_row = mid_row;
    nb_col = mid_col;
    mid_row = (nb_row + 1) / 2;
    mid_col = (nb_col + 1) / 2;
  }

  subbands[0] = mat_get_submatrix (wav, 0, 0, nb_row -1 , nb_col -1);
  return subbands;
}

mat it_wavelet2D_merge (mat * subbands, int nb_levels)
{
  int  mid_row, mid_col, nb_row, nb_col, l;
  mat  wav;

  assert (subbands);

  wav = mat_new (mat_height (subbands[nb_levels * 3])
	  + mat_height (subbands[nb_levels * 3 - 2]),
	  mat_width (subbands[nb_levels * 3])
	  + mat_width (subbands[nb_levels * 3 - 1]));

  nb_row = mat_height (wav);
  nb_col = mat_width (wav);
  mid_row = (nb_row + 1) / 2;
  mid_col = (nb_col + 1) / 2;

  for (l = nb_levels; l > 0; l--) {
    mat_set_submatrix (wav, subbands[l * 3 - 2], 0, mid_col);
    mat_set_submatrix (wav, subbands[l * 3 - 1], mid_row, 0);
    mat_set_submatrix (wav, subbands[l * 3], mid_row, mid_col);

    nb_row = mid_row;
    nb_col = mid_col;
    mid_row = (nb_row + 1) / 2;
    mid_col = (nb_col + 1) / 2;
  }

  mat_set_submatrix (wav, subbands[0], 0, 0);
  return wav;
}

mat it_dwt2D (mat m, it_wavelet_lifting_t const *lifting, int levels)
{
  it_wavelet2D_t *wavelet;
  mat  t;

  wavelet = it_wavelet2D_new (lifting, levels);

  t = (mat) it_transform2D (wavelet, m);

  it_delete (wavelet);

  return (t);
}

mat it_idwt2D (mat t, it_wavelet_lifting_t const *lifting, int levels)
{
  it_wavelet2D_t *wavelet;
  mat  m;

  wavelet = it_wavelet2D_new (lifting, levels);

  m = (mat) it_itransform2D (wavelet, t);

  it_delete (wavelet);

  return (m);
}
