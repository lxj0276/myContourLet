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
  Filtering
  Copyright (C) 2005 Vivien Chappelier.
*/

#include <it/types.h>
#include <it/filter.h>

/* X = x + d with proper bound checking and mirroring */
#define mirror_add(X, x, d, min, max) do {	\
    (X) = (x) + (d);				\
    if((X) < (min) || (X) >= (max))		\
      (X) = (x) - (d);				\
} while(0)

mat mat_filter_fir (mat input, mat filter, int px, int py)
{
  idx_t X, Y;
  idx_t x, y;
  idx_t dx, dy;
  idx_t wf, hf;
  idx_t w, h;
  mat  output;

  w = mat_width (input);
  h = mat_height (input);
  wf = mat_width (filter);
  hf = mat_height (filter);

  assert (px >= 0 && px < wf);
  assert (py >= 0 && py < hf);

  output = mat_new_zeros (h, w);

  for (x = 0; x < w; x++) {
    for (y = 0; y < h; y++) {

      /* filter around position x,y */
      for (dx = -px; dx < wf - px; dx++) {
	for (dy = -py; dy < hf - py; dy++) {
	  /* compute the input sample location */
	  mirror_add (X, x, dx, 0, w);
	  mirror_add (Y, y, dy, 0, h);

	  output[y][x] += filter[dy + py][dx + px] * input[Y][X];
	}
      }

    }
  }

  return (output);
}
