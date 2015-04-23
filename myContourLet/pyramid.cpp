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
#include "it/mat.h"
#include "it/io.h"
#include "it/distance.h"
#include <stdio.h>
#include <math.h>

/* X = x + d with proper bound checking and mirroring */
#define mirror_add(X, x, d, min, max) do {	\
    (X) = (x) + (d);				\
    if((X) < (min) || (X) >= (max))		\
      (X) = (x) - (d);				\
} while(0)

/* split an image into two laplacian pyramid bands */
void laplacian_pyramid_split(mat image, /* input band                 */
			     mat low,   /* low-frequency band         */
			     mat high,  /* high-frequency band        */
			     vec H0,    /* symmetric analysis  filter */
			     vec G0)    /* symmetric synthesis filter */
{
  int x, y;
  int w, h;
  idx_t dx, dy;
  idx_t wf, hf;
  idx_t wl, hl;
  int px, py;
  int X, Y;

  w = mat_width(image);
  h = mat_height(image);
  wl = (w + 1) / 2;
  hl = (h + 1) / 2;

  assert(mat_width(high)  == w);
  assert(mat_height(high) == h);
  assert(mat_width(low)  == wl);
  assert(mat_height(low) == hl);

  wf = hf = vec_length(H0);
  px = wf / 2;
  py = hf / 2;

  for(x = 0; x < wl; x++) {
    for(y = 0; y < hl; y++) {
      
      /* filter around position x,y */
      low[y][x] = 0;
      for(dx = -px; dx < wf-px; dx++) {
	for(dy = -py; dy < hf-py; dy++) {
	  /* compute the input sample location */
	  mirror_add(X, 2*x, dx, 0, w);
	  mirror_add(Y, 2*y, dy, 0, h);

	  low[y][x] += H0[dy+py] * H0[dx+px] * image[Y][X];
	}
      }

    }
  }

  wf = hf = vec_length(G0);
  px = wf / 2;
  py = hf / 2;

  for(x = 0; x < w; x++) {
    for(y = 0; y < h; y++) {
      
      /* filter around position x,y */
      high[y][x] = image[y][x];
      for(dx = -px; dx < wf-px; dx++) {
	for(dy = -py; dy < hf-py; dy++) {
	  /* compute the input sample location */
	  mirror_add(X, x, dx, 0, w);
	  mirror_add(Y, y, dy, 0, h);

	  if(!(X & 1) && !(Y & 1))
	    high[y][x] -= G0[dy+py] * G0[dx+px] * low[Y/2][X/2];
	}
      }

    }
  }
}


/* merge two laplacian pyramid bands into an image */
void laplacian_pyramid_merge(mat image, /* output band                */
			     mat low,   /* low-frequency band         */
			     mat high,  /* high-frequency band        */
			     vec H0,    /* symmetric analysis  filter */
			     vec G0)    /* symmetric synthesis filter */
{
  int x, y;
  int w, h;
  idx_t dx, dy;
  idx_t wf, hf;
  idx_t wl, hl;
  int px, py;
  int X, Y;

  w = mat_width(image);
  h = mat_height(image);
  wl = (w + 1) / 2;
  hl = (h + 1) / 2;

  assert(mat_width(high)  == w);
  assert(mat_height(high) == h);
  assert(mat_width(low)  == wl);
  assert(mat_height(low) == hl);

/* use pseudo inverse reconstruction */
/* this assumes the filters are orthogonal */
/* the 9/7 are close enough to orthogonality for this to work quite well */
#define DUAL
#ifdef DUAL
  wf = hf = vec_length(H0);
  px = wf / 2;
  py = hf / 2;

  for(x = 0; x < wl; x++) {
    for(y = 0; y < hl; y++) {
      
      /* filter around position x,y */
      for(dx = -px; dx < wf-px; dx++) {
	for(dy = -py; dy < hf-py; dy++) {
	  /* compute the input sample location */
	  mirror_add(X, 2*x, dx, 0, w);
	  mirror_add(Y, 2*y, dy, 0, h);

	  low[y][x] -= H0[dy+py] * H0[dx+px] * high[Y][X];
	}
      }

    }
  }
#endif

  wf = hf = vec_length(G0);
  px = wf / 2;
  py = hf / 2;

  for(x = 0; x < w; x++) {
    for(y = 0; y < h; y++) {
      
      /* filter around position x,y */
      image[y][x] = high[y][x];
      for(dx = -px; dx < wf-px; dx++) {
	for(dy = -py; dy < hf-py; dy++) {
	  /* compute the input sample location */
	  mirror_add(X, x, dx, 0, w);
	  mirror_add(Y, y, dy, 0, h);

	  if(!(X & 1) && !(Y & 1))
	    image[y][x] += G0[dy+py] * G0[dx+px] * low[Y/2][X/2];
	}
      }

    }
  }
}

