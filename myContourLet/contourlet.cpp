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
#include "dfb.h"
#include "contourlet.h"
#include "pyramid.h"

//#define MPEG
//#define VILLA_13_11
#define ANTON_9_7

contourlet_t *contourlet_new(int ct_levels, ivec dfb_levels)
{
  contourlet_t *ct;
  int l;

  ct = (contourlet_t *) malloc(sizeof(contourlet_t));
  ct->ct_levels = ct_levels;
  ct->dfb_levels = ivec_clone(dfb_levels);
  ct->high = (double ****)calloc(ct_levels, sizeof(mat *));
  ct->low = NULL;
  /* pyramid filters */
  ct->H0 = 0;
  ct->G0 = 0;
#ifdef MPEG
  ct->H0 = vec_new_string("2,  0, -4, -3,  5, 19,   26,   19,    5, -3, -4,  0,  2");
  vec_div_by(ct->H0, 64/sqrt(2));
  ct->G0 = vec_new_string("0,  1,  0, -5,  0, 20,   32,   20,    0, -5,  0,  1,  0");
  vec_div_by(ct->G0, 64/sqrt(2));
#endif
#ifdef ANTON_9_7
  ct->H0 = vec_new_string(" 3.782845550699535e-02, -2.384946501937986e-02, -1.106244044184226e-01, 3.774028556126536e-01, 8.526986790094022e-01, 3.774028556126537e-01, -1.106244044184226e-01, -2.384946501937986e-02, 3.782845550699535e-02");
  ct->G0 = vec_new_string("-6.453888262893856e-02, -4.068941760955867e-02, 4.180922732222124e-01, 7.884856164056651e-01, 4.180922732222124e-01, -4.068941760955867e-02, -6.453888262893856e-02");
#endif
#ifdef VILLA_13_11
  ct->H0 = vec_new_string("-8.472827741318157e-03, 3.759210316686883e-03, 4.728175282882753e-02, -3.347508104780150e-02, -6.887811419061032e-02, 3.832692613243884e-01, 7.672451593927493e-01, 3.832692613243889e-01, -6.887811419061045e-02, -3.347508104780156e-02, 4.728175282882753e-02, 3.759210316686883e-03, -8.472827741318157e-03");
  ct->G0 = vec_new_string(" 1.418215589126359e-02, 6.292315666859828e-03, -1.087373652243805e-01, -6.916271012030040e-02, 4.481085999263908e-01, 8.328475700934288e-01, 4.481085999263908e-01, -6.916271012030040e-02, -1.087373652243805e-01, 6.292315666859828e-03, 1.418215589126359e-02");
#endif
  assert(ct->H0);
  assert(ct->G0);

  /* allocate directionnal subbands */
  for(l = 0; l < ct_levels; l++) {
    /* allocate directionnal subbands for that level */
    ct->high[l] = (mat *) calloc(1 << ct->dfb_levels[l], sizeof(mat));
  }

  return(ct);
}

void contourlet_delete(contourlet_t *ct)
{
  int l;
  int ct_levels = ct->ct_levels;

  for(l = 0; l < ct_levels; l++)
    free(ct->high[l]);

  vec_delete(ct->H0);
  vec_delete(ct->G0);

  ivec_delete(ct->dfb_levels);
  free(ct->high);
}

/* contourlet analysis */
void contourlet_transform(contourlet_t *ct, /* contourlet decomposition */
			  mat image)        /* input image              */
{
  int w, h;
  int wh, hh;
  int wl, hl;
  int l;
  mat input, high, low;
  int ct_levels = ct->ct_levels;
  int dfb_levels;
  vec H0 = ct->H0;
  vec G0 = ct->G0;
  dfb_t *dfb;

  w = mat_width(image);
  h = mat_height(image);

  low = image;

  /* multiresolution analysis */
  for(l = 0; l < ct_levels; l++) {

    input = low;

    /* compute width and height of the current level subbands */
    wl = ct_shift_up(w, l+1);
    hl = ct_shift_up(h, l+1);
    wh = ct_shift_up(w, l);
    hh = ct_shift_up(h, l);
    dfb_levels = ct->dfb_levels[l];

    /* allocate low and high subbands */ 
    low = mat_new(hl, wl);
    high = mat_new(hh, wh);

    /* split with the laplacian pyramid */
    laplacian_pyramid_split(input, low, high, H0, G0);//拉普拉斯金字塔分离0422lxj

    // TEMP
    //    char buffer[256];
    //    sprintf(buffer, "low%d.pgm", l);
    //    mat_pgm_write(buffer, low);
    //    mat_incr(high, 128);
    //    sprintf(buffer, "high%d.pgm", l);
    //    mat_pgm_write(buffer, high);
    //    mat_decr(high, 128);

    /* directionnal analysis */
    if(dfb_levels) 
	{
      /* create the directional filter bank */
      dfb = dfb_new(wh, hh, dfb_levels);

      /* directional decomposition */
      dfb_transform(dfb, ct->high[l], high, dfb_levels);

      /* release resources */
      dfb_delete(dfb);
      mat_delete(high);
    } 
	else
	{
      ct->high[l][0] = high;
    }

    if(l) mat_delete(input);
  }

  ct->low = low;
}

/* contourlet synthesis */
void contourlet_itransform(contourlet_t *ct, /* contourlet decomposition */
			   mat image)        /* output image             */
{
  int w, h;
  int wh, hh;
  int wl, hl;
  int l;
  mat output, high, low;
  int ct_levels = ct->ct_levels;
  int dfb_levels;
  vec H0 = ct->H0;
  vec G0 = ct->G0;
  dfb_t *dfb;

  w = mat_width(image);
  h = mat_height(image);

  low = ct->low;

  /* multiresolution analysis */
  for(l = ct_levels - 1; l >= 0; l--) {

    /* compute width and height of the current level subbands */
    wl = ct_shift_up(w, l+1);
    hl = ct_shift_up(h, l+1);
    wh = ct_shift_up(w, l);
    hh = ct_shift_up(h, l);
    dfb_levels = ct->dfb_levels[l];

    /* allocate low and high subbands */ 
    output = mat_new(hh, wh);
    high = mat_new(hh, wh);

    /* directionnal synthesis */
    if(dfb_levels) {
      /* create the directional filter bank */
      dfb = dfb_new(wh, hh, dfb_levels);

      int i;
      // TEMP
      for(i = 0; i < (1 << dfb_levels); i++) {
	char buffer[256];
	sprintf(buffer, "rband%d_%d.pgm", l, i);
	mat_pgm_write(buffer, ct->high[l][i]);
      }

      /* directional decomposition */
      dfb_itransform(dfb, ct->high[l], high, dfb_levels);

      // TEMP
          char buffer[256];
         sprintf(buffer, "rlow%d.pgm", l);
         mat_pgm_write(buffer, low);
         mat_incr(high, 128);
         sprintf(buffer, "rhigh%d.pgm", l);
         mat_pgm_write(buffer, high);
         mat_decr(high, 128);
      dfb_delete(dfb);
    } else {
      mat_copy(high, ct->high[l][0]);
    }

    /* merge the laplacian pyramid */
    laplacian_pyramid_merge(output, low, high, H0, G0);

    /* release resources */
    mat_delete(high);
    mat_delete(low);

    low = output;
  }

  mat_copy(image, low);
}
