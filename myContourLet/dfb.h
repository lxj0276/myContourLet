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

#ifndef __DFB_H
#define __DFB_H

#include "it/mat.h"

typedef struct _slim_dfb_subband_ {
  int Q[2][2];                   /* downsampling matrix */
  int R[2][2];                   /* rotation matrix     */
  double norm;                   /* filter norm         */
  mat image;                     /* coefficients        */
  bmat mask;                     /* valid positions     */
} dfb_subband_t;

typedef struct _slim_dfb_ {
  int width;
  int height;
  int levels;
  ivec ox;                       /* subband x origin */
  ivec oy;                       /* subband y origin */
  ivec rx;                       /* subband x end    */
  ivec ry;                       /* subband y end    */
  dfb_subband_t **bands;
  double *filter_norm;
} dfb_t;

#define dfb_band(dfb, dfb_level, index) \
  dfb->bands[(1 << (dfb_level)) + index]

void direct_poly_filter_downsample(dfb_subband_t *hband,
				   dfb_subband_t *vband,
				   dfb_subband_t *input,
				   int *Q, int *R,
				   int w, int h);

void direct_poly_interpolate_upsample(dfb_subband_t *hband,
				      dfb_subband_t *vband,
				      dfb_subband_t *input,
				      int *Q, int *R,
				      int w, int h);

void resample_and_normalize(dfb_subband_t *input,
			    mat tmp,
			    int *B,
			    double norm,
			    int w, int h);

void dfb_init();

dfb_subband_t *dfb_subband_new();
void dfb_subband_delete(dfb_subband_t *subband);

dfb_t *dfb_new(int width, int height, int levels);
void dfb_delete(dfb_t *dfb);

void dfb_transform(dfb_t *dfb, mat *bands, mat image, int levels);
void dfb_itransform(dfb_t *dfb, mat *bands, mat image, int levels);

void dfb_flatten(mat *high, mat image, int levels);

#endif
