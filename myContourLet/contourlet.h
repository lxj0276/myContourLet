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

#ifndef __CONTOURLET_H
#define __CONTOURLET_H

#include "it/mat.h"

typedef struct {
  int ct_levels;   /* number of pyramid decomposition levels     */
  int wt_levels;   /* number of wavelet decomposition levels     */
  ivec dfb_levels; /* number of directional decomposition levels */
  mat low;         /* low-frequency band                         */
  mat **high;      /* high-frequency bands [level][orientation]  */
  mat *dwt;        /* wavelet subband decomposition of 'low' */ 
  vec H0;          /* pyramid analysis filter */
  vec G0;          /* pyramid synthesis filter */
} contourlet_t;

/* shift and round up */
#define ct_shift_up(x, l) (((x) + (1 << (l)) - 1) >> (l))

contourlet_t *contourlet_new(int levels, ivec dfb_levels);
void contourlet_delete(contourlet_t *ct);

void contourlet_transform(contourlet_t *ct, mat image);
void contourlet_itransform(contourlet_t *ct, mat image);   
#endif
