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
  Channel 
  Copyright (C) 2005 Vivien Chappelier, Herve Jegou
*/

#ifndef __it_channel_h
#define __it_channel_h

#include <it/vec.h>

#ifdef __cplusplus
extern "C" {
#endif


/* Basic BPSK modulation */
vec modulate_bpsk( bvec b );

/* Binary symmetrical channel */
bvec channel_bsc( bvec v, double crossover_proba );

/* Return a vector pertubated by a white gaussian additive noise of variance sigma */
vec channel_awgn( vec v, double sigma );



#ifdef __cplusplus
}
#endif

#endif

