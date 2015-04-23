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
  1D wavelet transform using lifting.
  Copyright (C) 2005 Vivien Chappelier, Herve Jegou
*/

#ifndef __it_wavelet_h
#define __it_wavelet_h

#include "it/types.h"
#include "it/mat.h"
#include "it/transform.h"

#ifdef __cplusplus
extern "C" {
#endif

/* the lifting structure */
typedef struct _wavelet_lifting_ {
  int const count;       /* number of lifting steps, must be even and > 0 */
  double const scale;    /* final scaling coefficient */
  double const *step;   /* steps in predict/update order */
} it_wavelet_lifting_t;

/* some predefined filters */
extern it_wavelet_lifting_t const *it_wavelet_lifting_97;
extern it_wavelet_lifting_t const *it_wavelet_lifting_53;

typedef struct _it_wavelet_ {
  it_extends(it_transform_t);

  void (* it_overloaded(destructor))(it_object_t *it_this);

  it_wavelet_lifting_t *lifting;
  size_t length;                 /* length of the vector to decompose */
  int level;                     /* current decomposition level    */
  int levels;                    /* number of decomposition levels */
  vec current;                   /* a double buffer of vectors     */
  vec next;                      /* for temporary processing       */

} it_wavelet_t;

#define IT_WAVELET(x) IT_CAST(it_wavelet_t, x)

it_instanciate(it_wavelet_t);

static inline it_wavelet_t *it_wavelet_new(it_wavelet_lifting_t const *lifting, int level)
{
  return(it_new_va(it_wavelet_t)(it_va, lifting, level));
}

/* 1D wavelet transform and inverse transform */
#define it_wavelet_transform(t, v) ((vec) it_transform(IT_WAVELET(t), (vec) v))
#define it_wavelet_itransform(t, v) ((vec) it_itransform(IT_WAVELET(t), (vec) v))
vec it_dwt(vec v, it_wavelet_lifting_t const *lifting, int levels);
vec it_idwt(vec t, it_wavelet_lifting_t const *lifting, int levels);

/*--------------------------------------------------------------------*/
/* Split a vector containing all wavelet coefficients into several vectors,
   one for each subband.
   All the memory is allocated inside.                                      */
vec * it_wavelet_split( vec wav, int nb_levels );
vec it_wavelet_merge( vec * subbands, int nb_level );

#ifdef __cplusplus
}
#endif /* extern "C" */
#endif
