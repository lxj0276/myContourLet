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

#ifndef __it_wavelet2D_h
#define __it_wavelet2D_h

#include "it/types.h"
#include "it/mat.h"
#include "it/wavelet.h"
#include "it/transform2D.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _it_wavelet2D_ {
  it_extends(it_transform2D_t);

  it_wavelet_lifting_t *lifting;
  int level;                     /* current decomposition level */
  int levels;                    /* number of decomposition levels */
  int width, height;             /* widthxheight of the original frame */
  vec buffer;

  void (* it_overloaded(destructor))(it_object_t *it_this);

  int (* copy)(struct _it_wavelet2D_ *it_this,
	       struct _it_wavelet2D_*source);

} it_wavelet2D_t;

#define IT_WAVELET2D(x) IT_CAST(it_wavelet2D_t, x)

it_instanciate(it_wavelet2D_t);

static inline it_wavelet2D_t *it_wavelet2D_new(it_wavelet_lifting_t const *lifting, int level)
{
  return(it_new_va(it_wavelet2D_t)(it_va, lifting, level));
}

#define it_wavelet2D_copy(a, b)  a->copy(a, b)

/* 2D wavelet transform and inverse transform */
#define it_wavelet2D_transform(t, m) ((mat) it_transform2D(IT_WAVELET2D(t), m))
#define it_wavelet2D_itransform(t, m) ((mat) it_itransform2D(IT_WAVELET2D(t), m))
/* same thing with internal object construction/destruction */
mat it_dwt2D(mat m, it_wavelet_lifting_t const *lifting, int levels);
mat it_idwt2D(mat t, it_wavelet_lifting_t const *lifting, int levels);

/*--------------------------------------------------------------------*/
/* Split a matrix containing all wavelet components into several matrices. 
   All the memory is allocated inside.                                      */
mat * it_wavelet2D_split( mat wav, int nb_levels );
mat it_wavelet2D_merge( mat * subbands, int nb_level );

#ifdef __cplusplus
}
#endif /* extern "C" */
#endif
