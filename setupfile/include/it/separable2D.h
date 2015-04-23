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
  Copyright (C) 2005 Vivien Chappelier
*/

#ifndef __it_separable2D_h
#define __it_separable2D_h

#include <it/types.h>
#include <it/mat.h>
#include <it/transform.h>
#include <it/transform2D.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _it_separable2D_ {
  it_extends(it_transform2D_t);

  it_transform_t *transform;

  void (* it_overloaded(destructor))(it_object_t *it_this);

} it_separable2D_t;

#define IT_SEPARABLE2D(x) IT_CAST(it_separable2D_t, x)

it_instanciate(it_separable2D_t);

#define it_separable2D_new(t) __it_separable2D_new(IT_TRANSFORM(t))
static inline it_separable2D_t *__it_separable2D_new(it_transform_t *transform)
{
  return(it_new_va(it_separable2D_t)(it_va, transform));
}

/* 2D separable transform and inverse transform */
#define it_separable2D_transform(t, m) ((mat) it_transform2D(IT_SEPARABLE2D(t), m))
#define it_separable2D_itransform(t, m) ((mat) it_itransform2D(IT_SEPARABLE2D(t), m))
/* same thing with internal object construction/destruction */
mat it_separable2D(mat m, it_transform_t);
mat it_iseparable2D(mat t, it_transform_t);

#ifdef __cplusplus
}
#endif /* extern "C" */
#endif
