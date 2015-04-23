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
  Transforms
  Copyright (C) 2005 Vivien Chappelier
*/

#ifndef __it_transform_h
#define __it_transform_h

#include "it/types.h"
#include "it/mat.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _it_transform_ {
  it_extends(it_object_t);

  Vec (* transform)(struct _it_transform_ *transform, Vec v);
  Vec (* itransform)(struct _it_transform_ *transform, Vec V);

  /* get the length of the output vector */
  /* this is equal to input_size for critically sampled transforms */
  void (* get_output_size)(struct _it_transform_ *transform, size_t *input_size);

  /* fix the length of the vectors on which to apply the transform */
  void (* set_size)(struct _it_transform_ *transform, size_t length);
  /* get the length of the vectors on which to apply the transform */
  /* 0 means automatic allocation/deletion during tranform/itransform */
  void (* get_size)(struct _it_transform_ *transform, size_t *length);

} it_transform_t;

#define IT_TRANSFORM(x) IT_CAST(it_transform_t, x)

static inline it_instanciate(it_transform_t)
{
  it_construct(it_object_t);
  it_set_magic(it_this, it_transform_t);
  return(it_this);
}

#define it_transform(t, v) __it_transform(IT_TRANSFORM(t), (Vec) v)
static inline Vec __it_transform(it_transform_t *t, Vec v)
{
  return(t->transform(t, v));
}

#define it_itransform(t, V) __it_itransform(IT_TRANSFORM(t), (Vec) V)
static inline Vec __it_itransform(it_transform_t *t, Vec V)
{
  return(t->itransform(t, V));
}

#define it_transform_get_output_size(t, l) \
            __it_transform_get_output_size(IT_TRANSFORM(t), l)
static inline void __it_transform_get_output_size(it_transform_t *t, size_t *l)
{
  t->get_output_size(t, l);
}

#define it_transform_set_size(t, l) __it_transform_set_size(IT_TRANSFORM(t), l)
static inline void __it_transform_set_size(it_transform_t *t, size_t l)
{
  t->set_size(t, l);
}

#define it_transform_get_size(t, l) __it_transform_get_size(IT_TRANSFORM(t), l)
static inline void __it_transform_get_size(it_transform_t *t, size_t *l)
{
  t->get_size(t, l);
}

#define it_transform_clear_size(t) it_transform_set_size(t, 0)

#ifdef __cplusplus
}
#endif /* extern "C" */
#endif
