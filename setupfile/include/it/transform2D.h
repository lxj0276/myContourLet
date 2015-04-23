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

#ifndef __it_transform2D_h
#define __it_transform2D_h

#include "it/types.h"
#include "it/mat.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _it_transform2D_ {
  it_extends(it_object_t);

  Mat (* transform)(struct _it_transform2D_ *transform, Mat image);
  Mat (* itransform)(struct _it_transform2D_ *transform, Mat image);

  /* get the size of the output matrix */
  /* this is equal to input_size for critically sampled transforms */
  void (* get_output_size)(struct _it_transform2D_ *transform, size_t *width, size_t *height);

  /* fix the length of the vectors on which to apply the transform */
  void (* set_size)(struct _it_transform2D_ *transform, size_t width, size_t height);
  /* get the length of the vectors on which to apply the transform */
  /* 0 means automatic allocation/deletion during tranform/itransform */
  void (* get_size)(struct _it_transform2D_ *transform, size_t *width, size_t *height);

} it_transform2D_t;

#define IT_TRANSFORM2D(x) IT_CAST(it_transform2D_t, x)

static inline it_instanciate(it_transform2D_t)
{
  it_construct(it_object_t);
  it_set_magic(it_this, it_transform2D_t);
  return(it_this);
}

#define it_transform2D(t, image) __it_transform2D(IT_TRANSFORM2D(t), (Mat) image)
static inline Mat __it_transform2D(it_transform2D_t *t, Mat image)
{
  return(t->transform(t, image));
}

#define it_itransform2D(t, image) __it_itransform2D(IT_TRANSFORM2D(t), (Mat) image)
static inline Mat __it_itransform2D(it_transform2D_t *t, Mat image)
{
  return(t->itransform(t, image));
}

#define it_transform2D_get_output_size(t, w, h) \
            __it_transform2D_get_output_size(IT_TRANSFORM2D(t, w, h))
static inline void __it_transform2D_get_output_size(it_transform2D_t *t,
						    size_t *w, size_t *h)
{
  t->get_output_size(t, w, h);
}

#define it_transform2D_set_size(t, w, h) \
            __it_transform2D_set_size(IT_TRANSFORM2D(t), w, h)
static inline void __it_transform2D_set_size(it_transform2D_t *t,
					     size_t w, size_t h)
{
  t->set_size(t, w, h);
}

#define it_transform2D_get_size(t, w, h) \
            __it_transform_get_size(IT_TRANSFORM2D(t), l)
static inline void __it_transform2D_get_size(it_transform2D_t *t, size_t *w, size_t *h)
{
  t->get_size(t, w, h);
}

#define it_transform2D_clear_size(t) it_transform2D_set_size(t, 0, 0)


#ifdef __cplusplus
}
#endif /* extern "C" */
#endif
