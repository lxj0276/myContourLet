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
  1D fourier transform using lifting.
  Copyright (C) 2005 Vivien Chappelier
*/

#ifndef __it_fourier_h
#define __it_fourier_h

#include <it/types.h>
#include <it/mat.h>
#include <it/transform.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _it_fourier_ {
  it_extends(it_transform_t);

  void (* it_overloaded(destructor))(it_object_t *it_this);

  size_t length;                    /* length of the vector to decompose */
  cvec roots;                    /* the unit roots                    */

} it_fourier_t;

#define IT_FOURIER(x) IT_CAST(it_fourier_t, x)

it_instanciate(it_fourier_t);

static inline it_fourier_t *it_fourier_new(void)
{
  return(it_new_va(it_fourier_t)(it_va));
}

/* 1D fourier transform and inverse transform */
#define it_fourier_transform(t, v) ((cvec) it_transform(IT_FOURIER(t), (cvec) v))
#define it_fourier_itransform(t, v) ((cvec) it_itransform(IT_FOURIER(t), (cvec) v))

static inline cvec it_dft(cvec v) {
  it_fourier_t *fourier;
  cvec vt;

  fourier = it_fourier_new();
  vt = it_fourier_transform(fourier, v);
  it_delete(fourier);
  return(vt);
}

static inline cvec it_idft(cvec v) {
  it_fourier_t *fourier;
  cvec vt;

  fourier = it_fourier_new();
  vt = it_fourier_itransform(fourier, v);
  it_delete(fourier);
  return(vt);
}

#ifdef __cplusplus
}
#endif /* extern "C" */
#endif
