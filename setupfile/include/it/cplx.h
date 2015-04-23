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
  Complex type
  Copyright (C) 2005 Vivien Chappelier
*/


#ifndef __it_cplx_h
#define __it_cplx_h

#include "it/types.h"
#include <math.h>
#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

/*---------------------------------------------------------------------------*/
/*               Complex type                                                */
/*---------------------------------------------------------------------------*/

typedef struct _cplx_ {
  double r; /* real part      */
  double i; /* imaginary part */
} cplx;

#define cplx(r, i) { r, i }
#define creal(c) ((c).r)
#define cimag(c) ((c).i)

static inline cplx cadd(cplx const a, cplx const b)
{
  cplx r;
  r.r = a.r + b.r;
  r.i = a.i + b.i;
  return(r);
}


static inline cplx csub(cplx const a, cplx const b)
{
  cplx r;
  r.r = a.r - b.r;
  r.i = a.i - b.i;
  return(r);
}


static inline cplx cmul(cplx const a, cplx const b)
{
  cplx r;
  r.r = a.r * b.r - a.i * b.i;
  r.i = a.r * b.i + a.i * b.r;
  return(r);
}


static inline cplx cdiv(cplx const a, cplx const b)
{
  cplx r;
  double n = b.r * b.r + b.i * b.i;

  assert(n);
  r.r = (a.r * b.r + a.i * b.i) / n;
  r.i = (a.i * b.r - a.r * b.i) / n;
  return(r);
}


static inline cplx cconj(cplx const a)
{
  cplx r;
  r.r =  a.r;
  r.i = -a.i;
  return( r );
}


static inline cplx cneg( cplx const a )
{
  cplx r;
  r.r = -a.r;
  r.i = -a.i;
  return( r );
}


static inline cplx cinv(cplx const a)
{
  double n = a.r * a.r + a.i * a.i;
  cplx r;
  r.r =  a.r / n;
  r.i = -a.i / n;
  return(r);
}


static inline int ceq( cplx const a, cplx const b )
{
  return( a.r == b.r && a.i == b.i );
}


static inline double cnorm( cplx const a )
{
  return( sqrt( a.r * a.r + a.i * a.i ) );
}

/* some constants */
extern cplx const cplx_0;
extern cplx const cplx_1;
extern cplx const cplx_I;
#define cplx_zero cplx_0
#define cplx_one cplx_1

#ifdef __cplusplus
}
#endif /* extern "C" */
#endif

