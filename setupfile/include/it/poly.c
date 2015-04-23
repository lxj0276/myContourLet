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

/*! \defgroup poly   Polynomials related functions
  Copyright (C) 2005-2007 Vivien Chappelier
@{
*/


#include <it/math.h>
#include <it/poly.h>
#include <it/io.h>

/*! \brief remove null factors from the polynomial 

  e.g. 0 X^2 + X + 1 becomes X + 1 */
void poly_normalize (vec v)
{
  int  i;

  for (i = vec_length (v) - 1; i >= 0; i--)
    if (v[i])
      break;
  vec_set_length (v, i + 1);
}

/* multiply the polynomial by X^shift.
   Factors with negative exponent are discarded
   e.g. 2X^2 + 3X + 1 shifted -1 becomes 2X + 3.
*/
void poly_shift (vec v, int shift)
{
  int  i;

  if (shift < 0) {
    for (i = -shift; i < vec_length (v); i++)
      v[i + shift] = v[i];
  }
  else {
    for (i = vec_length (v); i >= shift; i--)
      v[i] = v[i - shift];
    for (; i >= 0; i--)
      v[i] = 0;
  }
}

/* check if the polynomial is null. */
int poly_is_null (vec v)
{
  int  i;

  for (i = 0; i < vec_length (v); i++)
    if (v[i])
      return (0);

  return (1);
}

/* polynomial addition */
vec poly_add (vec a, vec b)
{
  int  i;

  if (vec_length (b) > vec_length (a)) {
    for (i = 0; i < vec_length (a); i++)
      a[i] += b[i];
    vec_set_length (a, vec_length (b));
    for (; i < vec_length (b); i++) {
      a[i] = b[i];
    }
  }
  else {
    for (i = 0; i < vec_length (b); i++)
      a[i] += b[i];
  }
  poly_normalize (a);
  return (a);
}

/* polynomial subtraction */
vec poly_sub (vec a, vec b)
{
  int  i;

  if (vec_length (b) > vec_length (a)) {
    for (i = 0; i < vec_length (a); i++)
      a[i] -= b[i];
    vec_set_length (a, vec_length (b));
    for (; i < vec_length (b); i++) {
      a[i] = -b[i];
    }
  }
  else {
    for (i = 0; i < vec_length (b); i++)
      a[i] -= b[i];
  }
  poly_normalize (a);
  return (a);
}

/* Laurent polynomial Euclidean division of a by b */
/* This finds Q and R such that A = B Q + R X^deg_x */
/* for classical polynomial division, deg_x = 0. */
/* a and b are assumed to be normalized. */
/* The returned remainder is not normalized. */
vec _lpoly_ediv (vec a, vec b, int deg_x, vec * _q)
{
  int  deg_a, deg_b, deg_q;
  vec  q;
  int  i, j;

  if (poly_deg (a) < poly_deg (b)) {
    *_q = vec_null;
    return (a);
  }

  deg_a = poly_deg (a);
  deg_b = poly_deg (b);
  deg_q = deg_a - deg_b;

  q = vec_new (deg_q + 1);

  /* set upper terms of A to zero recursively */
  for (i = 0; i <= (deg_q - deg_x); i++) {

    /* compute Q that cancels A[deg_a - i] */
    q[deg_q - i] = a[deg_a - i] / b[deg_b];

    /* subtract */
    for (j = 0; j <= deg_b; j++)
      a[deg_q - i + j] -= q[deg_q - i] * b[j];
  }

  /* set lower terms of A to zero recursively */
  for (i = 0; i < deg_x; i++) {

    /* compute Q that cancels A[i] */
    q[i] = a[i] / b[0];

    /* subtract */
    for (j = 0; j <= deg_b; j++)
      a[i + j] -= q[i] * b[j];
  }

  *_q = q;

  return (a);
}

vec lpoly_ediv (vec _a, vec _b, int deg_x, vec * _q)
{
  vec  b, a;

  /* normalize input polynomials */
  b = vec_clone (_b);
  a = vec_clone (_a);
  poly_normalize (b);
  poly_normalize (a);

  it_assert (poly_deg (b) >= 0, "Polynomial division by zero\n");

  a = _lpoly_ediv (a, b, deg_x, _q);

  /* normalize the remainder */
  poly_normalize (a);

  vec_delete (b);
  return (a);
}

vec poly_gcd (vec _a, vec _b)
{
  vec  b, a, t;

  /* normalize input polynomials */
  b = vec_clone (_b);
  a = vec_clone (_a);
  poly_normalize (b);
  poly_normalize (a);

  while (!poly_is_null (b)) {
    t = b;
    a = poly_mod (a, b);
    b = a;
    a = t;
  }

  vec_delete (b);
  return (a);
}

it_function (itf_polynomial)
{
  int  i;
  double v;
  vec  poly = it_this->poly;

  v = poly[poly_deg (poly)];
  for (i = poly_deg (poly) - 1; i >= 0; i--)
    v = v * x + poly[i];

  return (v);
}

/*! @} */
