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
  Vectors
  Copyright (C) 2005 Vivien Chappelier, Herve Jegou
*/


#include <it/math.h>
#include <it/vec.h>
#include <it/io.h>
#include <it/random.h>

/*---------------------------------------------------------------------------*/
/*                Constant vectors                                           */
/*---------------------------------------------------------------------------*/

/* empty vectors */
static Vec_header_t __vec_null = {
  0,				/* effective length of the vector (<= max_length)     */
  0,				/* amount of memory allocated for the vector elements */
  NULL,				/* memory block associated to this vector             */
  sizeof (double),		/* size of the stored elements                        */
};
static Vec_header_t __ivec_null = {
  0,				/* effective length of the vector (<= max_length)     */
  0,				/* amount of memory allocated for the vector elements */
  NULL,				/* memory block associated to this vector             */
  sizeof (int),			/* size of the stored elements                        */
};
static Vec_header_t __bvec_null = {
  0,				/* effective length of the vector (<= max_length)     */
  0,				/* amount of memory allocated for the vector elements */
  NULL,				/* memory block associated to this vector             */
  sizeof (byte),		/* size of the stored elements                        */
};
static Vec_header_t __cvec_null = {
  0,				/* effective length of the vector (<= max_length)     */
  0,				/* amount of memory allocated for the vector elements */
  NULL,				/* memory block associated to this vector             */
  sizeof (cplx),		/* size of the stored elements                        */
};


/* empty vectors for all types (+ 1 is to align at the end of the struct
   on the first element which is non-existent).                         */
vec const vec_null = (vec) (&__vec_null + 1);
ivec const ivec_null = (ivec) (&__ivec_null + 1);
bvec const bvec_null = (bvec) (&__bvec_null + 1);
cvec const cvec_null = (cvec) (&__cvec_null + 1);

/*---------------------------------------------------------------------------*/
/*                Vector allocation functions                                */
/*---------------------------------------------------------------------------*/

void *__Vec_new_alloc (size_t elem_size, idx_t length, idx_t length_max)
{
  Vec_header_t *hdr;
  char *ptr, *aligned;
  int  padding;

  /* allocate a vector of size 'length_max' with some extra room 
     for the header and padding                                  */
  ptr =
    (char *) malloc (sizeof (Vec_header_t) + length_max * elem_size +
		     IT_ALLOC_ALIGN);
  it_assert (ptr, "No enough memory to allocate the vector");

  /* make sure the first element is properly aligned */
  aligned = ptr + sizeof (Vec_header_t) + IT_ALLOC_ALIGN - 1;
  padding = ((int) (long) (aligned)) & (IT_ALLOC_ALIGN - 1);
  aligned -= padding;
  /* put the header before the first element and fill it */
  hdr = (Vec_header_t *) aligned - 1;
  hdr->length = length;
  hdr->length_max = length_max;
  hdr->ptr = ptr;
  hdr->element_size = elem_size;
  return (aligned);
}


vec vec_new_alloc (idx_t length, idx_t length_max) {
  assert (length_max >= length);
  return Vec_new_alloc (double, length, length_max);
}


ivec ivec_new_alloc (idx_t length, idx_t length_max) {
  assert (length_max >= length);
  return Vec_new_alloc (int, length, length_max);
}


bvec bvec_new_alloc (idx_t length, idx_t length_max) {
  assert (length_max >= length);
  return Vec_new_alloc (byte, length, length_max);
}


cvec cvec_new_alloc (idx_t length, idx_t length_max) {
  assert (length_max >= length);
  return Vec_new_alloc (cplx, length, length_max);
}



/* We basically had two choices to have a properly aligned realloc.
   One is to call realloc (which does a memory copy if the new block cannot
   fit at the location of the old one). In this case the returned memory block
   may not be aligned properly and another copy is required to ensure proper
   alignment. If the block could fit, no memory copy is required.
   The second method is to free the memory systematically and call a new malloc.
   Obviously this method is slower if the resized block could fit at the
   location of the old one, since it always requires a memory copy.
   However, due to the geometric reallocation procedure used throughout the
   code, a realloc has very high chances of returning a different pointer. 
   This has been verified in practice. 
   Therefore, we decided to use the free/malloc method.
*/
void *__Vec_new_realloc (void *V, size_t elem_size, idx_t length,
			 idx_t length_max)
{
  void *new_Vec;
  idx_t old_length;

  new_Vec = __Vec_new_alloc (elem_size, length, length_max);
  assert (new_Vec);
  if (V) {
    old_length = Vec_length (V);
    memcpy (new_Vec, V, old_length * elem_size);
    Vec_delete (V);
  }
  return (new_Vec);
}


Vec __Vec_new (size_t elem_size, idx_t N) 
{
  return (__Vec_new_alloc (elem_size, N, N));
}


vec vec_new (idx_t length) 
{
  return Vec_new (double, length);
}


ivec ivec_new (idx_t length) 
{
  return Vec_new (int, length);
}


bvec bvec_new (idx_t length) 
{
  return Vec_new (byte, length);
}


cvec cvec_new (idx_t length) 
{
  return Vec_new (cplx, length);
}


void vec_delete (vec v) 
{
  Vec_delete (v);
}


void ivec_delete (ivec v) 
{
  Vec_delete (v);
}


void bvec_delete (bvec v) 
{
  Vec_delete (v);
}


void cvec_delete (cvec v) 
{
  Vec_delete (v);
}


/*---------------------------------------------------------------------------*/
/** Length and max length vector operations                            */
/*---------------------------------------------------------------------------*/

idx_t vec_length (vec v) 
{
  return Vec_length (v);
}


idx_t ivec_length (ivec v) 
{
  return Vec_length (v);
}


idx_t bvec_length (bvec v) 
{
  return Vec_length (v);
}


idx_t cvec_length (cvec v) 
{
  return Vec_length (v);
}


idx_t vec_length_max (vec v) 
{
  return Vec_length_max (v);
}


idx_t ivec_length_max (ivec v) 
{
  return Vec_length_max (v);
}


idx_t bvec_length_max (bvec v) 
{
  return Vec_length_max (v);
}


idx_t cvec_length_max (cvec v) 
{
  return Vec_length_max (v);
}


vec _vec_set_length_max (vec v, idx_t N) 
{
  Vec_set_length_max (v, N);
  return v;
}


ivec _ivec_set_length_max (ivec v, idx_t N) 
{
  Vec_set_length_max (v, N);
  return v;
}


bvec _bvec_set_length_max (bvec v, idx_t N) 
{
  Vec_set_length_max (v, N);
  return v;
}


cvec _cvec_set_length_max (cvec v, idx_t N) 
{
  Vec_set_length_max (v, N);
  return v;
}


vec _vec_set_length (vec v, idx_t N) 
{
  Vec_set_length (v, N);
  return v;
}


ivec _ivec_set_length (ivec v, idx_t N) 
{
  Vec_set_length (v, N);
  return v;
}


bvec _bvec_set_length (bvec v, idx_t N) 
{
  Vec_set_length (v, N);
  return v;
}


cvec _cvec_set_length (cvec v, idx_t N) 
{
  Vec_set_length (v, N);
  return v;
}


/*---------------------------------------------------------------------------*/
/** @name Initialization of a vector or of some of its components            */
/*---------------------------------------------------------------------------*/

vec vec_init (vec v, double *buf, idx_t N) 
{
  return (vec) Vec_init (v, buf, N);
}


ivec ivec_init (ivec v, int *buf, idx_t N) 
{
  return (ivec) Vec_init (v, buf, N);
}


bvec bvec_init (bvec v, byte * buf, idx_t N) 
{
  return (bvec) Vec_init (v, buf, N);
}


cvec cvec_init (cvec v, cplx * buf, idx_t N) 
{
  return (cvec) Vec_init (v, buf, N);
}


vec vec_set (vec v, double val) 
{
  Vec_set (v, val);
  return v;
}


ivec ivec_set (ivec v, int val) 
{
  Vec_set (v, val);
  return v;
}


bvec bvec_set (bvec v, byte val) 
{
  Vec_set (v, val);
  return v;
}


cvec cvec_set (cvec v, cplx val) 
{
  Vec_set (v, val);
  return v;
}


vec vec_set_between (vec v, idx_t i1, idx_t i2, double val) 
{
  Vec_set_between (v, i1, i2, val);
  return v;
}


ivec ivec_set_between (ivec v, idx_t i1, idx_t i2, int val) 
{
  Vec_set_between (v, i1, i2, val);
  return v;
}


bvec bvec_set_between (bvec v, idx_t i1, idx_t i2, byte val) 
{
  Vec_set_between (v, i1, i2, val);
  return v;
}


cvec cvec_set_between (cvec v, idx_t i1, idx_t i2, cplx val) 
{
  Vec_set_between (v, i1, i2, val);
  return v;
}


void vec_set_subvector (vec v, vec s, idx_t idx) 
{
  VEC_END_PARAM (v, idx);
  Vec_set_subvector (v, s, idx);
}


void ivec_set_subvector (ivec v, ivec s, idx_t idx) 
{
  VEC_END_PARAM (v, idx);
  Vec_set_subvector (v, s, idx);
}


void bvec_set_subvector (bvec v, bvec s, idx_t idx) 
{
  VEC_END_PARAM (v, idx);
  Vec_set_subvector (v, s, idx);
}


void cvec_set_subvector (cvec v, cvec s, idx_t idx) 
{
  VEC_END_PARAM (v, idx);
  Vec_set_subvector (v, s, idx);
}


vec vec_get_subvector (vec v, idx_t i1, idx_t i2) 
{
  vec  s;
  VEC_END_PARAM (v, i1);
  VEC_END_PARAM (v, i2);
  s = vec_new (i2 - i1 + 1);
  vec_init (s, v + i1, i2 - i1 + 1);
  return s;
}


ivec ivec_get_subvector (ivec v, idx_t i1, idx_t i2) 
{
  ivec s;
  VEC_END_PARAM (v, i1);
  VEC_END_PARAM (v, i2);
  s = ivec_new (i2 - i1 + 1);
  ivec_init (s, v + i1, i2 - i1 + 1);
  return s;
}


bvec bvec_get_subvector (bvec v, idx_t i1, idx_t i2) 
{
  bvec s;
  VEC_END_PARAM (v, i1);
  VEC_END_PARAM (v, i2);
  s = bvec_new (i2 - i1 + 1);
  bvec_init (s, v + i1, i2 - i1 + 1);
  return s;
}


cvec cvec_get_subvector (cvec v, idx_t i1, idx_t i2) 
{
  cvec s;
  VEC_END_PARAM (v, i1);
  VEC_END_PARAM (v, i2);
  s = cvec_new (i2 - i1 + 1);
  cvec_init (s, v + i1, i2 - i1 + 1);
  return s;
}



/*------------------------------------------------------------------------------*/
/*                Copy and Conversions Functions                                */
/*------------------------------------------------------------------------------*/

Vec __Vec_copy (Vec v1, Vec v2) 
{
  assert (v1);
  assert (v2);
  assert (Vec_element_size (v1) == Vec_element_size (v2));
  assert (Vec_length (v1) == Vec_length (v2));
  memcpy (v1, v2, Vec_element_size (v1) * Vec_length (v1));
  return (v1);
}


void vec_copy (vec dest, vec orig) 
{
  Vec_copy (dest, orig);
}


void ivec_copy (ivec dest, ivec orig) 
{
  Vec_copy (dest, orig);
}


void bvec_copy (bvec dest, bvec orig) 
{
  Vec_copy (dest, orig);
}


void cvec_copy (cvec dest, cvec orig) 
{
  Vec_copy (dest, orig);
}


Vec Vec_clone (Vec v) 
{
  assert (v);
  return (Vec_copy (__Vec_new (Vec_element_size (v), Vec_length (v)), v));
}


vec vec_clone (vec v) 
{
  assert (v);
  return ((vec) Vec_copy (vec_new (vec_length (v)), v));
}


ivec ivec_clone (ivec v) 
{
  assert (v);
  return ((ivec) Vec_copy (ivec_new (ivec_length (v)), v));
}


bvec bvec_clone (bvec v) 
{
  assert (v);
  return ((bvec) Vec_copy (bvec_new (bvec_length (v)), v));
}


cvec cvec_clone (cvec v) 
{
  assert (v);
  return ((cvec) Vec_copy (cvec_new (cvec_length (v)), v));
}


void vec_copy_from_ivec (vec dest, ivec orig)
{
  int i;
  assert (vec_length(dest) == ivec_length (orig));

  for (i = 0 ; i < vec_length (dest) ; i++ )
    dest[i] = orig[i];
}


void vec_copy_from_bvec (vec dest, bvec orig) 
{
  int i;
  assert (vec_length(dest) == bvec_length (orig));

  for (i = 0 ; i < vec_length (dest) ; i++ )
    dest[i] = orig[i];
}


void vec_copy_from_cvec (vec dest, cvec orig) 
{
  int i;
  assert (vec_length(dest) == cvec_length (orig));

  for (i = 0 ; i < vec_length (dest) ; i++ )
    dest[i] = creal (orig[i]);
}


void ivec_copy_from_vec (ivec dest, vec orig) 
{
  int i;
  assert (ivec_length(dest) == vec_length (orig));

  for (i = 0 ; i < ivec_length (dest) ; i++ )
    dest[i] = (int) orig[i];
}


void ivec_copy_from_bvec (ivec dest, bvec orig) 
{
  int i;
  assert (ivec_length(dest) == bvec_length (orig));

  for (i = 0 ; i < ivec_length (dest) ; i++ )
    dest[i] = orig[i];
}


void ivec_copy_from_cvec (ivec dest, cvec orig) 
{
  int i;
  assert (ivec_length(dest) == cvec_length (orig));

  for (i = 0 ; i < ivec_length (dest) ; i++ )
    dest[i] = (int) creal (orig[i]);
}


void bvec_copy_from_vec (bvec dest, vec orig) 
{
  int i;
  assert (bvec_length(dest) == vec_length (orig));

  for (i = 0 ; i < bvec_length (dest) ; i++ )
    dest[i] = (byte) orig[i];
}


void bvec_copy_from_ivec (bvec dest, ivec orig) 
{
  int i;
  assert (bvec_length(dest) == ivec_length (orig));

  for (i = 0 ; i < bvec_length (dest) ; i++ )
    dest[i] = (byte) orig[i];
}


void bvec_copy_from_cvec (bvec dest, cvec orig) 
{
  int i;
  assert (bvec_length(dest) == cvec_length (orig));

  for (i = 0 ; i < bvec_length (dest) ; i++ )
    dest[i] = (byte) creal (orig[i]);
}



void cvec_copy_from_vec (cvec dest, vec orig) 
{
  int i;
  assert (cvec_length(dest) == vec_length (orig));

  for (i = 0 ; i < cvec_length (dest) ; i++ ) {
    dest[i].r = orig[i];
    dest[i].i = 0;
  }
}


void cvec_copy_from_ivec (cvec dest, ivec orig) 
{
  int i;
  assert (cvec_length(dest) == ivec_length (orig));

  for (i = 0 ; i < cvec_length (dest) ; i++ ) {
    dest[i].r = (int) orig[i];
    dest[i].i = 0;
  }
}


void cvec_copy_from_bvec (cvec dest, bvec orig) 
{
  int i;
  assert (cvec_length(dest) == bvec_length (orig));

  for (i = 0 ; i < cvec_length (dest) ; i++ ) {
    dest[i].r = (byte) orig[i];
    dest[i].i = 0;
  }
}




/*------------------------------------------------------------------------------*/
void vec_copy_mem (double *buf, vec v)
{
  idx_t i;
  assert (v);
  assert (buf);
  for (i = 0; i < vec_length (v); i++)
    buf[i] = v[i];
}


/*------------------------------------------------------------------------------*/
void ivec_copy_mem (int *buf, ivec v)
{
  idx_t i;
  assert (v);
  assert (buf);
  for (i = 0; i < ivec_length (v); i++)
    buf[i] = v[i];
}


/*------------------------------------------------------------------------------*/
void bvec_copy_mem (byte * buf, bvec v)
{
  idx_t i;
  assert (v);
  assert (buf);
  for (i = 0; i < bvec_length (v); i++)
    buf[i] = v[i];
}


/*--------------------------------------------------------------------*/
void cvec_copy_mem (cplx * buf, cvec v)
{
  idx_t i;
  assert (v);
  assert (buf);
  for (i = 0; i < cvec_length (v); i++)
    buf[i] = v[i];
}

/*--------------------------------------------------------------------*/
void bvec_pack (byte * buf, bvec v)
{
  idx_t i, j;
  idx_t l = bvec_length (v);
  byte t;
  assert (v);
  assert (buf);
  for (i = 0; i < l / 8; i++) {
    t = 0;
    for (j = 0; j < 8; j++)
      t |= v[8 * i + j] << (7 - j);
    buf[i] = t;
  }
  t = 0;
  for (j = 0; j < l % 8; j++)
    t |= v[8 * i + j] << (7 - j);
  if (l % 8)
    buf[i] = t;
}

void bvec_unpack (bvec v, byte * buf)
{
  idx_t i, j;
  idx_t l = bvec_length (v);
  byte t = 0;
  assert (v);
  assert (buf);
  for (i = 0; i < l / 8; i++) {
    t = buf[i];
    for (j = 0; j < 8; j++)
      v[8 * i + j] = (t >> (7 - j)) & 1;
  }
  if (l % 8)
    t = buf[i];

  for (j = 0; j < l % 8; j++)
    v[8 * i + j] = (t >> (7 - j)) & 1;
}


/*--------------------------------------------------------------------*/
vec ivec_to_vec (ivec v)
{
  idx_t i;
  idx_t l = ivec_length (v);
  vec  r = vec_new (l);

  for (i = 0; i < l; i++)
    r[i] = (double) v[i];

  return (r);
}


bvec ivec_to_bvec (ivec v)
{
  idx_t i;
  idx_t l = ivec_length (v);
  bvec r = bvec_new (l);

  for (i = 0; i < l; i++)
    r[i] = (byte) v[i];

  return (r);
}

cvec ivec_to_cvec (ivec v)
{
  idx_t i;
  idx_t l = ivec_length (v);
  cvec r = cvec_new (l);

  for (i = 0; i < l; i++) {
    creal (r[i]) = v[i];
    cimag (r[i]) = 0;
  }

  return (r);
}


/*--------------------------------------------------------------------*/
ivec bvec_to_ivec (bvec v)
{
  idx_t i;
  idx_t l = bvec_length (v);
  ivec r = ivec_new (l);

  for (i = 0; i < l; i++)
    r[i] = (int) v[i];

  return (r);
}


vec bvec_to_vec (bvec v)
{
  idx_t i;
  idx_t l = bvec_length (v);
  vec  r = vec_new (l);

  for (i = 0; i < l; i++)
    r[i] = (double) v[i];

  return (r);
}

cvec bvec_to_cvec (bvec v)
{
  idx_t i;
  idx_t l = bvec_length (v);
  cvec r = cvec_new (l);

  for (i = 0; i < l; i++) {
    creal (r[i]) = v[i];
    cimag (r[i]) = 0;
  }

  return (r);
}

/*--------------------------------------------------------------------*/
ivec vec_to_ivec (vec v)
{
  idx_t i;
  idx_t l = vec_length (v);
  ivec r = ivec_new (l);

  for (i = 0; i < l; i++)
    r[i] = (int) v[i];

  return (r);
}


bvec vec_to_bvec (vec v)
{
  idx_t i;
  idx_t l = vec_length (v);
  bvec r = bvec_new (l);

  for (i = 0; i < l; i++)
    r[i] = (byte) v[i];

  return (r);
}

cvec vec_to_cvec (vec v)
{
  idx_t i;
  idx_t l = vec_length (v);
  cvec r = cvec_new (l);

  for (i = 0; i < l; i++) {
    creal (r[i]) = v[i];
    cimag (r[i]) = 0;
  }

  return (r);
}


/*------------------------------------------------------------------------------*/
/* Stack operations                                                             */
/*------------------------------------------------------------------------------*/

/* Set operations */
vec vec_del (vec v, idx_t pos) 
{
  Vec_del (v, pos);
  return v;
}


ivec ivec_del (ivec v, idx_t pos) 
{
  Vec_del (v, pos);
  return v;
}


bvec bvec_del (bvec v, idx_t pos) 
{
  Vec_del (v, pos);
  return v;
}


cvec cvec_del (cvec v, idx_t pos) 
{
  Vec_del (v, pos);
  return v;
}


/* To be used as v = ivec_ins( v, pos, elt ) */
vec _vec_ins (vec v, idx_t pos, double elt) 
{
  Vec_ins (v, pos, elt);
  return v;
}


ivec _ivec_ins (ivec v, idx_t pos, int elt) 
{
  Vec_ins (v, pos, elt);
  return v;
}


bvec _bvec_ins (bvec v, idx_t pos, byte elt) 
{
  Vec_ins (v, pos, elt);
  return v;
}


cvec _cvec_ins (cvec v, idx_t pos, cplx elt) 
{
  Vec_ins (v, pos, elt);
  return v;
}


vec _vec_push (vec v, double elt) 
{
  vec_ins (v, vec_length (v), elt);
  return v;
}


ivec _ivec_push (ivec v, int elt) 
{
  ivec_ins (v, ivec_length (v), elt);
  return v;
}


bvec _bvec_push (bvec v, byte elt) 
{
  bvec_ins (v, bvec_length (v), elt);
  return v;
}


cvec _cvec_push (cvec v, cplx elt) {
  cvec_ins (v, cvec_length (v), elt);
  return v;
}


vec vec_pop (vec v) 
{
  Vec_pop (v);
  return v;
}


ivec ivec_pop (ivec v) 
{
  Vec_pop (v);
  return v;
}


bvec bvec_pop (bvec v) 
{
  Vec_pop (v);
  return v;
}


cvec cvec_pop (cvec v) 
{
  Vec_pop (v);
  return v;
}


double vec_head (vec v) 
{
  return Vec_head (v);
}


int ivec_head (ivec v) 
{
  return Vec_head (v);
}


byte bvec_head (bvec v) 
{
  return Vec_head (v);
}


cplx cvec_head (const cvec v) 
{
  return Vec_head (v);
}


/*--------------------------------------------------------------------*/
/* Safe vector access                                                 */
/*--------------------------------------------------------------------*/

double *__vec (vec v, idx_t i) 
{
  assert (v);
  VEC_END_PARAM (v, i);
  assert (i >= 0 && i < vec_length (v));
  return (&v[i]);
}


int *__ivec (ivec v, idx_t i) 
{
  assert (v);
  VEC_END_PARAM (v, i);
  assert (i >= 0 && i < ivec_length (v));
  return (&v[i]);
}


byte *__bvec (bvec v, idx_t i) 
{
  assert (v);
  VEC_END_PARAM (v, i);
  assert (i >= 0 && i < bvec_length (v));
  return (&v[i]);
}


cplx *__cvec (cvec v, idx_t i) 
{
  assert (v);
  VEC_END_PARAM (v, i);
  assert (i >= 0 && i < cvec_length (v));
  return (&v[i]);
}


/*--------------------------------------------------------------------*/
/*                Comparisons functions                               */
/*--------------------------------------------------------------------*/

int vec_eq (vec v1, vec v2)
{
  idx_t i;
  assert (v1);
  assert (v2);

  if (vec_length (v1) != vec_length (v2))
    return 0;

  for (i = 0; i < vec_length (v1); i++)
    if (v1[i] != v2[i])
      return 0;
  return 1;
}


/*-----------------------------------------------------------------*/
int ivec_eq (ivec v1, ivec v2)
{
  idx_t i;
  assert (v1);
  assert (v2);

  if (ivec_length (v1) != ivec_length (v2))
    return 0;

  for (i = 0; i < ivec_length (v1); i++)
    if (v1[i] != v2[i])
      return 0;
  return 1;
}


/*-----------------------------------------------------------------*/
int bvec_eq (bvec v1, bvec v2)
{
  idx_t i;
  assert (v1);
  assert (v2);

  if (bvec_length (v1) != bvec_length (v2))
    return 0;

  for (i = 0; i < bvec_length (v1); i++)
    if (v1[i] != v2[i])
      return 0;
  return 1;
}


/*-----------------------------------------------------------------*/
int cvec_eq (cvec v1, cvec v2)
{
  idx_t i;
  assert (v1);
  assert (v2);

  if (cvec_length (v1) != cvec_length (v2))
    return 0;

  for (i = 0; i < cvec_length (v1); i++)
    if (!ceq (v1[i], v2[i]))
      return 0;
  return 1;
}


/*-----------------------------------------------------------------*/
int vec_geq (vec v1, vec v2)
{
  idx_t i, minl;
  assert (v1);
  assert (v2);
  minl =
    (vec_length (v1) > vec_length (v2) ? vec_length (v2) : vec_length (v1));

  for (i = 0; i < minl; i++)
    if (v1[i] > v2[i])
      return 0;

  if (vec_length (v1) >= vec_length (v2))
    return 1;
  else
    return 0;
}


/*-----------------------------------------------------------------*/
int ivec_geq (ivec v1, ivec v2)
{
  idx_t i, minl;
  assert (v1);
  assert (v2);
  minl =
    (ivec_length (v1) >
     ivec_length (v2) ? ivec_length (v2) : ivec_length (v1));

  for (i = 0; i < minl; i++)
    if (v1[i] > v2[i])
      return 0;

  if (ivec_length (v1) >= ivec_length (v2))
    return 1;
  else
    return 0;
}


/*------------------------------------------------------------------------------*/
int bvec_geq (bvec v1, bvec v2)
{
  idx_t i, minl;
  assert (v1);
  assert (v2);
  minl =
    (bvec_length (v1) >
     bvec_length (v2) ? bvec_length (v2) : bvec_length (v1));

  for (i = 0; i < minl; i++)
    if (v1[i] > v2[i])
      return 0;

  if (bvec_length (v1) >= bvec_length (v2))
    return 1;
  else
    return 0;
}




/*------------------------------------------------------------------------------*/
/*                Arithmetic functions                                          */
/*------------------------------------------------------------------------------*/


/* Operations with a scalar value                                               */
void vec_incr (vec v, double a)
{
  idx_t i;
  assert (v);
  for (i = 0; i < vec_length (v); i++)
    v[i] += a;
}


void vec_decr (vec v, double a)
{
  idx_t i;
  assert (v);
  for (i = 0; i < vec_length (v); i++)
    v[i] -= a;
}


void vec_mul_by (vec v, double a)
{
  idx_t i;
  assert (v);
  for (i = 0; i < vec_length (v); i++)
    v[i] *= a;
}


void vec_div_by (vec v, double a)
{
  idx_t i;
  assert (v);
  assert (a);
  a = 1 / a;
  for (i = 0; i < vec_length (v); i++)
    v[i] *= a;
}


/* Operations with a scalar value                                               */
void ivec_incr (ivec v, int a)
{
  idx_t i;
  assert (v);
  for (i = 0; i < ivec_length (v); i++)
    v[i] += a;
}


void ivec_decr (ivec v, int a)
{
  idx_t i;
  assert (v);
  for (i = 0; i < ivec_length (v); i++)
    v[i] -= a;
}


void ivec_mul_by (ivec v, int a)
{
  idx_t i;
  assert (v);
  for (i = 0; i < ivec_length (v); i++)
    v[i] *= a;
}


void ivec_div_by (ivec v, int a)
{
  idx_t i;
  assert (v);
  assert (a);
  for (i = 0; i < ivec_length (v); i++)
    v[i] /= a;
}


/* Operations with a scalar value                                               */
void cvec_incr_real (cvec v, double a)
{
  idx_t i;
  assert (v);
  for (i = 0; i < cvec_length (v); i++)
    creal (v[i]) += a;
}


void cvec_decr_real (cvec v, double a)
{
  idx_t i;
  assert (v);
  for (i = 0; i < cvec_length (v); i++)
    creal (v[i]) -= a;
}


void cvec_mul_by_real (cvec v, double a)
{
  idx_t i;
  assert (v);
  for (i = 0; i < cvec_length (v); i++) {
    v[i].r *= a;
    v[i].i *= a;
  }
}


void cvec_div_by_real (cvec v, double a)
{
  idx_t i;
  assert (v);
  assert (a);
  a = 1 / a;
  for (i = 0; i < cvec_length (v); i++) {
    v[i].r *= a;
    v[i].i *= a;
  }
}


void cvec_conj (cvec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < cvec_length (v); i++)
    cconj (v[i]);
}

void cvec_incr (cvec v, cplx a)
{
  idx_t i;
  assert (v);
  for (i = 0; i < cvec_length (v); i++)
    v[i] = cadd (v[i], a);
}


void cvec_decr (cvec v, cplx a)
{
  idx_t i;
  assert (v);
  for (i = 0; i < cvec_length (v); i++)
    v[i] = csub (v[i], a);
}


void cvec_mul_by (cvec v, cplx a)
{
  idx_t i;
  assert (v);
  for (i = 0; i < cvec_length (v); i++)
    v[i] = cmul (v[i], a);
}


void cvec_div_by (cvec v, cplx a)
{
  idx_t i;
  assert (v);
  assert (!ceq (a, cplx_0));
  a = cinv (a);
  for (i = 0; i < cvec_length (v); i++)
    v[i] = cmul (v[i], a);
}


/*------------------------------------------------------------------------------*/
/* Components per components operations (vectors must be of same size)    */
void vec_add (vec v1, vec v2)
{
  idx_t i;
  assert (v1);
  assert (v2);
  assert (vec_length (v1) == vec_length (v2));
  for (i = 0; i < vec_length (v1); i++)
    v1[i] += v2[i];
}


void vec_sub (vec v1, vec v2)
{
  idx_t i;
  assert (v1);
  assert (v2);
  assert (vec_length (v1) == vec_length (v2));
  for (i = 0; i < vec_length (v1); i++)
    v1[i] -= v2[i];
}


void vec_mul (vec v1, vec v2)
{
  idx_t i;
  assert (v1);
  assert (v2);
  assert (vec_length (v1) == vec_length (v2));
  for (i = 0; i < vec_length (v1); i++)
    v1[i] *= v2[i];
}


void vec_div (vec v1, vec v2)
{
  idx_t i;
  assert (v1);
  assert (v2);
  assert (vec_length (v1) == vec_length (v2));
  for (i = 0; i < vec_length (v1); i++)
    v1[i] /= v2[i];
}


void ivec_add (ivec v1, ivec v2)
{
  idx_t i;
  assert (v1);
  assert (v2);
  assert (ivec_length (v1) == ivec_length (v2));
  for (i = 0; i < ivec_length (v1); i++)
    v1[i] += v2[i];
}


void ivec_sub (ivec v1, ivec v2)
{
  idx_t i;
  assert (v1);
  assert (v2);
  assert (ivec_length (v1) == ivec_length (v2));
  for (i = 0; i < ivec_length (v1); i++)
    v1[i] -= v2[i];
}


void ivec_mul (ivec v1, ivec v2)
{
  idx_t i;
  assert (v1);
  assert (v2);
  assert (ivec_length (v1) == ivec_length (v2));
  for (i = 0; i < ivec_length (v1); i++)
    v1[i] *= v2[i];
}


void ivec_div (ivec v1, ivec v2)
{
  idx_t i;
  assert (v1);
  assert (v2);
  assert (ivec_length (v1) == ivec_length (v2));
  for (i = 0; i < ivec_length (v1); i++)
    v1[i] /= v2[i];
}


void cvec_add (cvec v1, cvec v2)
{
  idx_t i;
  assert (v1);
  assert (v2);
  assert (cvec_length (v1) == cvec_length (v2));
  for (i = 0; i < cvec_length (v1); i++)
    v1[i] = cadd (v1[i], v2[i]);
}


void cvec_sub (cvec v1, cvec v2)
{
  idx_t i;
  assert (v1);
  assert (v2);
  assert (cvec_length (v1) == cvec_length (v2));
  for (i = 0; i < cvec_length (v1); i++)
    v1[i] = csub (v1[i], v2[i]);
}


void cvec_mul (cvec v1, cvec v2)
{
  idx_t i;
  assert (v1);
  assert (v2);
  assert (cvec_length (v1) == cvec_length (v2));
  for (i = 0; i < cvec_length (v1); i++)
    v1[i] = cmul (v1[i], v2[i]);
}


void cvec_div (cvec v1, cvec v2)
{
  idx_t i;
  assert (v1);
  assert (v2);
  assert (cvec_length (v1) == cvec_length (v2));
  for (i = 0; i < cvec_length (v1); i++)
    v1[i] = cdiv (v1[i], v2[i]);
}

/* v1 .* v2' */
void cvec_conj_mul (cvec v1, cvec v2)
{
  idx_t i;
  assert (v1);
  assert (v2);
  assert (cvec_length (v1) == cvec_length (v2));
  for (i = 0; i < cvec_length (v1); i++)
    v1[i] = cmul (v1[i], cconj (v2[i]));
}

/*------------------------------------------------------------------------------*/
vec vec_new_add (vec v1, vec v2)
{
  vec  r = vec_clone (v1);
  vec_add (r, v2);
  return r;
}


vec vec_new_sub (vec v1, vec v2)
{
  vec  r = vec_clone (v1);
  vec_sub (r, v2);
  return r;
}


vec vec_new_mul (vec v1, vec v2)
{
  vec  r = vec_clone (v1);
  vec_mul (r, v2);
  return r;
}


vec vec_new_div (vec v1, vec v2)
{
  vec  r = vec_clone (v1);
  vec_div (r, v2);
  return r;
}


ivec ivec_new_add (ivec v1, ivec v2)
{
  ivec r = ivec_clone (v1);
  ivec_add (r, v2);
  return r;
}


ivec ivec_new_sub (ivec v1, ivec v2)
{
  ivec r = ivec_clone (v1);
  ivec_sub (r, v2);
  return r;
}


ivec ivec_new_mul (ivec v1, ivec v2)
{
  ivec r = ivec_clone (v1);
  ivec_mul (r, v2);
  return r;
}


ivec ivec_new_div (ivec v1, ivec v2)
{
  ivec r = ivec_clone (v1);
  ivec_div (r, v2);
  return r;
}


/*------------------------------------------------------------------------------*/
double vec_inner_product (vec v1, vec v2)
{
  double p = 0;
  idx_t i;
  assert (v1);
  assert (v2);
  assert (vec_length (v1) == vec_length (v2));
  for (i = 0; i < vec_length (v1); i++)
    p += v1[i] * v2[i];
  return p;
}

/* Robust version with Kahan's method */
double vec_inner_product_robust (vec v1, vec v2)
{
  idx_t i; 
  double p; 
  double c = 0., t, y;
  assert (v1);
  assert (v2);
  assert (vec_length (v1) == vec_length (v2)); 
  p  = v1[0]*v2[0];
  for ( i= 1; i< vec_length (v1); i++ )
    {
      y = v1[i]*v2[i]-c;
      t = p+y;
      c = (t-p) - y;
      p = t;
    }  
  return p;
}


/*------------------------------------------------------------------------------*/
int ivec_inner_product (ivec v1, ivec v2)
{
  int  p = 0;
  idx_t i;
  assert (v1);
  assert (v2);
  assert (ivec_length (v1) == ivec_length (v2));
  for (i = 0; i < ivec_length (v1); i++)
    p += v1[i] * v2[i];
  return p;
}


/*------------------------------------------------------------------------------*/
double vecivec_inner_product (vec v1, ivec v2)
{
  double p = 0;
  idx_t i;
  assert (v1);
  assert (v2);
  assert (vec_length (v1) == ivec_length (v2));
  for (i = 0; i < vec_length (v1); i++)
    p += v1[i] * v2[i];
  return p;
}


/*------------------------------------------------------------------------------*/
int bvecivec_inner_product (bvec v1, ivec v2)
{
  int  p = 0;
  idx_t i;
  assert (v1);
  assert (v2);
  assert (bvec_length (v1) == ivec_length (v2));
  for (i = 0; i < bvec_length (v1); i++)
    p += v1[i] * v2[i];
  return p;
}


/*------------------------------------------------------------------------------*/
void vec_neg (vec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < vec_length (v); i++)
    v[i] = -v[i];
}


void ivec_neg (ivec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < ivec_length (v); i++)
    v[i] = -v[i];
}


void cvec_neg (cvec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < cvec_length (v); i++)
    v[i] = cneg (v[i]);
}

/* elementwise |v_i|^2 */
void cvec_abssqr (cvec v)
{
  idx_t i;
  double re, im;
  assert (v);
  for (i = 0; i < cvec_length (v); i++) {
    re = creal (v[i]);
    im = cimag (v[i]);
    creal (v[i]) = re * re + im * im;
    cimag (v[i]) = 0;
  }
}

void vec_sqr (vec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < vec_length (v); i++)
    v[i] *= v[i];
}


void ivec_sqr (ivec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < ivec_length (v); i++)
    v[i] *= v[i];
}


void vec_sqrt (vec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < Vec_length (v); i++)
    v[i] = v[i]<0?v[i]:sqrt (v[i]); /* negative test added by FC for fastica */
}


void vec_log (vec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < Vec_length (v); i++)
    v[i] = log (v[i]);
}


void vec_log10 (vec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < Vec_length (v); i++)
    v[i] = log10 (v[i]);
}


void vec_exp (vec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < Vec_length (v); i++)
    v[i] = exp (v[i]);
}


void vec_abs (vec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < vec_length (v); i++)
    if (v[i] < 0)
      v[i] = -(v[i]);
}


void ivec_abs (ivec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < ivec_length (v); i++)
    if (v[i] < 0)
      v[i] = -(v[i]);
}


vec vec_new_abs (vec v)
{
  vec  r = vec_new (vec_length (v));
  vec_abs (r);
  return r;
}


ivec ivec_new_abs (ivec v)
{
  ivec r = ivec_new (ivec_length (v));
  ivec_abs (r);
  return r;
}


vec cvec_new_abs (cvec v)
{
  idx_t i;
  vec  va;
  assert (v);
  va = vec_new (cvec_length (v));
  for (i = 0; i < cvec_length (v); i++)
    va[i] = cnorm (v[i]);
  return va;
}


void vec_pow (vec v, double a)
{
  idx_t i;
  assert (v);
  for (i = 0; i < Vec_length (v); i++)
    v[i] = pow (v[i], a);
}


void vec_normalize (vec v, double nr)
{
  idx_t i;
  double s = 0;
  assert (v);

  /* For optimization purpose, the norm 1 is treated separately */
  if (nr == 1)
    s = vec_sum (v);
  else {
    for (i = 0; i < vec_length (v); i++)
      s += pow (fabs (v[i]), nr);
    s = pow (s, 1.0 / nr);
  }

  if (s == 0)
    return;

  for (i = 0; i < vec_length (v); i++)
    v[i] /= s;
}


vec vec_new_pow (vec v, double a)
{
  idx_t i;
  assert (v);
  vec r = vec_new (vec_length (v));

  for (i = 0; i < Vec_length (v); i++)
    r[i] = pow (v[i], a);

  return r;
}


vec vec_new_normalize (vec v, double nr)
{
  idx_t i;
  double s = 0;
  assert (v);
  vec r = vec_new (vec_length (v));

  /* For optimization purpose, the norm 1 is treated separately */
  if (nr == 1)
    s = vec_sum (v);
  else {
    for (i = 0; i < vec_length (v); i++)
      s += pow (fabs (v[i]), nr);
    s = pow (s, 1.0 / nr);
  }

  if (s == 0)
    return vec_new_set (1. / vec_length (v), vec_length (v));

  for (i = 0; i < vec_length (v); i++)
    r[i] = v[i] / s;

  return r;
}


int ivec_min (ivec v)
{
  idx_t i;
  int  m = INT_MAX;
  assert (v);
  for (i = 0; i < Vec_length (v); i++)
    if (v[i] < m)
      m = v[i];
  return m;
}


int ivec_max (ivec v)
{
  idx_t i;
  int  m = INT_MIN;
  assert (v);
  for (i = 0; i < Vec_length (v); i++)
    if (v[i] > m)
      m = v[i];
  return m;
}

idx_t ivec_min_index (ivec v)
{
  idx_t i, mi = NULL_INDEX;
  int  m = INT_MAX;
  assert (v);
  for (i = 0; i < ivec_length (v); i++)
    if (v[i] < m) {
      m = v[i];
      mi = i;
    }
  return mi;
}


idx_t ivec_max_index (ivec v)
{
  idx_t i, mi = NULL_INDEX;
  int  m = INT_MIN;
  assert (v);
  for (i = 0; i < ivec_length (v); i++)
    if (v[i] > m) {
      m = v[i];
      mi = i;
    }
  return mi;
}


double ivec_mean (ivec v)
{
  assert (v);
  return ivec_sum (v) / (double) ivec_length (v);
}


/*------------------------------------------------------------------------------*/
double vec_sum (vec v)
{
  idx_t i;
  double s = 0;
  assert (v);
  for (i = 0; i < Vec_length (v); i++)
    s += v[i];
  return s;
}

/* Robust version with Kahan's method */
double vec_sum_robust(vec v)
{

  idx_t i; 
  double s; 
  double c = 0., t, y;
  assert (v);
  s  = v[0];
  for ( i= 1; i< Vec_length (v); i++ )
    {
      y = v[i]-c;
      t = s+y;
      c = (t-s) - y;
      s = t;
    }  
  return s;
}


int ivec_sum (ivec v)
{
  idx_t i;
  int  s = 0;
  assert (v);
  for (i = 0; i < ivec_length (v); i++) {
    if (v[i] < 0 && s < INT_MIN - v[i])
      it_warning ("underflow in ivec_sum");
    if (v[i] > 0 && s > INT_MAX - v[i])
      it_warning ("overflow in ivec_sum");
    s += v[i];
  }
  return s;
}


cplx cvec_sum (cvec v)
{
  idx_t i;
  cplx s = cplx_0;
  assert (v);
  for (i = 0; i < cvec_length (v); i++)
    s = cadd (s, v[i]);
  return s;
}


vec vec_cum_sum (vec v)
{
  vec  cs = vec_new (vec_length (v));
  idx_t i;
  assert (v);

  if (vec_length (v) > 0)
    cs[0] = v[0];

  for (i = 1; i < Vec_length (v); i++)
    cs[i] = cs[i - 1] + v[i];

  return cs;
}


ivec ivec_cum_sum (ivec v)
{
  ivec cs = ivec_new (ivec_length (v));
  idx_t i;
  assert (v);

  if (ivec_length (v) > 0)
    cs[0] = v[0];

  for (i = 1; i < Vec_length (v); i++) {
    if (v[i] < 0 && cs[i - 1] < INT_MIN - v[i])
      it_warning ("underflow in vec_cum_sum");
    if (v[i] > 0 && cs[i - 1] > INT_MAX - v[i])
      it_warning ("overflow in vec_cum_sum");
    cs[i] = cs[i - 1] + v[i];
  }
  return cs;
}


cvec cvec_cum_sum (cvec v)
{
  cvec cs = cvec_new (cvec_length (v));
  idx_t i;
  assert (v);

  if (cvec_length (v) > 0)
    cs[0] = v[0];

  for (i = 1; i < cvec_length (v); i++) {
    cs[i] = cadd (cs[i - 1], v[i]);
  }
  return cs;
}


double vec_sum_sqr (vec v)
{
  idx_t i;
  double s = 0;
  assert (v);
  for (i = 0; i < Vec_length (v); i++)
    s += v[i] * v[i];
  return s;
}


double vec_sum_between (vec v, idx_t i1, idx_t i2)
{
  idx_t i;
  double s = 0;
  VEC_END_PARAM (v, i2);
  assert (v);
  for (i = i1; i <= i2; i++)
    s += v[i];
  return s;
}


int ivec_sum_between (ivec v, idx_t i1, idx_t i2)
{
  idx_t i;
  int  s = 0;
  VEC_END_PARAM (v, i2);
  assert (v);
  for (i = i1; i <= i2; i++) {
    if (v[i] < 0 && s < INT_MIN - v[i])
      it_warning ("underflow in ivec_sum");
    if (v[i] > 0 && s > INT_MAX - v[i])
      it_warning ("overflow in ivec_sum");
    s += v[i];
  }
  return s;
}


cplx cvec_sum_between (cvec v, idx_t i1, idx_t i2)
{
  idx_t i;
  cplx s = cplx_0;
  VEC_END_PARAM (v, i2);
  assert (v);
  for (i = i1; i <= i2; i++)
    s = cadd (s, v[i]);
  return s;
}


double vec_min (vec v)
{
  idx_t i;
  double m = HUGE_VAL;
  assert (v);
  for (i = 0; i < Vec_length (v); i++)
    if (v[i] < m)
      m = v[i];
  return m;
}


double vec_max (vec v)
{
  idx_t i;
  double m = -HUGE_VAL;
  assert (v);
  for (i = 0; i < Vec_length (v); i++)
    if (v[i] > m)
      m = v[i];
  return m;
}


idx_t vec_min_index (vec v)
{
  idx_t i, mi = NULL_INDEX;
  double m = HUGE_VAL;
  assert (v);
  for (i = 0; i < vec_length (v); i++)
    if (v[i] < m) {
      m = v[i];
      mi = i;
    }
  return mi;
}


idx_t vec_max_index (vec v)
{
  idx_t i, mi = NULL_INDEX;
  double m = -HUGE_VAL;
  assert (v);
  for (i = 0; i < vec_length (v); i++)
    if (v[i] > m) {
      m = v[i];
      mi = i;
    }
  return mi;
}


/* Find the k smallest values using the merge sort */
vec vec_k_min_between (vec v, int k, idx_t a, idx_t b)
{
  vec  r, v1, v2;
  idx_t i1 = 0, i2 = 0;
  int  lmax;

  if (a > b)
    return vec_new_void();

  if (a == b) {
    r = vec_new (1);
    r[0] = v[a];
    return r;
  }

  v1 = vec_k_min_between (v, k, a, (a + b) / 2);
  v2 = vec_k_min_between (v, k, (a + b) / 2 + 1, b);

  lmax = vec_length (v1) + vec_length (v2);
  if (lmax > k)
    lmax = k;

  r = vec_new_alloc (0, lmax);

  /* Proceed with the merge. Note that values are sorted in increasing order */
  while (i1 + i2 < lmax) {
    if (i1 == vec_length (v1))
      vec_push (r, v2[i2++]);
    else if (i2 == vec_length (v2))
      vec_push (r, v1[i1++]);
    else if (v1[i1] < v2[i2])
      vec_push (r, v1[i1++]);
    else
      vec_push (r, v2[i2++]);
  }

  vec_delete (v1);
  vec_delete (v2);
  return r;
}


vec vec_k_min (vec v, int k)
{
  return vec_k_min_between (v, k, 0, vec_length (v) - 1);
}


/* Find the k highest values using the merge sort */
vec vec_k_max_between (vec v, int k, idx_t a, idx_t b)
{
  vec  r, v1, v2;
  idx_t i1 = 0, i2 = 0;
  int  lmax;

  if (a > b)
    return vec_new_void();

  if (a == b) {
    r = vec_new (1);
    r[0] = v[a];
    return r;
  }

  v1 = vec_k_max_between (v, k, a, (a + b) / 2);
  v2 = vec_k_max_between (v, k, (a + b) / 2 + 1, b);

  lmax = vec_length (v1) + vec_length (v2);
  if (lmax > k)
    lmax = k;

  r = vec_new_alloc (0, lmax);

  /* Proceed with the merge. Note that values are sorted in increasing order */
  while (i1 + i2 < lmax) {
    if (i1 >= vec_length (v1))
      vec_push (r, v2[i2++]);
    else if (i2 >= vec_length (v2))
      vec_push (r, v1[i1++]);
    else if (v1[i1] > v2[i2])
      vec_push (r, v1[i1++]);
    else
      vec_push (r, v2[i2++]);
  }
  vec_delete (v1);
  vec_delete (v2);
  return r;
}


vec vec_k_max (vec v, int k)
{
  return vec_k_max_between (v, k, 0, vec_length (v) - 1);
}


/* Find the k smallest values using the merge sort */
ivec ivec_k_min_between (ivec v, int k, idx_t a, idx_t b)
{
  ivec r, v1, v2;
  idx_t i1 = 0, i2 = 0;
  int  lmax;

  if (a > b)
    return ivec_new_void();

  if (a == b) {
    r = ivec_new (1);
    r[0] = v[a];
    return r;
  }

  v1 = ivec_k_min_between (v, k, a, (a + b) / 2);
  v2 = ivec_k_min_between (v, k, (a + b) / 2 + 1, b);

  lmax = ivec_length (v1) + ivec_length (v2);
  if (lmax > k)
    lmax = k;

  r = ivec_new_alloc (0, lmax);

  /* Proceed with the merge. Note that values are sorted in increasing order */
  while (i1 + i2 < lmax) {
    if (i1 == ivec_length (v1))
      ivec_push (r, v2[i2++]);
    else if (i2 == ivec_length (v2))
      ivec_push (r, v1[i1++]);
    else if (v1[i1] < v2[i2])
      ivec_push (r, v1[i1++]);
    else
      ivec_push (r, v2[i2++]);
  }

  ivec_delete (v1);
  ivec_delete (v2);
  return r;
}


ivec ivec_k_min (ivec v, int k)
{
  return ivec_k_min_between (v, k, 0, ivec_length (v) - 1);
}


/* Find the k highest values using the merge sort */
ivec ivec_k_max_between (ivec v, int k, idx_t a, idx_t b)
{
  ivec r, v1, v2;
  idx_t i1 = 0, i2 = 0;
  int  lmax;

  if (a > b)
    return ivec_new_void();

  if (a == b) {
    r = ivec_new (1);
    r[0] = v[a];
    return r;
  }

  v1 = ivec_k_max_between (v, k, a, (a + b) / 2);
  v2 = ivec_k_max_between (v, k, (a + b) / 2 + 1, b);

  lmax = k;
  if (lmax > k)
    lmax = k;

  r = ivec_new_alloc (0, lmax);

  /* Proceed with the merge. Note that values are sorted in increasing order */
  while (i1 + i2 < lmax) {
    if (i1 >= ivec_length (v1))
      ivec_push (r, v2[i2++]);
    else if (i2 >= ivec_length (v2))
      ivec_push (r, v1[i1++]);
    else if (v1[i1] > v2[i2])
      ivec_push (r, v1[i1++]);
    else
      ivec_push (r, v2[i2++]);
  }
  ivec_delete (v1);
  ivec_delete (v2);
  return r;
}


ivec ivec_k_max (ivec v, int k)
{
  return ivec_k_max_between (v, k, 0, ivec_length (v) - 1);
}


/* Find the k smallest values using the merge sort */
ivec vec_k_min_index_between (vec v, int k, idx_t a, idx_t b)
{
  ivec r, v1, v2;
  idx_t i1 = 0, i2 = 0;
  int  lmax;

  if (a > b)
    return ivec_new_void();

  if (a == b) {
    r = ivec_new (1);
    r[0] = a;
    return r;
  }

  v1 = vec_k_min_index_between (v, k, a, (a + b) / 2);
  v2 = vec_k_min_index_between (v, k, (a + b) / 2 + 1, b);

  lmax = ivec_length (v1) + ivec_length (v2);
  if (lmax > k)
    lmax = k;

  r = ivec_new_alloc (0, lmax);

  /* Proceed with the merge. Note that values are sorted in increasing order */
  while (i1 + i2 < lmax) {
    if (i1 == ivec_length (v1))
      ivec_push (r, v2[i2++]);
    else if (i2 == ivec_length (v2))
      ivec_push (r, v1[i1++]);
    else if (v[v1[i1]] < v[v2[i2]])
      ivec_push (r, v1[i1++]);
    else
      ivec_push (r, v2[i2++]);
  }

  ivec_delete (v1);
  ivec_delete (v2);
  return r;
}


ivec vec_k_min_index (vec v, int k)
{
  return vec_k_min_index_between (v, k, 0, vec_length (v) - 1);
}


/* Find the k smallest values using the merge sort */
ivec vec_k_max_index_between (vec v, int k, idx_t a, idx_t b)
{
  ivec r, v1, v2;
  idx_t i1 = 0, i2 = 0;
  int  lmax;

  if (a > b)
    return ivec_new_void();

  if (a == b) {
    r = ivec_new (1);
    r[0] = a;
    return r;
  }

  v1 = vec_k_max_index_between (v, k, a, (a + b) / 2);
  v2 = vec_k_max_index_between (v, k, (a + b) / 2 + 1, b);

  lmax = ivec_length (v1) + ivec_length (v2);
  if (lmax > k)
    lmax = k;

  r = ivec_new_alloc (0, lmax);

  /* Proceed with the merge. Note that values are sorted in increasing order */
  while (i1 + i2 < lmax) {
    if (i1 == ivec_length (v1))
      ivec_push (r, v2[i2++]);
    else if (i2 == ivec_length (v2))
      ivec_push (r, v1[i1++]);
    else if (v[v1[i1]] > v[v2[i2]])
      ivec_push (r, v1[i1++]);
    else
      ivec_push (r, v2[i2++]);
  }

  ivec_delete (v1);
  ivec_delete (v2);
  return r;
}


ivec vec_k_max_index (vec v, int k)
{
  return vec_k_max_index_between (v, k, 0, vec_length (v) - 1);
}


/* Find the k smallest values using the merge sort */
ivec ivec_k_min_index_between (ivec v, int k, idx_t a, idx_t b)
{
  ivec r, v1, v2;
  idx_t i1 = 0, i2 = 0;
  int  lmax;

  if (a > b)
    return ivec_new_void();

  if (a == b) {
    r = ivec_new (1);
    r[0] = a;
    return r;
  }

  v1 = ivec_k_min_index_between (v, k, a, (a + b) / 2);
  v2 = ivec_k_min_index_between (v, k, (a + b) / 2 + 1, b);

  lmax = ivec_length (v1) + ivec_length (v2);
  if (lmax > k)
    lmax = k;

  r = ivec_new_alloc (0, lmax);

  /* Proceed with the merge. Note that values are sorted in increasing order */
  while (i1 + i2 < lmax) {
    if (i1 == ivec_length (v1))
      ivec_push (r, v2[i2++]);
    else if (i2 == ivec_length (v2))
      ivec_push (r, v1[i1++]);
    else if (v[v1[i1]] < v[v2[i2]])
      ivec_push (r, v1[i1++]);
    else
      ivec_push (r, v2[i2++]);
  }

  ivec_delete (v1);
  ivec_delete (v2);
  return r;
}


ivec ivec_k_min_index (ivec v, int k)
{
  return ivec_k_min_index_between (v, k, 0, ivec_length (v) - 1);
}


/* Find the k smallest values using the merge sort */
ivec ivec_k_max_index_between (ivec v, int k, idx_t a, idx_t b)
{
  ivec r, v1, v2;
  idx_t i1 = 0, i2 = 0;
  int  lmax;

  if (a > b)
    return ivec_new_void();

  if (a == b) {
    r = ivec_new (1);
    r[0] = a;
    return r;
  }

  v1 = ivec_k_max_index_between (v, k, a, (a + b) / 2);
  v2 = ivec_k_max_index_between (v, k, (a + b) / 2 + 1, b);

  lmax = ivec_length (v1) + ivec_length (v2);
  if (lmax > k)
    lmax = k;

  r = ivec_new_alloc (0, lmax);

  /* Proceed with the merge. Note that values are sorted in increasing order */
  while (i1 + i2 < lmax) {
    if (i1 == ivec_length (v1))
      ivec_push (r, v2[i2++]);
    else if (i2 == ivec_length (v2))
      ivec_push (r, v1[i1++]);
    else if (v[v1[i1]] > v[v2[i2]])
      ivec_push (r, v1[i1++]);
    else
      ivec_push (r, v2[i2++]);
  }

  ivec_delete (v1);
  ivec_delete (v2);
  return r;
}


ivec ivec_k_max_index (ivec v, int k)
{
  return ivec_k_max_index_between (v, k, 0, ivec_length (v) - 1);
}


double vec_mean (vec v)
{
  assert (v);
  return vec_sum (v) / Vec_length (v);
}

double vec_mean_robust (vec v)
{
  assert (v); 
  return vec_sum_robust (v) / Vec_length (v);
}


int ivec_median (ivec v)
{
  ivec c;
  int  m;
  assert (v);
  if (ivec_length (v) == 0) {
    it_warning ("Undefined median value for vector of size 0. Return 0.\n");
    return 0;
  }

  if (ivec_length (v) == 1)
    return v[0];

  c = ivec_clone (v);
  ivec_sort (c);
  
  m = (c[(ivec_length (v) - 1) / 2] + c[ivec_length (v) / 2]) / 2;

  ivec_delete (c);
  return m;
}


double vec_median (vec v)
{
  vec  c;
  double m;
  assert (v);
  if (vec_length (v) == 0) {
    it_warning ("Undefined median value for vector of size 0. Return 0.\n");
    return 0;
  }

  if (vec_length (v) == 1)
    return v[0];

  c = vec_clone (v);
  vec_sort (c);

  m = (c[(vec_length (v) - 1) / 2] + c[vec_length (v) / 2]) / 2;

  vec_delete (c);
  return m;
}


double vec_variance (vec v)
{
  idx_t i, l;
  double sum = 0;
  double var = 0;
  assert (v);
  l = Vec_length (v);
  assert (l > 1);		/* otherwise the unbiased variance is not defined */
  for (i = 0; i < l; i++) {
    sum += v[i];
    var += v[i] * v[i];
  }

  return (var - sum * sum / l) / (l - 1);
}

/* Knuth's two-pass algorithm for damn-robust variance computation */
double vec_variance_robust (vec v) 
{
  double s= 0., r, c, m= 0.; 
  idx_t i;
  size_t l;
  
  assert (v); 
  l = Vec_length( v );
  /* First pass: compute mean with Kahan's method. */
  r = vec_mean_robust (v); 
  /* r is the mean estimate: use residuals for Knuth's algorithm. */ 
  for ( i= 0; i< l; i++ )
    {
      c = ( v[i]-r ) - m;
      m += c/(double)(i+1);
      s+= c*( ( v[i]-r ) - m );
    }  
  return (s/(double)(l-1) );
}

double vec_cov (vec v1, vec v2)
{
  double n = vec_length (v1);
  double m1 = vec_mean (v1);
  double m2 = vec_mean (v2);
  double cov;
  vec v1_s, v2_s;

  vec v1_mean = vec_new_set (m1, n);
  vec v2_mean = vec_new_set (m2, n);


  v1_s = vec_clone (v1);
  v2_s = vec_clone (v2);

  vec_sub (v1_s, v1_mean);
  vec_sub (v2_s, v2_mean);



  vec_mul (v1_s, v2_s);


  cov = (1 / n) * vec_sum (v1_s);

  vec_delete (v1_s);
  vec_delete (v2_s);
  vec_delete (v1_mean);
  vec_delete (v2_mean);

  return cov;

}

double vec_norm (vec v, double n)
{
  idx_t i;
  double nr = 0;
  assert (v);
  assert (n > 0);
  for (i = 0; i < vec_length (v); i++)
    nr += pow (fabs (v[i]), n);

  return pow (nr, 1.0 / n);
}


/* integer norm is not optimized at all */
double ivec_norm (ivec v, int n)
{
  idx_t i;
  double nr = 0;
  assert (v);
  assert (n > 0);
  for (i = 0; i < ivec_length (v); i++)
    nr += pow (fabs ((double) v[i]), (double) n);

  return pow ((double) nr, 1.0 / (double) n);
}


/* General function                                                             */
void vec_apply_function (vec v, it_function_t function, it_args_t args)
{
  idx_t i;
  idx_t l;
  assert (v);
  l = vec_length (v);
  for (i = 0; i < l; i++)
    v[i] = function (v[i], args);
}


vec vec_new_apply_function (vec v, it_function_t function, it_args_t args)
{
  vec  r;
  idx_t i;
  idx_t l;
  assert (v);
  l = vec_length (v);
  r = vec_new (l);
  for (i = 0; i < l; i++)
    r[i] = function (v[i], args);
  return (r);
}


ivec ivec_apply_function (ivec v, it_ifunction_t function, it_args_t args)
{
  idx_t i;
  idx_t l;
  assert (v);
  l = ivec_length (v);
  for (i = 0; i < l; i++)
    v[i] = function (v[i], args);
  return v;
}


ivec ivec_new_apply_function (ivec v, it_ifunction_t function, it_args_t args)
{
  ivec r;
  idx_t i;
  idx_t l;
  assert (v);
  l = ivec_length (v);
  r = ivec_new (l);
  for (i = 0; i < l; i++)
    r[i] = function (v[i], args);
  return (r);
}


/*------------------------------------------------------------------------------*/
void vec_reverse (vec v)
{
  idx_t i, j, m;
  double tmp;
  assert (v);

  m = vec_length (v) / 2;
  for (i = 0, j = vec_length (v) - 1; i < m; i++, j--) {
    tmp = v[i];
    v[i] = v[j];
    v[j] = tmp;
  }
}


void ivec_reverse (ivec v)
{
  idx_t i, j, m;
  int  tmp;
  assert (v);

  m = ivec_length (v) / 2;
  for (i = 0, j = ivec_length (v) - 1; i < m; i++, j--) {
    tmp = v[i];
    v[i] = v[j];
    v[j] = tmp;
  }
}


void bvec_reverse (bvec v)
{
  idx_t i, j, m;
  byte tmp;
  assert (v);

  m = bvec_length (v) / 2;
  for (i = 0, j = bvec_length (v) - 1; i < m; i++, j--) {
    tmp = v[i];
    v[i] = v[j];
    v[j] = tmp;
  }
}


void cvec_reverse (cvec v)
{
  idx_t i, j, m;
  double tmp_real, tmp_imag;
  assert (v);

  m = cvec_length (v) / 2;
  for (i = 0, j = cvec_length (v) - 1; i < m; i++, j--) {
    tmp_real = creal (v[i]);
    tmp_imag = cimag (v[i]);
    creal (v[i]) = creal (v[j]);
    cimag (v[i]) = cimag (v[j]);
    creal (v[j]) = tmp_real;
    creal (v[j]) = tmp_imag;
  }
}


/*------------------------------------------------------------------------------*/
vec vec_new_reverse (vec v)
{
  idx_t i;
  vec  cl;
  assert (v);

  cl = vec_new (vec_length (v));
  if (cl == NULL)
    return NULL;

  for (i = 0; i < vec_length (v); i++)
    cl[i] = v[vec_length (v) - 1 - i];

  return cl;
}

ivec ivec_new_reverse (ivec v)
{
  idx_t i;
  ivec cl;
  assert (v);

  cl = ivec_new (ivec_length (v));
  if (cl == NULL)
    return NULL;

  for (i = 0; i < ivec_length (v); i++)
    cl[i] = v[ivec_length (v) - 1 - i];

  return cl;
}

bvec bvec_new_reverse (bvec v)
{
  idx_t i;
  bvec cl;
  assert (v);

  cl = bvec_new (bvec_length (v));
  if (cl == NULL)
    return NULL;

  for (i = 0; i < bvec_length (v); i++)
    cl[i] = v[bvec_length (v) - 1 - i];

  return cl;
}


cvec cvec_new_reverse (cvec v)
{
  idx_t i;
  cvec cl;
  assert (v);

  cl = cvec_new (cvec_length (v));
  if (cl == NULL)
    return NULL;

  for (i = 0; i < cvec_length (v); i++) {
    creal (cl[i]) = creal (v[cvec_length (v) - 1 - i]);
    cimag (cl[i]) = cimag (v[cvec_length (v) - 1 - i]);
  }

  return cl;
}


/*-----------------------------------------------------------------*/
/* Return the first position where the value a occurs              */
idx_t vec_find_first (vec v, double a)
{
  idx_t i;
  assert (v);
  for (i = 0; i < Vec_length (v); i++)
    if (v[i] == a)
      return i;
  return NULL_INDEX;
}


idx_t ivec_find_first (ivec v, int a)
{
  idx_t i;
  assert (v);
  for (i = 0; i < ivec_length (v); i++)
    if (v[i] == a)
      return i;
  return NULL_INDEX;
}


idx_t bvec_find_first (bvec v, byte a)
{
  idx_t i;
  assert (v);
  for (i = 0; i < bvec_length (v); i++)
    if (v[i] == a)
      return i;
  return NULL_INDEX;
}


idx_t cvec_find_first (cvec v, cplx a)
{
  idx_t i;
  assert (v);
  for (i = 0; i < cvec_length (v); i++)
    if (ceq (v[i], a))
      return i;
  return NULL_INDEX;
}


/*-----------------------------------------------------------------*/
ivec vec_find (vec v, double a)
{
  idx_t count = vec_count (v, a);
  ivec pos = ivec_new (count);
  idx_t i, j;
  assert (v);

  for (i = 0, j = 0; i < Vec_length (v); i++)
    if (v[i] == a)
      pos[j++] = i;

  return pos;
}


ivec ivec_find (ivec v, int a)
{
  idx_t count = ivec_count (v, a);
  ivec pos = ivec_new (count);
  idx_t i, j;
  assert (v);

  for (i = 0, j = 0; i < ivec_length (v); i++)
    if (v[i] == a)
      pos[j++] = i;

  return pos;
}


ivec bvec_find (bvec v, byte a)
{
  idx_t count = bvec_count (v, a);
  ivec pos = ivec_new (count);
  idx_t i, j;
  assert (v);

  for (i = 0, j = 0; i < bvec_length (v); i++)
    if (v[i] == a)
      pos[j++] = i;

  return pos;
}


ivec cvec_find (cvec v, cplx a)
{
  idx_t count = cvec_count (v, a);
  ivec pos = ivec_new (count);
  idx_t i, j;
  assert (v);

  for (i = 0, j = 0; i < cvec_length (v); i++)
    if (ceq (v[i], a))
      pos[j++] = i;

  return pos;
}


/* Assuming that the vector is sorted, return the position where it can be found
   or return the position where is should be inserted to maintain the 
   vector sorted                                                                */
idx_t vec_find_sorted (vec v, double a)
{
  idx_t left = 0;
  idx_t right = vec_length (v) - 1;
  idx_t m;

  if (a > v[right] )
    return vec_length (v);

  while (left < right - 1) {
    m = ( left + right ) / 2;
    
    if (a < v[m])
      right = m;
    else
      left = m;
  }
  
  if (a <= v[left])
    return left;      
  else
    return right;      
}


idx_t ivec_find_sorted (ivec v, int a)
{
  idx_t left = 0;
  idx_t right = ivec_length (v) - 1;
  idx_t m;

  if (a > v[right] )
    return ivec_length (v);

  while (left < right - 1) {
    m = ( left + right ) / 2;
    
    if (a < v[m])
      right = m;
    else
      left = m;
  }
  
  if (a <= v[left])
    return left;      
  else
    return right;
}


idx_t bvec_find_sorted (bvec v, byte a)
{
  idx_t left = 0;
  idx_t right = bvec_length (v) - 1;
  idx_t m;

  if (a > v[right] )
    return bvec_length (v);

  while (left < right - 1) {
    m = ( left + right ) / 2;
    
    if (a < v[m])
      right = m;
    else
      left = m;
  }
  
  if (a <= v[left])
    return left;      
  else
    return right;      
}



/*-----------------------------------------------------------------*/
ivec vec_replace (vec v, double a, double b)
{
  idx_t count = vec_count (v, a);
  ivec pos = ivec_new (count);
  idx_t i, j;
  assert (v);

  for (i = 0, j = 0; i < Vec_length (v); i++)
    if (v[i] == a) {
      pos[j++] = i;
      v[i] = b;
    }

  return pos;
}


ivec ivec_replace (ivec v, int a, int b)
{
  idx_t count = ivec_count (v, a);
  ivec pos = ivec_new (count);
  idx_t i, j;
  assert (v);

  for (i = 0, j = 0; i < ivec_length (v); i++)
    if (v[i] == a) {
      pos[j++] = i;
      v[i] = b;
    }

  return pos;
}


ivec bvec_replace (bvec v, byte a, byte b)
{
  idx_t count = bvec_count (v, a);
  ivec pos = ivec_new (count);
  idx_t i, j;
  assert (v);

  for (i = 0, j = 0; i < bvec_length (v); i++)
    if (v[i] == a) {
      pos[j++] = i;
      v[i] = b;
    }

  return pos;
}


ivec cvec_replace (cvec v, cplx a, cplx b)
{
  idx_t count = cvec_count (v, a);
  ivec pos = ivec_new (count);
  idx_t i, j;
  assert (v);

  for (i = 0, j = 0; i < cvec_length (v); i++)
    if (ceq (v[i], a)) {
      pos[j++] = i;
      v[i] = b;
    }

  return pos;
}


/*-----------------------------------------------------------------*/
idx_t vec_count (vec v, double a)
{
  idx_t c = 0;
  idx_t i;
  assert (v);
  for (i = 0; i < Vec_length (v); i++)
    if (v[i] == a)
      c++;
  return c;
}


idx_t bvec_count (bvec v, byte a)
{
  idx_t c = 0;
  idx_t i;
  assert (v);
  for (i = 0; i < bvec_length (v); i++)
    if (v[i] == a)
      c++;
  return c;
}


idx_t ivec_count (ivec v, int a)
{
  idx_t c = 0;
  idx_t i;
  assert (v);
  for (i = 0; i < ivec_length (v); i++)
    if (v[i] == a)
      c++;
  return c;
}


idx_t cvec_count (cvec v, cplx a)
{
  idx_t c = 0;
  idx_t i;
  assert (v);
  for (i = 0; i < cvec_length (v); i++)
    if (ceq (v[i], a))
      c++;
  return c;
}


/*----------------------------------------------------------------*/
/* Return the set of positions of value a                                       */
/*----------------------------------------------------------------*/

vec vec_new_concat (vec v1, vec v2)
{
  idx_t i, j;
  vec  v;
  assert (v1);
  assert (v2);
  v = vec_new (vec_length (v1) + vec_length (v2));
  for (i = 0; i < vec_length (v1); i++)
    v[i] = v1[i];

  for (j = 0; j < vec_length (v2); j++, i++)
    v[i] = v2[j];
  return v;
}


ivec ivec_new_concat (ivec v1, ivec v2)
{
  idx_t i, j;
  ivec v;
  assert (v1);
  assert (v2);
  v = ivec_new (ivec_length (v1) + ivec_length (v2));
  for (i = 0; i < ivec_length (v1); i++)
    v[i] = v1[i];

  for (j = 0; j < ivec_length (v2); j++, i++)
    v[i] = v2[j];
  return v;
}


bvec bvec_new_concat (bvec v1, bvec v2)
{
  idx_t i, j;
  bvec v;
  assert (v1);
  assert (v2);
  v = bvec_new (bvec_length (v1) + bvec_length (v2));
  for (i = 0; i < bvec_length (v1); i++)
    v[i] = v1[i];

  for (j = 0; j < bvec_length (v2); j++, i++)
    v[i] = v2[j];
  return v;
}


cvec cvec_new_concat (cvec v1, cvec v2)
{
  idx_t i, j;
  cvec v;
  assert (v1);
  assert (v2);
  v = cvec_new (cvec_length (v1) + cvec_length (v2));
  for (i = 0; i < cvec_length (v1); i++)
    v[i] = v1[i];

  for (j = 0; j < cvec_length (v2); j++, i++)
    v[i] = v2[j];
  return v;
}


/*----------------------------------------------------------------*/

vec vec_new_unique (vec v)
{
  idx_t i;
  vec  vsorted, vtmp = vec_new (0);

  assert (v);
  if (vec_length (v) == 0)
    return vtmp;

  vsorted = vec_clone (v);
  vec_sort (vsorted);

  vec_push (vtmp, vsorted[0]);
  for (i = 1; i < vec_length (v); i++)
    if (vsorted[i] != vsorted[i - 1])
      vec_push (vtmp, vsorted[i]);

  vec_delete (vsorted);
  return vtmp;
}


ivec ivec_new_unique (ivec v)
{
  idx_t i;
  ivec vsorted, vtmp = ivec_new (0);

  assert (v);
  if (ivec_length (v) == 0)
    return vtmp;

  vsorted = ivec_clone (v);
  ivec_sort (vsorted);

  ivec_push (vtmp, vsorted[0]);
  for (i = 1; i < ivec_length (v); i++)
    if (vsorted[i] != vsorted[i - 1])
      ivec_push (vtmp, vsorted[i]);

  ivec_delete (vsorted);
  return vtmp;
}


bvec bvec_new_unique (bvec v)
{
  idx_t i;
  bvec vsorted = bvec_new (0);	/* The vector that will be generated */
  ivec vtmp = ivec_new_zeros (1 << (8 * sizeof (byte)));

  assert (v);

  for (i = 0; i < bvec_length (v); i++)
    vtmp[v[i]]++;

  for (i = 0; i < ivec_length (vtmp); i++)
    if (vtmp[i] > 0)
      bvec_push (vsorted, (byte) i);

  ivec_delete (vtmp);
  return vsorted;
}


vec vec_new_union (vec v1, vec v2)
{
  idx_t i1, i2;
  vec  vu1 = vec_new_unique (v1);
  vec  vu2 = vec_new_unique (v2);
  vec  vu = vec_new_alloc (0, vec_length (vu1) + vec_length (vu2));

  for (i1 = 0, i2 = 0; i1 < vec_length (vu1) && i2 < vec_length (vu2);)
    if (vu1[i1] < vu2[i2])
      vec_push (vu, vu1[i1++]);
    else if (vu1[i1] > vu2[i2])
      vec_push (vu, vu2[i2++]);
    else {
      vec_push (vu, vu1[i1]);
      i1++;
      i2++;
    };

  /* Put the remaining elements */
  while (i1 < vec_length (vu1))
    vec_push (vu, vu1[i1++]);

  while (i2 < vec_length (vu2))
    vec_push (vu, vu2[i2++]);

  vec_delete (vu1);
  vec_delete (vu2);
  return vu;
}


ivec ivec_new_union (ivec v1, ivec v2)
{
  idx_t i1, i2;
  ivec vu1 = ivec_new_unique (v1);
  ivec vu2 = ivec_new_unique (v2);
  ivec vu = ivec_new_alloc (0, ivec_length (vu1) + ivec_length (vu2));

  for (i1 = 0, i2 = 0; i1 < ivec_length (vu1) && i2 < ivec_length (vu2);)
    if (vu1[i1] < vu2[i2])
      ivec_push (vu, vu1[i1++]);
    else if (vu1[i1] > vu2[i2])
      ivec_push (vu, vu2[i2++]);
    else {
      ivec_push (vu, vu1[i1]);
      i1++;
      i2++;
    };

  /* Put the remaining elements */
  while (i1 < ivec_length (vu1))
    ivec_push (vu, vu1[i1++]);

  while (i2 < ivec_length (vu2))
    ivec_push (vu, vu2[i2++]);

  ivec_delete (vu1);
  ivec_delete (vu2);
  return vu;
}


bvec bvec_new_union (bvec v1, bvec v2)
{
  idx_t i;
  bvec vsorted = bvec_new (0);	/* The vector that will be generated */
  ivec vtmp = ivec_new_zeros (1 << (8 * sizeof (byte)));

  assert (v1);
  assert (v2);

  for (i = 0; i < bvec_length (v1); i++)
    vtmp[v1[i]]++;

  for (i = 0; i < bvec_length (v2); i++)
    vtmp[v2[i]]++;

  for (i = 0; i < ivec_length (vtmp); i++)
    if (vtmp[i] > 0)
      bvec_push (vsorted, (byte) i);

  ivec_delete (vtmp);
  return vsorted;
}


vec vec_new_intersection (vec v1, vec v2)
{
  idx_t i1, i2;
  vec  vu1 = vec_new_unique (v1);
  vec  vu2 = vec_new_unique (v2);
  vec  vu = vec_new_alloc (0, vec_length (vu1) > vec_length (vu2) ?
			   vec_length (vu2) : vec_length (vu1));

  for (i1 = 0, i2 = 0; i1 < vec_length (vu1) && i2 < vec_length (vu2);
       i1++, i2++) {
    while (vu1[i1] < vu2[i2] && i1 < vec_length (vu1))
      i1++;

    if (i1 == vec_length (vu1))
      break;

    while (vu2[i2] < vu1[i1] && i2 < vec_length (vu2))
      i2++;

    if (i2 == vec_length (vu2))
      break;

    vec_push (vu, vu1[i1]);
  }

  vec_delete (vu1);
  vec_delete (vu2);
  return vu;
}


ivec ivec_new_intersection (ivec v1, ivec v2)
{
  idx_t i1, i2;
  ivec vu1 = ivec_new_unique (v1);
  ivec vu2 = ivec_new_unique (v2);
  ivec vu = ivec_new_alloc (0, ivec_length (vu1) > ivec_length (vu2) ?
			    ivec_length (vu2) : ivec_length (vu1));

  for (i1 = 0, i2 = 0; i1 < ivec_length (vu1) && i2 < ivec_length (vu2);
       i1++, i2++) {
    while (vu1[i1] < vu2[i2] && i1 < ivec_length (vu1))
      i1++;

    if (i1 == ivec_length (vu1))
      break;

    while (vu2[i2] < vu1[i1] && i2 < ivec_length (vu2))
      i2++;

    if (i2 == ivec_length (vu2))
      break;

    ivec_push (vu, vu1[i1]);
  }

  ivec_delete (vu1);
  ivec_delete (vu2);
  return vu;
}


bvec bvec_new_intersection (bvec v1, bvec v2)
{
  idx_t i;
  bvec vsorted = bvec_new (0);	/* The vector that will be generated */
  ivec vtmp = ivec_new_zeros (1 << (8 * sizeof (byte)));

  assert (v1);
  assert (v2);

  for (i = 0; i < bvec_length (v1); i++)
    vtmp[v1[i]]++;

  for (i = 0; i < bvec_length (v2); i++)
    if (vtmp[v2[i]] > 0)
      vtmp[v2[i]] = -1;

  for (i = 0; i < ivec_length (vtmp); i++)
    if (vtmp[i] == -1)
      bvec_push (vsorted, (byte) i);

  ivec_delete (vtmp);
  return vsorted;
}


/*----------------------------------------------------------------*/
/* Return the vector composed of the elements of v indexed by idx */
vec vec_index_by (vec v, ivec idx)
{
  vec  r = vec_new (ivec_length (idx));
  idx_t i;
  assert (v);
  assert (idx);
  for (i = 0; i < ivec_length (idx); i++)
    r[i] = v[idx[i]];
  return r;
}


ivec ivec_index_by (ivec v, ivec idx)
{
  ivec r = ivec_new (ivec_length (idx));
  idx_t i;
  assert (v);
  assert (idx);
  for (i = 0; i < ivec_length (idx); i++)
    r[i] = v[idx[i]];
  return r;
}


bvec bvec_index_by (bvec v, ivec idx)
{
  bvec r = bvec_new (ivec_length (idx));
  idx_t i;
  assert (v);
  assert (idx);
  for (i = 0; i < ivec_length (idx); i++)
    r[i] = v[idx[i]];
  return r;
}


cvec cvec_index_by (cvec v, ivec idx)
{
  cvec r = cvec_new (ivec_length (idx));
  idx_t i;
  assert (v);
  assert (idx);
  for (i = 0; i < ivec_length (idx); i++)
    r[i] = v[idx[i]];
  return r;
}


/*----------------------------------------------------------------*/
/* Sorting                                                        */

void Vec_sort (Vec v, int (*elem_leq) (const void *, const void *))
{
  assert (v);
  qsort (v, Vec_length (v), Vec_element_size (v), elem_leq);
}

ivec Vec_sort_index (Vec v, int (*elem_leq_idx) (const void *, const void *))
{
  /* The trick here is to create a vector of pointers to the elements
     to sort and view it jointly as a vector of integers of the form
     index * sizeof(element_type) + offset. The sort is done on the
     pointed elements, actually sorting the pointers. By doing some
     pointer arithmetic to remove the offset and divide properly by
     the element size, the indexes are retrieved.
     On platforms where sizeof(void *) > sizeof(int) it will return
     a pointer with unused allocated memory at the end.
  */
  int  i;
  size_t elem_size = Vec_element_size (v);
  void **ptr = Vec_new (void *, Vec_length (v));
  ivec idx = (ivec) ptr;

  /* create a vector of pointers to each element of v */
  for (i = 0; i < Vec_length (v); i++)
    ptr[i] = (void *) ((char *) v + i * elem_size);

  /* sort it */
  qsort (ptr, Vec_length (ptr), sizeof (void *), elem_leq_idx);

  /* now subtract the address of v from each element */
  for (i = 0; i < Vec_length (v); i++)
    idx[i] = (int) ((char *) ptr[i] - (char *) v) / elem_size;

  /* correct the vector size (ultimately ugly) */
  if (sizeof (void *) != sizeof (int)) {
    Vec_element_size (idx) = sizeof (int);
    Vec_length_max (idx) *= sizeof (void *) / sizeof (int);
  }
  return (idx);
}


/*----------------------------------------------------------------*/
static int double_leq (const void *d1, const void *d2)
{
  if (*((double *) d1) > *((double *) d2))
    return 1;
  else if (*((double *) d1) == *((double *) d2))
    return 0;
  return -1;
}


static int double_leq_index (const void *d1, const void *d2)
{
  return (double_leq (*(void **) d1, *(void **) d2));
}

void vec_sort (vec v)
{
  Vec_sort (v, double_leq);
}


ivec vec_sort_index (vec v)
{
  return (Vec_sort_index (v, double_leq_index));
}


/*----------------------------------------------------------------*/
static int int_leq (const void *d1, const void *d2)
{
  if (*((int *) d1) > *((int *) d2))
    return 1;
  else if (*((int *) d1) == *((int *) d2))
    return 0;
  return -1;
}


static int int_leq_index (const void *d1, const void *d2)
{
  return (int_leq (*(void **) d1, *(void **) d2));
}

void ivec_sort (ivec v)
{
  Vec_sort (v, int_leq);
}


ivec ivec_sort_index (ivec v)
{
  return (Vec_sort_index (v, int_leq_index));
}


/*----------------------------------------------------------------*/
static int byte_leq (const void *d1, const void *d2)
{
  if (*((byte *) d1) > *((byte *) d2))
    return 1;
  else if (*((byte *) d1) == *((byte *) d2))
    return 0;
  return -1;
}


static int byte_leq_index (const void *d1, const void *d2)
{
  return (byte_leq (*(void **) d1, *(void **) d2));
}


void bvec_sort (bvec v)
{
  Vec_sort (v, byte_leq);
}


ivec bvec_sort_index (bvec v)
{
  return (Vec_sort_index (v, byte_leq_index));
}


/*---------------------------------------------------------------------------*/
/*                Special Vectors                                            */
/*---------------------------------------------------------------------------*/

/* The following functions proceeds with the modification of already allocated
   vector.                                                                   */

void vec_void (vec v)
{
  Vec_void (v);
}


void ivec_void (ivec v)
{
  Vec_void (v);
}


void bvec_void (bvec v)
{
  Vec_void (v);
}


void cvec_void (cvec v)
{
  Vec_void (v);
}


/*------------------------------------------------------------------------------*/
void vec_zeros (vec v)
{
  Vec_set (v, 0);
}


void ivec_zeros (ivec v)
{
  Vec_set (v, 0);
}


void bvec_zeros (bvec v)
{
  Vec_set (v, 0);
}


void cvec_zeros (cvec v)
{
  Vec_set (v, cplx_0);
}


vec vec_new_zeros (idx_t N) {
  return (vec_new_set (0., N));
}


ivec ivec_new_zeros (idx_t N) {
  return (ivec_new_set (0, N));
}


bvec bvec_new_zeros (idx_t N) {
  return (bvec_new_set (0, N));
}


cvec cvec_new_zeros (idx_t N) {
  return (cvec_new_set (cplx_0, N));
}


/*------------------------------------------------------------------------------*/
void vec_ones (vec v)
{
  vec_set (v, 1);
}


void ivec_ones (ivec v)
{
  ivec_set (v, 1);
}


void bvec_ones (bvec v)
{
  bvec_set (v, 1);
}


void cvec_ones (cvec v)
{
  cvec_set (v, cplx_1);
}


vec vec_new_ones (idx_t N) 
{
  return (vec_new_set (1., N));
}


ivec ivec_new_ones (idx_t N) 
{
  return (ivec_new_set (1, N));
}


bvec bvec_new_ones (idx_t N) 
{
  return (bvec_new_set (1, N));
}


cvec cvec_new_ones (idx_t N) 
{
  return (cvec_new_set (cplx_1, N));
}

/*------------------------------------------------------------------------------*/
void vec_range (vec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < vec_length (v); i++)
    v[i] = i;
}


void ivec_range (ivec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < ivec_length (v); i++)
    v[i] = i;
}


void bvec_range (bvec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < bvec_length (v); i++)
    v[i] = i;
}


void cvec_range (cvec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < cvec_length (v); i++) {
    v[i].r = i;
    v[i].i = 0;
  }
}


/*------------------------------------------------------------------------------*/
void vec_1N (vec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < vec_length (v); i++)
    v[i] = i + 1;
}


void ivec_1N (ivec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < ivec_length (v); i++)
    v[i] = i + 1;
}


void bvec_1N (bvec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < bvec_length (v); i++)
    v[i] = i + 1;
}


void cvec_1N (cvec v)
{
  idx_t i;
  assert (v);
  for (i = 0; i < cvec_length (v); i++) {
    v[i].r = i + 1;
    v[i].i = 0;
  }
}


/*------------------------------------------------------------------------------*/
void vec_arithm (vec v, double start, double incr)
{
  idx_t i;
  double value = start;
  assert (v);
  for (i = 0; i < vec_length (v); i++, value += incr)
    v[i] = value;
}


void ivec_arithm (ivec v, int start, int incr)
{
  idx_t i;
  int  value = start;
  assert (v);
  for (i = 0; i < ivec_length (v); i++, value += incr)
    v[i] = value;
}


void bvec_arithm (bvec v, byte start, byte incr)
{
  idx_t i;
  byte value = start;
  assert (v);
  for (i = 0; i < bvec_length (v); i++, value += incr)
    v[i] = value;
}


void cvec_arithm (cvec v, cplx start, cplx incr)
{
  idx_t i;
  cplx value = start;
  assert (v);
  for (i = 0; i < cvec_length (v); i++, value = cadd (value, incr))
    v[i] = value;
}


/*------------------------------------------------------------------------------*/
void vec_geom (vec v, double start, double r)
{
  idx_t i;
  double value = start;
  assert (v);
  for (i = 0; i < vec_length (v); i++, value *= r)
    v[i] = value;
}


void ivec_geom (ivec v, int start, int r)
{
  idx_t i;
  int  value = start;
  assert (v);
  for (i = 0; i < ivec_length (v); i++, value *= r)
    v[i] = value;
}


void bvec_geom (bvec v, byte start, byte r)
{
  idx_t i;
  byte value = start;
  assert (v);
  for (i = 0; i < bvec_length (v); i++, value *= r)
    v[i] = value;
}


void cvec_geom (cvec v, cplx start, cplx r)
{
  idx_t i;
  cplx value = start;
  assert (v);
  for (i = 0; i < cvec_length (v); i++, value = cmul (value, r))
    v[i] = value;
}


/*------------------------------------------------------------------------------*/
/* Note: the functions returning a vector pointer allocate memory that must     */
/* be free afterwards                                                           */

vec vec_new_void ()
{
  return vec_new (0);
}


ivec ivec_new_void ()
{
  return ivec_new (0);
}


bvec bvec_new_void ()
{
  return bvec_new (0);
}

cvec cvec_new_void ()
{
  return cvec_new (0);
}


/*---------------------------------------------------------------------------*/
vec vec_new_set (double val, idx_t N)
{
  vec  v = vec_new (N);
  vec_set (v, val);
  return (v);
}

ivec ivec_new_set (int val, idx_t N)
{
  ivec v = ivec_new (N);
  ivec_set (v, val);
  return (v);
}

bvec bvec_new_set (byte val, idx_t N)
{
  bvec v = bvec_new (N);
  bvec_set (v, val);
  return (v);
}

cvec cvec_new_set (cplx val, idx_t N)
{
  cvec v = cvec_new (N);
  cvec_set (v, val);
  return (v);
}



/*------------------------------------------------------------------------------*/
vec vec_new_1N (idx_t N)
{
  idx_t i;
  vec  v;
  v = vec_new (N);
  for (i = 0; i < N; i++)
    v[i] = (double) (i + 1);
  return v;
}


ivec ivec_new_1N (idx_t N)
{
  idx_t i;
  ivec v;
  v = ivec_new (N);
  for (i = 0; i < N; i++)
    v[i] = i + 1;
  return v;
}


bvec bvec_new_1N (idx_t N)
{
  idx_t i;
  bvec v;
  v = bvec_new (N);
  for (i = 0; i < N; i++)
    v[i] = i + 1;
  return v;
}


cvec cvec_new_1N (idx_t N)
{
  idx_t i;
  cvec v;
  v = cvec_new (N);
  for (i = 0; i < N; i++) {
    creal (v[i]) = i + 1;
    cimag (v[i]) = 0;
  }
  return v;
}

vec vec_new_range (idx_t N)
{
  idx_t i;
  vec  v;
  v = vec_new (N);
  for (i = 0; i < N; i++)
    v[i] = (double) (i);
  return v;
}


ivec ivec_new_range (idx_t N)
{
  idx_t i;
  ivec v;
  v = ivec_new (N);
  for (i = 0; i < N; i++)
    v[i] = i;
  return v;
}


bvec bvec_new_range (idx_t N)
{
  idx_t i;
  bvec v;
  v = bvec_new (N);
  for (i = 0; i < N; i++)
    v[i] = i;
  return v;
}

cvec cvec_new_range (idx_t N)
{
  idx_t i;
  cvec v;
  v = cvec_new (N);
  for (i = 0; i < N; i++) {
    creal (v[i]) = i;
    cimag (v[i]) = 0;
  }
  return v;
}


/*------------------------------------------------------------------------------*/
vec vec_new_arithm (double start, double incr, idx_t N)
{
  idx_t i;
  double value = start;
  vec  v;
  v = vec_new (N);
  for (i = 0; i < N; i++, value += incr)
    v[i] = value;
  return v;
}


ivec ivec_new_arithm (int start, int incr, idx_t N)
{
  idx_t i;
  int  value = start;
  ivec v;
  v = ivec_new (N);
  for (i = 0; i < N; i++, value += incr)
    v[i] = value;
  return v;
}


bvec bvec_new_arithm (byte start, byte incr, idx_t N)
{
  idx_t i;
  int  value = start;
  bvec v;
  v = bvec_new (N);
  for (i = 0; i < N; i++, value += incr)
    v[i] = value;
  return v;
}

cvec cvec_new_arithm (cplx start, cplx incr, idx_t N)
{
  idx_t i;
  cplx value = start;
  cvec v;
  v = cvec_new (N);
  for (i = 0; i < N; i++) {
    v[i] = value;
    value = cadd (value, incr);
  }
  return v;
}


/*------------------------------------------------------------------------------*/
vec vec_new_geom (double start, double r, idx_t N)
{
  idx_t i;
  double value = start;
  vec  v;
  it_assert (N >= 0, "Sequence should not be of negative length");
  v = vec_new (N);
  for (i = 0; i < N; i++, value *= r)
    v[i] = value;
  return v;
}


ivec ivec_new_geom (int start, int r, idx_t N)
{
  idx_t i;
  int  value = start;
  ivec v;
  it_assert (N >= 0, "Sequence should not be of negative length");
  v = ivec_new (N);
  for (i = 0; i < N; i++, value *= r)
    v[i] = value;
  return v;
}


bvec bvec_new_geom (byte start, byte r, idx_t N)
{
  idx_t i;
  byte value = start;
  bvec v;
  it_assert (N >= 0, "Sequence should not be of negative length");
  v = bvec_new (N);
  for (i = 0; i < N; i++, value *= r)
    v[i] = value;
  return v;
}

cvec cvec_new_geom (cplx start, cplx r, idx_t N)
{
  idx_t i;
  cplx value = start;
  cvec v;
  v = cvec_new (N);
  for (i = 0; i < N; i++) {
    v[i] = value;
    value = cmul (value, r);
  }
  return v;
}

/*------------------------------------------------------------------------------*/
cvec cvec_new_unit_roots (idx_t N)
{
  idx_t k;
  cvec v;

  v = cvec_new (N);

  /* generate e^{2i\pi k / N} = cos(2 k \pi / N) + i sin(2 k \pi / N) */
  for (k = 0; k < N; k++) {
    creal (v[k]) = cos (2 * k * M_PI / N);
    cimag (v[k]) = sin (2 * k * M_PI / N);
  }

  return (v);
}


/*------------------------------------------------------------------------------*/
vec vec_conv (vec v1, vec v2)
{
  int  i1, i2;
  vec  v = vec_new_zeros (vec_length (v1) + vec_length (v2) - 1);

  for (i1 = 0; i1 < vec_length (v1); i1++)
    for (i2 = 0; i2 < vec_length (v2); i2++)
      v[i1 + i2] += v1[i1] * v2[i2];

  return v;
}


/*------------------------------------------------------------------------------*/
ivec ivec_conv (ivec v1, ivec v2)
{
  int  i1, i2;
  ivec v = ivec_new_zeros (ivec_length (v1) + ivec_length (v2) - 1);

  for (i1 = 0; i1 < ivec_length (v1); i1++)
    for (i2 = 0; i2 < ivec_length (v2); i2++)
      v[i1 + i2] += v1[i1] * v2[i2];

  return v;
}


/*---------------------------------------------------------------------------*/
void vec_rand (vec v)
{
  int  i;
  for (i = 0; i < vec_length (v); i++)
    v[i] = it_rand ();
}


void vec_randn (vec v)
{
  int  i;
  for (i = 0; i < vec_length (v); i++)
    v[i] = it_randn ();
}


/*---------------------------------------------------------------------------*/
vec vec_new_rand (idx_t n)
{
  vec  v = vec_new (n);
  vec_rand (v);
  return v;
}


vec vec_new_randn (idx_t n)
{
  vec  v = vec_new (n);
  vec_randn (v);
  return v;
}

/*
Knuth (or Fisher-Yates) shuffle.

See: R. A. Fisher and F. Yates, Example 12,
Statistical Tables, London, 1938.

Generates uniformly random permutations.

Input  : - len  : size of vector
- seed : to initialize the PRNG
Output : a len-long ivec containing a shuffle
between 0 and len-1.

  Allocate memory : YES (ivec[len])
 */
ivec ivec_new_perm (size_t len)
{

  ivec perm;
  idx_t i = 0;
  idx_t p = 0;
  int  n = 0;

  it_assert (len > 1, "Permutation of a scalar is crazy");

  perm = ivec_new_arithm (0, 1, len);

  for (i = 0; i < len; i++) {
    p = i + (unsigned int) (it_rand () * (len - i));
    n = perm[p];
    perm[p] = perm[i];
    perm[i] = n;
  }

  return (perm);
}


/*-----------------------------------------------------------------*/
/* Sparse vectors                                                  */
/*-----------------------------------------------------------------*/

/* convert a vector into a sparse vector */
void vec_to_spvec (vec v, ivec * svi_out, vec * sv_out)
{
  int  n = vec_length (v);
  int  nz = 0, i, j;
  ivec svi;
  vec sv;
  for (i = 0; i < n; i++)
    if (v[i] != 0)
      nz++;
  svi = ivec_new (nz);
  sv = vec_new (nz);

  j = 0;
  for (i = 0; i < n; i++)
    if (v[i] != 0) {
      svi[j] = i;
      sv[j] = v[i];
      j++;
    }
  *svi_out = svi;
  *sv_out = sv;
}


void ivec_to_spivec (ivec v, ivec * svi_out, ivec * sv_out)
{
  int  n = ivec_length (v);
  int  nz = 0, i, j;
  ivec svi, sv;
  for (i = 0; i < n; i++)
    if (v[i] != 0)
      nz++;
  svi = ivec_new (nz);
  sv = ivec_new (nz);
  j = 0;
  for (i = 0; i < n; i++)
    if (v[i] != 0) {
      svi[j] = i;
      sv[j] = v[i];
      j++;
    }
  *svi_out = svi;
  *sv_out = sv;
}


/* convert a sparse vector into a vector */
vec spvec_to_vec (ivec svi, vec sv, int n)
{
  int  nz = ivec_length (svi), i;
  vec v;

  assert (vec_length (sv) == nz);

  if (n < 0) {			/* guess size */
    int  i;
    for (i = 0; i < nz; i++)
      if (svi[i] + 1 > n)
	n = svi[i] + 1;
  }

  v = vec_new_zeros (n);

  for (i = 0; i < nz; i++)
    v[svi[i]] = sv[i];

  return v;
}


ivec spivec_to_ivec (ivec svi, ivec sv, int n)
{
  int  nz = ivec_length (svi), i;
  ivec v;

  assert (ivec_length (sv) == nz);

  if (n < 0) {			/* guess size */
    int  i;
    for (i = 0; i < nz; i++)
      if (svi[i] + 1 > n)
	n = svi[i] + 1;
  }

  v = ivec_new_zeros (n);

  for (i = 0; i < nz; i++)
    v[svi[i]] = sv[i];

  return v;
}


/* Inner product between two sparse vectors */
double spvec_inner_product (ivec svi1, vec sv1, ivec svi2, vec sv2)
{
  double s = 0;
  idx_t i1 = 0, i2 = 0;

  while (1) {

    if (i1 == ivec_length (svi1)) 
      break;

    if (i2 == ivec_length (svi2)) 
      break;

    if (svi1[i1] == svi2[i2]) {
      s += sv1[i1] * sv2[i2];
      i1++;
      i2++;
    }

    else {
      if (svi1[i1] < svi2[i2])
	i1++;

      else
	i2++;
    }
  }

  return s;
}


/* Inner product between two sparse vectors */
int spivec_inner_product (ivec svi1, ivec sv1, ivec svi2, ivec sv2)
{
  int s = 0;
  idx_t i1 = 0, i2 = 0;

  while (1) {

    if (i1 == ivec_length (svi1)) 
      break;

    if (i2 == ivec_length (svi2)) 
      break;

    if (svi1[i1] == svi2[i2]) {
      s += sv1[i1] * sv2[i2];
      i1++;
      i2++;
    }

    else {
      if (svi1[i1] < svi2[i2])
	i1++;

      else
	i2++;
    }
  }

  return s;
}
