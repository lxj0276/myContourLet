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
  Matrices functions
  Copyright (C) 2005 Vivien Chappelier, Herve Jegou
*/

#include <it/mat.h>
#include <it/math.h>
#include <it/io.h>
#include <it/random.h>

#define __mat_c_my_min(x,y)  (x>y?y:x)


/*---------------------------------------------------------------------------*/
/* Basic vector allocation, initialization and deletion                      */
/*---------------------------------------------------------------------------*/

Mat __Mat_new_alloc (size_t elem_size, idx_t h, idx_t w,
		     idx_t hmax, idx_t wmax) {
  idx_t i;
  Mat  m;

  if   (hmax <= 0)		/* The maximum height is at least 1 */
    hmax = 1;
  m = Vec_new_alloc (Vec, h, hmax);
  for  (i = 0; i < hmax; i++)
    m[i] = __Vec_new_alloc (elem_size, w, wmax);
  return (m);
}


Mat __Mat_new (size_t elem_size, idx_t h, idx_t w) {
  return (__Mat_new_alloc (elem_size, h, w, h, w));
}


void __Mat_delete (Mat m) {
  idx_t i;
  for (i = 0; i < Mat_height_max (m); i++)
    Vec_delete (m[i]);
  Vec_delete (m);
}


mat mat_new (idx_t h, idx_t w) {
  return (Mat_new (double, h, w));
}


imat imat_new (idx_t h, idx_t w) {
  return (Mat_new (int, h, w));
}


bmat bmat_new (idx_t h, idx_t w) {
  return (Mat_new (byte, h, w));
}


cmat cmat_new (idx_t h, idx_t w) {
  return (Mat_new (cplx, h, w));
}


mat mat_new_alloc (idx_t h, idx_t w, idx_t hmax, idx_t wmax) {
  return (Mat_new_alloc (double, h, w, hmax, wmax));
}


imat imat_new_alloc (idx_t h, idx_t w, idx_t hmax, idx_t wmax) {
  return (Mat_new_alloc (int, h, w, hmax, wmax));
}


bmat bmat_new_alloc (idx_t h, idx_t w, idx_t hmax, idx_t wmax) {
  return (Mat_new_alloc (byte, h, w, hmax, wmax));
}


cmat cmat_new_alloc (idx_t h, idx_t w, idx_t hmax, idx_t wmax) {
  return (Mat_new_alloc (cplx, h, w, hmax, wmax));
}


void mat_delete (mat m) {
  Mat_delete (m);
}


void imat_delete (imat m) {
  Mat_delete (m);
}


void bmat_delete (bmat m) {
  Mat_delete (m);
}


void cmat_delete (cmat m) {
  Mat_delete (m);
}


/*---------------------------------------------------------------------------*/
/* Width and Height functions                                                */
/*---------------------------------------------------------------------------*/

idx_t mat_height (mat m) {
  return Mat_height (m);
}


idx_t mat_height_max (mat m) {
  return Mat_height_max (m);
}


idx_t mat_width (mat m) {
  return Mat_width (m);
}


idx_t mat_width_max (mat m) {
  return Mat_width_max (m);
}


idx_t imat_height (imat m) {
  return Mat_height (m);
}


idx_t imat_height_max (imat m) {
  return Mat_height_max (m);
}


idx_t imat_width (imat m) {
  return Mat_width (m);
}


idx_t imat_width_max (imat m) {
  return Mat_width_max (m);
}


idx_t bmat_height (bmat m) {
  return Mat_height (m);
}


idx_t bmat_height_max (bmat m) {
  return Mat_height_max (m);
}


idx_t bmat_width (bmat m) {
  return Mat_width (m);
}


idx_t bmat_width_max (bmat m) {
  return Mat_width_max (m);
}


idx_t cmat_height (cmat m) {
  return Mat_height (m);
}


idx_t cmat_height_max (cmat m) {
  return Mat_height_max (m);
}


idx_t cmat_width (cmat m) {
  return Mat_width (m);
}


idx_t cmat_width_max (cmat m) {
  return Mat_width_max (m);
}


mat _mat_set_height (mat m, idx_t h) {
  Mat_set_height (m, h);
  return m;
}


mat _mat_set_height_max (mat m, idx_t hmax) {
  Mat_set_height_max (m, hmax);
  return m;
}


imat _imat_set_height (imat m, idx_t h) {
  Mat_set_height (m, h);
  return m;
}


imat _imat_set_height_max (imat m, idx_t hmax) {
  Mat_set_height_max (m, hmax);
  return m;
}


bmat _bmat_set_height (bmat m, idx_t h) {
  Mat_set_height (m, h);
  return m;
}


bmat _bmat_set_height_max (bmat m, idx_t hmax) {
  Mat_set_height_max (m, hmax);
  return (m);
}


cmat _cmat_set_height (cmat m, idx_t h) {
  Mat_set_height (m, h);
  return m;
}


cmat _cmat_set_height_max (cmat m, idx_t hmax) {
  Mat_set_height_max (m, hmax);
  return (m);
}


/*---------------------------------------------------------------------------*/
/* Basic vector initialization/modification                                  */
/*---------------------------------------------------------------------------*/

mat _mat_init (mat m, double *buf, idx_t w, idx_t h) {
  Mat_init (m, buf, w, h);
  return (m);
}


imat _imat_init (imat m, int *buf, idx_t w, idx_t h) {
  Mat_init (m, buf, w, h);
  return (m);
}


bmat _bmat_init (bmat m, byte * buf, idx_t w, idx_t h) {
  Mat_init (m, buf, w, h);
  return (m);
}


cmat _cmat_init (cmat m, cplx * buf, idx_t w, idx_t h) {
  Mat_init (m, buf, w, h);
  return (m);
}

mat mat_set (mat m, double val) {
  Mat_set (m, val);
  return m;
}


imat imat_set (imat m, int val) {
  Mat_set (m, val);
  return m;
}


bmat bmat_set (bmat m, byte val) {
  Mat_set (m, val);
  return m;
}


cmat cmat_set (cmat m, cplx val) {
  Mat_set (m, val);
  return m;
}


mat mat_set_between (mat m, idx_t r1, idx_t c1, idx_t r2,
		     idx_t c2, double val) {
  Mat_set_between (m, r1, c1, r2, c2, val);
  return m;
}


imat imat_set_between (imat m, idx_t r1, idx_t c1, idx_t r2,
		       idx_t c2, int val) {
  Mat_set_between (m, r1, c1, r2, c2, val);
  return m;
}


bmat bmat_set_between (bmat m, idx_t r1, idx_t c1, idx_t r2,
		       idx_t c2, byte val) {
  Mat_set_between (m, r1, c1, r2, c2, val);
  return m;
}


cmat cmat_set_between (cmat m, idx_t r1, idx_t c1, idx_t r2,
		       idx_t c2, cplx val) {
  Mat_set_between (m, r1, c1, r2, c2, val);
  return m;
}


void mat_void (mat m)
{
  mat_set_height (m, 0);
}


void imat_void (imat m)
{
  imat_set_height (m, 0);
}


void bmat_void (bmat m)
{
  bmat_set_height (m, 0);
}


void cmat_void (cmat m)
{
  cmat_set_height (m, 0);
}


/*----------------------------------------------------------------------------*/
void mat_zeros (mat m)
{
  mat_set (m, 0);
}


void imat_zeros (imat m)
{
  imat_set (m, 0);
}


void bmat_zeros (bmat m)
{
  bmat_set (m, 0);
}


void cmat_zeros (cmat m)
{
  cmat_set (m, cplx_0);
}


mat mat_new_zeros (idx_t h, idx_t w) 
{
  return (mat_new_set (0., h, w));
}


imat imat_new_zeros (idx_t h, idx_t w) {
  return (imat_new_set (0, h, w));
}


bmat bmat_new_zeros (idx_t h, idx_t w) {
  return (bmat_new_set (0, h, w));
}


cmat cmat_new_zeros (idx_t h, idx_t w) {
  return (cmat_new_set (cplx_0, h, w));
}



/*----------------------------------------------------------------------------*/
void mat_eye (mat m)
{
  idx_t i, mi;
  mat_set (m, 0);
  mi = (mat_height (m) > mat_width (m) ? mat_width (m) : mat_height (m));
  for (i = 0; i < mi; i++)
    m[i][i] = 1;
}


void imat_eye (imat m)
{
  idx_t i, mi;
  imat_set (m, 0);
  mi = (imat_height (m) > imat_width (m) ? imat_width (m) : imat_height (m));
  for (i = 0; i < mi; i++)
    m[i][i] = 1;
}


void bmat_eye (bmat m)
{
  idx_t i, mi;
  bmat_set (m, 0);
  mi = (bmat_height (m) > bmat_width (m) ? bmat_width (m) : bmat_height (m));
  for (i = 0; i < mi; i++)
    m[i][i] = 1;
}


void cmat_eye (cmat m)
{
  idx_t i, mi;
  cmat_set (m, cplx_0);
  mi = (cmat_height (m) > cmat_width (m) ? cmat_width (m) : cmat_height (m));
  for (i = 0; i < mi; i++)
    m[i][i] = cplx_1;
}


/*----------------------------------------------------------------------------*/
void mat_ones (mat m)
{
  mat_set (m, 1);
}


void imat_ones (imat m)
{
  imat_set (m, 1);
}


void bmat_ones (bmat m)
{
  bmat_set (m, 1);
}


void cmat_ones (cmat m)
{
  cmat_set (m, cplx_1);
}


mat mat_new_ones (idx_t h, idx_t w)
{
  return (mat_new_set (1., h, w));
}


imat imat_new_ones (idx_t h, idx_t w)
{
  return (imat_new_set (1, h, w));
}


bmat bmat_new_ones (idx_t h, idx_t w) 
{
  return (bmat_new_set (1, h, w));
}


cmat cmat_new_ones (idx_t h, idx_t w) 
{
  return (cmat_new_set (cplx_1, h, w));
}


/*----------------------------------------------------------------------------*/
mat mat_new_void ()
{
  mat  m = mat_new (0, 0);
  return m;
}


imat imat_new_void ()
{
  imat m = imat_new (0, 0);
  return m;
}


bmat bmat_new_void ()
{
  bmat m = bmat_new (0, 0);
  return m;
}


cmat cmat_new_void ()
{
  cmat m = cmat_new (0, 0);
  return m;
}


/*---------------------------------------------------------------------------*/
mat mat_new_set (double val, idx_t h, idx_t w)
{
  mat  m = mat_new (h, w);
  mat_set (m, val);
  return m;
}


imat imat_new_set (int val, idx_t h, idx_t w)
{
  imat m = imat_new (h, w);
  imat_set (m, val);
  return m;
}


bmat bmat_new_set (byte val, idx_t h, idx_t w)
{
  bmat m = bmat_new (h, w);
  bmat_set (m, val);
  return m;
}


cmat cmat_new_set (cplx val, idx_t h, idx_t w)
{
  cmat m = cmat_new (h, w);
  cmat_set (m, val);
  return m;
}


/*---------------------------------------------------------------------------*/
mat mat_new_eye (idx_t n)
{
  mat  m = mat_new (n, n);
  mat_eye (m);
  return m;
}


imat imat_new_eye (idx_t n)
{
  imat m = imat_new (n, n);
  imat_eye (m);
  return m;
}


bmat bmat_new_eye (idx_t n)
{
  bmat m = bmat_new (n, n);
  bmat_eye (m);
  return m;
}


cmat cmat_new_eye (idx_t n)
{
  cmat m = cmat_new (n, n);
  cmat_eye (m);
  return m;
}

/*---------------------------------------------------------------------------*/
mat mat_new_diag (vec v)
{
  mat  m = mat_new (vec_length (v), vec_length (v));
  mat_diag (m, v);
  return m;
}


imat imat_new_diag (ivec v)
{
  imat m = imat_new (ivec_length (v), ivec_length (v));
  imat_diag (m, v);
  return m;
}


bmat bmat_new_diag (bvec v)
{
  bmat m = bmat_new (bvec_length (v), bvec_length (v));
  bmat_diag (m, v);
  return m;
}


cmat cmat_new_diag (cvec v)
{
  cmat m = cmat_new (cvec_length (v), cvec_length (v));
  cmat_diag (m, v);
  return m;
}


/*---------------------------------------------------------------------------*/
/*               Copy and Conversions Functions                              */
/*---------------------------------------------------------------------------*/

Mat __Mat_copy (Mat m1, Mat m2) {
  int  i;
  assert (m1);
  assert (m2);
  assert (Mat_element_size (m1) == Mat_element_size (m2));
  assert (Mat_height (m1) >= Mat_height (m2));

  for (i = 0; i < Mat_height (m2); i++)
    Vec_copy (m1[i], m2[i]);
  return (m1);
}


void mat_copy (mat dest, mat orig)
{
  idx_t i, j;
  assert (dest);
  assert (orig);
  assert (mat_width_max (dest) >= mat_width (orig));
  assert (mat_height_max (dest) >= mat_height (orig));
  for (i = 0; i < mat_height (orig); i++)
    for (j = 0; j < mat_width (orig); j++)
      dest[i][j] = orig[i][j];
}


void imat_copy (imat dest, imat orig)
{
  idx_t i, j;
  assert (dest);
  assert (orig);
  assert (imat_width_max (dest) >= imat_width (orig));
  assert (imat_height_max (dest) >= imat_height (orig));
  for (i = 0; i < imat_height (orig); i++)
    for (j = 0; j < imat_width (orig); j++)
      dest[i][j] = orig[i][j];
}


void bmat_copy (bmat dest, bmat orig)
{
  idx_t i, j;
  assert (dest);
  assert (orig);
  assert (bmat_width_max (dest) >= bmat_width (orig));
  assert (bmat_height_max (dest) >= bmat_height (orig));
  for (i = 0; i < bmat_height (orig); i++)
    for (j = 0; j < bmat_width (orig); j++)
      dest[i][j] = orig[i][j];
}

void cmat_copy (cmat dest, cmat orig)
{
  idx_t i, j;
  assert (dest);
  assert (orig);
  assert (cmat_width_max (dest) >= cmat_width (orig));
  assert (cmat_height_max (dest) >= cmat_height (orig));
  for (i = 0; i < cmat_height (orig); i++)
    for (j = 0; j < cmat_width (orig); j++)
      dest[i][j] = orig[i][j];
}


/*------------------------------------------------------------------*/
Mat Mat_clone (Mat m)
{
  Mat  cl;
  assert (m);
  cl = __Mat_new (Mat_element_size (m), Mat_height (m), Mat_width (m));
  if (cl == NULL)
    return NULL;
  Mat_copy (cl, m);
  return cl;
}

mat mat_clone (mat m)
{
  mat  cl;
  assert (m);
  cl = mat_new (mat_height (m), mat_width (m));
  if (cl == NULL)
    return NULL;
  mat_copy (cl, m);
  return cl;
}


imat imat_clone (imat m)
{
  imat cl;
  assert (m);
  cl = imat_new (imat_height (m), imat_width (m));
  if (cl == NULL)
    return NULL;
  imat_copy (cl, m);
  return cl;

}


bmat bmat_clone (bmat m)
{
  bmat cl;
  assert (m);
  cl = bmat_new (bmat_height (m), bmat_width (m));
  if (cl == NULL)
    return NULL;
  bmat_copy (cl, m);
  return cl;
}


cmat cmat_clone (cmat m)
{
  cmat cl;
  assert (m);
  cl = cmat_new (cmat_height (m), cmat_width (m));
  if (cl == NULL)
    return NULL;
  cmat_copy (cl, m);
  return cl;
}


vec mat_to_vec (mat m)
{
  idx_t i, j;
  vec  v = vec_new (mat_height (m) * mat_width (m));
  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      v[i * mat_width (m) + j] = m[i][j];
  return v;
}

ivec imat_to_ivec (imat m)
{
  idx_t i, j;
  ivec v = ivec_new (imat_height (m) * imat_width (m));
  for (i = 0; i < imat_height (m); i++)
    for (j = 0; j < imat_width (m); j++)
      v[i * imat_width (m) + j] = m[i][j];
  return v;
}

bvec bmat_to_bvec (bmat m)
{
  idx_t i, j;
  bvec v = bvec_new (bmat_height (m) * bmat_width (m));
  for (i = 0; i < bmat_height (m); i++)
    for (j = 0; j < bmat_width (m); j++)
      v[i * bmat_width (m) + j] = m[i][j];
  return v;
}

cvec cmat_to_cvec (cmat m)
{
  idx_t i, j;
  cvec v = cvec_new (cmat_height (m) * cmat_width (m));
  for (i = 0; i < cmat_height (m); i++)
    for (j = 0; j < cmat_width (m); j++)
      v[i * cmat_width (m) + j] = m[i][j];
  return v;
}

/*----------------------------------------------------------------------------*/
mat vec_to_mat (vec v, idx_t width)
{
  mat  m;
  idx_t i, j, k = 0;
  it_assert (vec_length (v) % width == 0,
	     "Vector dimension does match the required matrix width\n");

  m = mat_new (vec_length (v) / width, width);
  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      m[i][j] = v[k++];
  return m;
}

imat ivec_to_imat (ivec v, idx_t width)
{
  imat m;
  idx_t i, j, k = 0;
  it_assert (ivec_length (v) % width == 0,
	     "Vector dimension does match the required matrix width\n");

  m = imat_new (ivec_length (v) / width, width);
  for (i = 0; i < imat_height (m); i++)
    for (j = 0; j < imat_width (m); j++)
      m[i][j] = v[k++];
  return m;
}

bmat bvec_to_bmat (bvec v, idx_t width)
{
  bmat m;
  idx_t i, j, k = 0;
  it_assert (bvec_length (v) % width == 0,
	     "Vector dimension does match the required matrix width\n");

  m = bmat_new (bvec_length (v) / width, width);
  for (i = 0; i < bmat_height (m); i++)
    for (j = 0; j < bmat_width (m); j++)
      m[i][j] = v[k++];
  return m;
}

cmat cvec_to_cmat (cvec v, idx_t width)
{
  cmat m;
  idx_t i, j, k = 0;
  it_assert (cvec_length (v) % width == 0,
	     "Vector dimension does match the required matrix width\n");

  m = cmat_new (cvec_length (v) / width, width);
  for (i = 0; i < cmat_height (m); i++)
    for (j = 0; j < cmat_width (m); j++)
      m[i][j] = v[k++];
  return m;
}

/*----------------------------------------------------------------------------*/
mat imat_to_mat (imat m)
{
  idx_t i, j;
  mat  cl;

  assert (m);
  cl = mat_new (imat_height (m), imat_width (m));

  for (i = 0; i < imat_height (m); i++)
    for (j = 0; j < imat_width (m); j++)
      cl[i][j] = m[i][j];

  return cl;
}



mat bmat_to_mat (bmat m)
{
  idx_t i, j;
  mat  cl;

  assert (m);
  cl = mat_new (bmat_height (m), bmat_width (m));

  for (i = 0; i < bmat_height (m); i++)
    for (j = 0; j < bmat_width (m); j++)
      cl[i][j] = m[i][j];

  return cl;
}


imat bmat_to_imat (bmat m)
{
  idx_t i, j;
  imat cl;

  assert (m);
  cl = imat_new (bmat_height (m), bmat_width (m));

  for (i = 0; i < bmat_height (m); i++)
    for (j = 0; j < bmat_width (m); j++)
      cl[i][j] = m[i][j];

  return cl;
}


imat mat_to_imat (mat m)
{
  idx_t i, j;
  imat cl;

  assert (m);
  cl = imat_new (mat_height (m), mat_width (m));

  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      cl[i][j] = (int) m[i][j];

  return cl;
}


/*---------------------------------------------------------------------------*/
/* Comparison functions                                                      */
/*---------------------------------------------------------------------------*/

int Mat_eq (Mat m1, Mat m2)
{
  idx_t i;
  assert (m1);
  assert (m2);

  if (Mat_height (m1) != Mat_height (m2))
    return 0;

  for (i = 0; i < Mat_height (m1); i++)
    if (!Vec_eq (m1[i], m2[i]))
      return 0;

  return 1;
}


int mat_eq (mat m1, mat m2)
{
  idx_t i;
  assert (m1);
  assert (m2);

  if (mat_height (m1) != mat_height (m2))
    return 0;

  for (i = 0; i < mat_height (m1); i++)
    if (!vec_eq (m1[i], m2[i]))
      return 0;
  return 1;
}

int imat_eq (imat m1, imat m2)
{
  idx_t i;
  assert (m1);
  assert (m2);

  if (imat_height (m1) != imat_height (m2))
    return 0;

  for (i = 0; i < imat_height (m1); i++)
    if (!ivec_eq (m1[i], m2[i]))
      return 0;
  return 1;
}

int bmat_eq (bmat m1, bmat m2)
{
  idx_t i;
  assert (m1);
  assert (m2);

  if (bmat_height (m1) != bmat_height (m2))
    return 0;

  for (i = 0; i < bmat_height (m1); i++)
    if (!bvec_eq (m1[i], m2[i]))
      return 0;
  return 1;
}

int cmat_eq (cmat m1, cmat m2)
{
  idx_t i;
  assert (m1);
  assert (m2);

  if (cmat_height (m1) != cmat_height (m2))
    return 0;

  for (i = 0; i < cmat_height (m1); i++)
    if (!cvec_eq (m1[i], m2[i]))
      return 0;
  return 1;
}


int mat_is_void (mat m)
{
  return mat_width (m) == 0 || mat_height (m) == 0;
}


int imat_is_void (imat m)
{
  return imat_width (m) == 0 || imat_height (m) == 0;
}


int bmat_is_void (bmat m)
{
  return bmat_width (m) == 0 || bmat_height (m) == 0;
}


int cmat_is_void (cmat m)
{
  return cmat_width (m) == 0 || cmat_height (m) == 0;
}


/*---------------------------------------------------------------------------*/
/** @name Row/column manipulation                                            */
/*---------------------------------------------------------------------------*/

mat __mat_del_row (mat m, idx_t pos) {
  Mat_del_row (m, pos);
  return m;
}
imat __imat_del_row (imat m, idx_t pos) {
  Mat_del_row (m, pos);
  return m;
}
bmat __bmat_del_row (bmat m, idx_t pos) {
  Mat_del_row (m, pos);
  return m;
}
cmat __cmat_del_row (cmat m, idx_t pos) {
  Mat_del_row (m, pos);
  return m;
}

mat __mat_ins_row (mat m, idx_t pos, vec v) {
  Mat_ins_row (m, pos, v);
  return m;
}
imat __imat_ins_row (imat m, idx_t pos, ivec v) {
  Mat_ins_row (m, pos, v);
  return m;
}
bmat __bmat_ins_row (bmat m, idx_t pos, bvec v) {
  Mat_ins_row (m, pos, v);
  return m;
}
cmat __cmat_ins_row (cmat m, idx_t pos, cvec v) {
  Mat_ins_row (m, pos, v);
  return m;
}

/* Add a vector at the end row of the matrix (Use matrix as a vector stack)  */
mat __mat_push_row (mat m, vec v) {
  mat_ins_row (m, mat_height (m), v);
  return m;
}
imat __imat_push_row (imat m, ivec v) {
  imat_ins_row (m, imat_height (m), v);
  return m;
}
bmat __bmat_push_row (bmat m, bvec v) {
  bmat_ins_row (m, bmat_height (m), v);
  return m;
}
cmat __cmat_push_row (cmat m, cvec v) {
  cmat_ins_row (m, cmat_height (m), v);
  return m;
}

/* Retrieve the last element at the end of the vector (Use vector as a stack)*/
mat __mat_pop_row (mat m) {
  Mat_pop_row (m);
  return m;
}
imat __imat_pop_row (imat m) {
  Mat_pop_row (m);
  return m;
}
bmat __bmat_pop_row (bmat m) {
  Mat_pop_row (m);
  return m;
}
cmat __cmat_pop_row (cmat m) {
  Mat_pop_row (m);
  return m;
}


/*----------------------------------------------------------------------------*/

void mat_swap_rows (mat m, idx_t i, idx_t j)
{

  vec r;

  it_assert (m, "Please use an existing matrix");
  it_assert (0 <= i && i <= mat_height (m) &&
	     0 <= j && j <= mat_height (m), "Please specify existing rows");

  r = m[i];

  m[i] = m[j];
  m[j] = r;

  return;

}


void mat_swap_cols (mat m, idx_t i, idx_t j)
{

  idx_t r;
  double t;

  it_assert (m, "Please use an existing matrix");
  it_assert (0 <= i && i <= mat_width (m) &&
	     0 <= j && j <= mat_width (m), "Please specify existing columns");

  for (r = 0; r < mat_height (m); r++) {
    t = m[r][i];
    m[r][i] = m[r][j];
    m[r][j] = t;
  }

  return;

}


void imat_swap_rows (imat m, idx_t i, idx_t j)
{
  ivec r;

  it_assert (m, "Please use an existing matrix");
  it_assert (0 <= i && i <= imat_height (m) &&
	     0 <= j && j <= imat_height (m), "Please specify existing rows");

  r = m[i];

  m[i] = m[j];
  m[j] = r;
}


void imat_swap_cols (imat m, idx_t i, idx_t j)
{
  idx_t r;
  int  t;

  it_assert (m, "Please use an existing matrix");
  it_assert (0 <= i && i <= imat_width (m) &&
	     0 <= j && j <= imat_width (m), "Please specify existing columns");

  for (r = 0; r < imat_height (m); r++) {
    t = m[r][i];
    m[r][i] = m[r][j];
    m[r][j] = t;
  }
}


void bmat_swap_rows (bmat m, idx_t i, idx_t j)
{
  bvec r;

  it_assert (m, "Please use an existing matrix");
  it_assert (0 <= i && i <= bmat_height (m) &&
	     0 <= j && j <= bmat_height (m), "Please specify existing rows");

  r = m[i];

  m[i] = m[j];
  m[j] = r;
}


void bmat_swap_cols (bmat m, idx_t i, idx_t j)
{
  idx_t r;
  byte t;

  it_assert (m, "Please use an existing matrix");
  it_assert (0 <= i && i <= bmat_width (m) &&
	     0 <= j && j <= bmat_width (m), "Please specify existing columns");

  for (r = 0; r < bmat_height (m); r++) {
    t = m[r][i];
    m[r][i] = m[r][j];
    m[r][j] = t;
  }
}


void cmat_swap_rows (cmat m, idx_t i, idx_t j)
{
  cvec r;

  it_assert (m, "Please use an existing matrix");
  it_assert (0 <= i && i <= cmat_height (m) &&
	     0 <= j && j <= cmat_height (m), "Please specify existing rows");

  r = m[i];

  m[i] = m[j];
  m[j] = r;
}


void cmat_swap_cols (cmat m, idx_t i, idx_t j)
{
  idx_t r;
  cplx t;

  it_assert (m, "Please use an existing matrix");
  it_assert (0 <= i && i <= cmat_width (m) &&
	     0 <= j && j <= cmat_width (m), "Please specify existing columns");

  for (r = 0; r < cmat_height (m); r++) {
    t = m[r][i];
    m[r][i] = m[r][j];
    m[r][j] = t;
  }
}


mat _mat_transpose (mat m)
{
  idx_t i, j;
  double val;

  it_assert (mat_height (m) == mat_width (m),
	     "This function work with symmetric matrices only");

  for (i = 0; i < mat_width (m); i++)
    for (j = 0; j < mat_width (m); j++) {
      val = m[i][j];
      m[i][j] = m[j][i];
      m[j][i] = val;
    }

  return m;
}


imat _imat_transpose (imat m)
{
  idx_t i, j;
  int  val;

  it_assert (imat_height (m) == imat_width (m),
	     "This function work with symmetric matrices only");

  for (i = 0; i < imat_width (m); i++)
    for (j = 0; j < imat_width (m); j++) {
      val = m[i][j];
      m[i][j] = m[j][i];
      m[j][i] = val;
    }

  return m;
}


bmat _bmat_transpose (bmat m)
{
  idx_t i, j;
  byte val;

  it_assert (bmat_height (m) == bmat_width (m),
	     "This function work with symmetric matrices only");

  for (i = 0; i < bmat_width (m); i++)
    for (j = 0; j < bmat_width (m); j++) {
      val = m[i][j];
      m[i][j] = m[j][i];
      m[j][i] = val;
    }

  return m;
}


/*----------------------------------------------------------------------------*/
mat mat_new_transpose (mat m)
{
  idx_t i, j;
  mat  t = mat_new (mat_width (m), mat_height (m));

  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      t[j][i] = m[i][j];
  return t;
}


imat imat_new_transpose (imat m)
{
  idx_t i, j;
  imat t = imat_new (imat_width (m), imat_height (m));

  for (i = 0; i < imat_height (m); i++)
    for (j = 0; j < imat_width (m); j++)
      t[j][i] = m[i][j];
  return t;
}


bmat bmat_new_transpose (bmat m)
{
  idx_t i, j;
  bmat t = bmat_new (bmat_width (m), bmat_height (m));

  for (i = 0; i < bmat_height (m); i++)
    for (j = 0; j < bmat_width (m); j++)
      t[j][i] = m[i][j];
  return t;
}


/*------------------------------------------------------------------------------*/
/* Mean, median, min, max, sum and elementary statistic indicators              */
/*------------------------------------------------------------------------------*/


double mat_sum (mat m)
{
  idx_t i, j;
  double s = 0;
  assert (m);
  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      s += m[i][j];
  return s;
}


int imat_sum (imat m)
{
  idx_t i, j;
  int  s = 0;
  assert (m);
  for (i = 0; i < imat_height (m); i++)
    for (j = 0; j < imat_width (m); j++)
      s += m[i][j];
  return s;
}


cplx cmat_sum (cmat m)
{
  idx_t i, j;
  cplx s = cplx_0;
  assert (m);
  for (i = 0; i < cmat_height (m); i++)
    for (j = 0; j < cmat_width (m); j++)
      s = cadd (s, m[i][j]);
  return s;
}


/*----------------------------------------------------------------------------*/
double mat_row_sum (mat m, idx_t c)
{
  MAT_END_ROW_PARAM (m, c);
  return vec_sum (m[c]);
}


double mat_col_sum (mat m, idx_t c)
{
  idx_t i;
  double s = 0;
  MAT_END_COL_PARAM (m, c);
  assert (m);
  for (i = 0; i < mat_height (m); i++)
    s += m[i][c];
  return s;
}


vec mat_rows_sum (mat m)
{
  vec  s = vec_new (mat_height (m));
  idx_t i;
  for (i = 0; i < vec_length (s); i++)
    s[i] = mat_row_sum (m, i);
  return s;
}


vec mat_cols_sum (mat m)
{
  vec  s = vec_new (mat_width (m));
  idx_t i;
  for (i = 0; i < vec_length (s); i++)
    s[i] = mat_col_sum (m, i);
  return s;
}


int imat_row_sum (imat m, idx_t c)
{
  MAT_END_ROW_PARAM (m, c);
  return ivec_sum (m[c]);
}


int imat_col_sum (imat m, idx_t c)
{
  idx_t i;
  int  s = 0;
  MAT_END_COL_PARAM (m, c);
  assert (m);
  for (i = 0; i < imat_height (m); i++)
    s += m[i][c];
  return s;
}


ivec imat_rows_sum (imat m)
{
  ivec s = ivec_new (imat_height (m));
  idx_t i;
  for (i = 0; i < ivec_length (s); i++)
    s[i] = imat_row_sum (m, i);
  return s;
}


ivec imat_cols_sum (imat m)
{
  ivec s = ivec_new (imat_width (m));
  idx_t i;
  for (i = 0; i < ivec_length (s); i++)
    s[i] = imat_col_sum (m, i);
  return s;
}


cplx cmat_row_sum (cmat m, idx_t c)
{
  MAT_END_ROW_PARAM (m, c);
  return cvec_sum (m[c]);
}


cplx cmat_col_sum (cmat m, idx_t c)
{
  idx_t i;
  cplx s = cplx_0;
  MAT_END_COL_PARAM (m, c);
  assert (m);
  for (i = 0; i < cmat_height (m); i++)
    s = cadd (s, m[i][c]);
  return s;
}


cvec cmat_rows_sum (cmat m)
{
  cvec s = cvec_new (cmat_height (m));
  idx_t i;
  for (i = 0; i < cvec_length (s); i++)
    s[i] = cmat_row_sum (m, i);
  return s;
}


cvec cmat_cols_sum (cmat m)
{
  cvec s = cvec_new (cmat_width (m));
  idx_t i;
  for (i = 0; i < cvec_length (s); i++)
    s[i] = cmat_col_sum (m, i);
  return s;
}


double mat_diag_sum (mat m)
{
  double res = 0;
  int i, n = __mat_c_my_min (mat_height (m), mat_width (m));

  for (i = 0 ; i < n ; i++)
    res += m[i][i];

  return res;
}


int imat_diag_sum (imat m)
{
  int res = 0;
  int i, n = __mat_c_my_min (imat_height (m), imat_width (m));

  for (i = 0 ; i < n ; i++)
    res += m[i][i];

  return res;
}


int bmat_diag_sum (bmat m)
{
  int res = 0;
  int i, n = __mat_c_my_min (bmat_height (m), bmat_width (m));

  for (i = 0 ; i < n ; i++)
    res += m[i][i];

  return res;
}


cplx cmat_diag_sum (cmat m)
{
  cplx res = cplx_0;
  int i, n = __mat_c_my_min (cmat_height (m), cmat_width (m));

  for (i = 0 ; i < n ; i++)
    cadd (res, m[i][i]);

  return res;
}




/*----------------------------------------------------------------------------*/
double mat_mean (mat m)
{
  assert (m);
  return mat_sum (m) / (mat_width (m) * mat_height (m));
}


double imat_mean (imat m)
{
  assert (m);
  return imat_sum (m) / (double) (imat_width (m) * imat_height (m));
}


cplx cmat_mean (cmat m)
{
  cplx s;
  assert (m);
  s = cmat_sum (m);
  creal (s) /= (cmat_width (m) * cmat_height (m));
  cimag (s) /= (cmat_width (m) * cmat_height (m));
  return s;
}


/*----------------------------------------------------------------------------*/
double mat_max_index_submatrix (mat m, idx_t rmin, idx_t rmax, idx_t cmin,
				idx_t cmax, idx_t * r, idx_t * c)
{

  double max = -HUGE_VAL;
  idx_t i = 0;
  idx_t j = 0;

  it_assert (m, "Please use an existing matrix");
  it_assert (rmin >= 0 && rmin <= rmax
	     && rmax <= mat_height (m), "Wrong row indexes");
  it_assert (cmin >= 0 && cmin <= cmax
	     && cmax <= mat_width (m), "Wrong column indexes");

  for (j = cmin; j < cmax; j++)
    for (i = rmin; i < rmax; i++)
      if (max < m[i][j]) {
	*r = i;
	*c = j;
	max = m[i][j];
      }

  return max;

}

double mat_min_index_submatrix (mat m, idx_t rmin, idx_t rmax, idx_t cmin,
				idx_t cmax, idx_t * r, idx_t * c)
{

  double min = HUGE_VAL;
  idx_t i = 0;
  idx_t j = 0;

  it_assert (m, "Please use an existing matrix");
  it_assert (rmin >= 0 && rmin <= rmax
	     && rmax <= mat_height (m), "Wrong row indexes");
  it_assert (cmin >= 0 && cmin <= cmax
	     && cmax <= mat_width (m), "Wrong column indexes");

  for (j = cmin; j < cmax; j++)
    for (i = rmin; i < rmax; i++)
      if (min > m[i][j]) {
	*r = i;
	*c = j;
	min = m[i][j];
      }

  return min;

}

int imat_max_index_submatrix (imat m, idx_t rmin, idx_t rmax, idx_t cmin,
			      idx_t cmax, idx_t * r, idx_t * c)
{

  int  max = INT_MIN;
  idx_t i = 0;
  idx_t j = 0;

  it_assert (m, "Please use an existing matrix");
  it_assert (rmin >= 0 && rmin <= rmax
	     && rmax <= imat_height (m), "Wrong row indexes");
  it_assert (cmin >= 0 && cmin <= cmax
	     && cmax <= imat_width (m), "Wrong column indexes");

  for (j = cmin; j < cmax; j++)
    for (i = rmin; i < rmax; i++)
      if (max < m[i][j]) {
	*r = i;
	*c = j;
	max = m[i][j];
      }

  return max;

}

int imat_min_index_submatrix (imat m, idx_t rmin, idx_t rmax, idx_t cmin,
			      idx_t cmax, idx_t * r, idx_t * c)
{

  int  min = INT_MAX;
  idx_t i = 0;
  idx_t j = 0;

  it_assert (m, "Please use an existing matrix");
  it_assert (rmin >= 0 && rmin <= rmax
	     && rmax <= imat_height (m), "Wrong row indexes");
  it_assert (cmin >= 0 && cmin <= cmax
	     && cmax <= imat_width (m), "Wrong column indexes");

  for (j = cmin; j < cmax; j++)
    for (i = rmin; i < rmax; i++)
      if (min > m[i][j]) {
	*r = i;
	*c = j;
	min = m[i][j];
      }

  return min;

}

double mat_max (mat m)
{
  idx_t r, c;

  return mat_max_index_submatrix (m, 0, mat_height (m), 0, mat_width (m),
				  &r, &c);

}


int imat_max (imat m)
{
  idx_t r, c;

  return imat_max_index_submatrix (m, 0, imat_height (m), 0,
				   imat_width (m), &r, &c);

}


double mat_min (mat m)
{
  idx_t r, c;

  return mat_min_index_submatrix (m, 0, mat_height (m), 0, mat_width (m),
				  &r, &c);

}


int imat_min (imat m)
{
  idx_t r, c;

  return imat_min_index_submatrix (m, 0, imat_height (m), 0,
				   imat_width (m), &r, &c);

}

double mat_max_index (mat m, idx_t * r, idx_t * c)
{

  return mat_max_index_submatrix (m, 0, mat_height (m), 0, mat_width (m),
				  r, c);

}

double mat_max_col_index (mat m, idx_t c, idx_t * r)
{

  idx_t i;

  return mat_max_index_submatrix (m, 0, mat_height (m), c, c, r, &i);

}

double mat_max_row_index (mat m, idx_t r, idx_t * c)
{

  idx_t i;

  return mat_max_index_submatrix (m, r, r, 0, mat_width (m), &i, c);

}

double mat_min_index (mat m, idx_t * r, idx_t * c)
{

  return mat_min_index_submatrix (m, 0, mat_height (m), 0, mat_width (m),
				  r, c);

}

double mat_min_col_index (mat m, idx_t c, idx_t * r)
{

  idx_t i;

  return mat_min_index_submatrix (m, 0, mat_height (m), c, c, r, &i);

}

double mat_min_row_index (mat m, idx_t r, idx_t * c)
{

  idx_t i;

  return mat_min_index_submatrix (m, r, r, 0, mat_width (m), &i, c);

}

int imat_max_index (imat m, idx_t * r, idx_t * c)
{

  return imat_max_index_submatrix (m, 0, imat_height (m), 0,
				   imat_width (m), r, c);

}

int imat_max_col_index (imat m, idx_t c, idx_t * r)
{

  idx_t i;

  return imat_max_index_submatrix (m, 0, imat_height (m), c, c, r, &i);

}

int imat_max_row_index (imat m, idx_t r, idx_t * c)
{

  idx_t i;

  return imat_max_index_submatrix (m, r, r, 0, imat_width (m), &i, c);

}

int imat_min_index (imat m, idx_t * r, idx_t * c)
{

  return imat_min_index_submatrix (m, 0, imat_height (m), 0,
				   imat_width (m), r, c);

}

int imat_min_col_index (imat m, idx_t c, idx_t * r)
{

  idx_t i;

  return imat_min_index_submatrix (m, 0, imat_height (m), c, c, r, &i);

}

int imat_min_row_index (imat m, idx_t r, idx_t * c)
{

  idx_t i;

  return imat_min_index_submatrix (m, r, r, 0, imat_width (m), &i, c);

}

/*----------------------------------------------------------------------------*/
double mat_norm_1 (mat m)
{
  idx_t i, j;
  double r = -1, s;

  for (i = 0; i < mat_height (m); i++) {
    s = 0;
    for (j = 0; j < mat_width (m); j++)
      s += fabs (m[i][j]);

    if (s > r)
      r = s;
  }
  return r;
}


double mat_norm_inf (mat m)
{
  idx_t i, j;
  double r = -1, s;

  for (j = 0; j < mat_width (m); j++) {
    s = 0;
    for (i = 0; i < mat_height (m); i++)
      s += fabs (m[i][j]);

    if (s > r)
      r = s;
  }
  return r;
}


/*----------------------------------------------------------------------------*/
double mat_variance (mat m)
{
  idx_t i, j;
  idx_t w, h;
  double sum = 0;
  double var = 0;
  assert (m);
  h = Mat_height (m);
  w = Mat_width (m);
  assert (w * h > 1);		/* otherwise the unbiased variance is not defined */
  for (i = 0; i < h; i++)
    for (j = 0; j < w; j++) {
      sum += m[i][j];
      var += m[i][j] * m[i][j];
    }

  return (var - sum * sum / (w * h)) / (w * h - 1);
}


/*----------------------------------------------------------------------------*/
void mat_normalize (mat m)
{
  idx_t i, j;
  double s = mat_sum (m);
  assert (m);
  for (j = 0; j < Mat_height (m); j++)
    for (i = 0; i < Mat_width (m); i++)
      m[j][i] /= s;
}


/*----------------------------------------------------------------------------*/
/* Column Normalization */
void mat_cols_normalize (mat m, double nr)
{
  idx_t i, j;
  vec  sums = mat_cols_sum (m);
  double vecnorm;

  assert (m);
  for (i = 0; i < Mat_width (m); i++) {

    /* compute the norm */
    vecnorm = 0;
    for (j = 0; j < mat_height (m); j++)
      vecnorm += pow (fabs (m[j][i]), nr);

    vecnorm = pow (vecnorm, 1.0 / nr);

    /* Normalize the corresponding vector */
    for (j = 0; j < mat_height (m); j++)
      m[j][i] /= vecnorm;
  }

  vec_delete (sums);
}


/* Column Normalization */
void mat_rows_normalize (mat m, double nr)
{
  idx_t i;
  assert (m);

  for (i = 0; i < Mat_height (m); i++) 
    vec_normalize (m[i], nr);
}



/*------------------------------------------------------------------------------*/
/** @name Basic arithmetic operations                                           */
/*------------------------------------------------------------------------------*/


/* Operations with a scalar value                                             */
void mat_incr (mat m, double a)
{
  idx_t i, j;
  assert (m);
  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      m[i][j] += a;
}


void mat_decr (mat m, double a)
{
  idx_t i, j;
  assert (m);
  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      m[i][j] -= a;
}


void mat_mul_by (mat m, double a)
{
  idx_t i, j;
  assert (m);
  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      m[i][j] *= a;
}


void mat_div_by (mat m, double a)
{
  idx_t i, j;
  assert (m);
  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      m[i][j] /= a;
}


void imat_incr (imat m, int a)
{
  idx_t i, j;
  assert (m);
  for (i = 0; i < imat_height (m); i++)
    for (j = 0; j < imat_width (m); j++)
      m[i][j] += a;
}


void imat_decr (imat m, int a)
{
  idx_t i, j;
  assert (m);
  for (i = 0; i < imat_height (m); i++)
    for (j = 0; j < imat_width (m); j++)
      m[i][j] -= a;
}


void imat_mul_by (imat m, int a)
{
  idx_t i, j;
  assert (m);
  for (i = 0; i < imat_height (m); i++)
    for (j = 0; j < imat_width (m); j++)
      m[i][j] *= a;
}


void imat_div_by (imat m, int a)
{
  idx_t i, j;
  assert (m);
  for (i = 0; i < imat_height (m); i++)
    for (j = 0; j < imat_width (m); j++)
      m[i][j] /= a;
}


void cmat_incr (cmat m, cplx a)
{
  idx_t i, j;
  assert (m);
  for (i = 0; i < cmat_height (m); i++)
    for (j = 0; j < cmat_width (m); j++)
      m[i][j] = cadd (m[i][j], a);
}


void cmat_decr (cmat m, cplx a)
{
  idx_t i, j;
  assert (m);
  for (i = 0; i < cmat_height (m); i++)
    for (j = 0; j < cmat_width (m); j++)
      m[i][j] = csub (m[i][j], a);
}


void cmat_mul_by (cmat m, cplx a)
{
  idx_t i, j;
  assert (m);
  for (i = 0; i < cmat_height (m); i++)
    for (j = 0; j < cmat_width (m); j++)
      m[i][j] = cmul (m[i][j], a);
}


void cmat_div_by (cmat m, cplx a)
{
  idx_t i, j;
  assert (m);
  for (i = 0; i < cmat_height (m); i++)
    for (j = 0; j < cmat_width (m); j++)
      m[i][j] = cdiv (m[i][j], a);
}


/*------------------------------------------------------------------------------*/
/** @name Basic arithmetic operations on rows/columns                           */
/*------------------------------------------------------------------------------*/

void mat_col_set (mat m, idx_t col, double a)
{
  idx_t i;
  assert (m);
  MAT_END_COL_PARAM (m, col);
  assert (col >= 0 && col < mat_width (m));
  for (i = 0; i < mat_height (m); i++)
    m[i][col] = a;
}


void mat_col_incr (mat m, idx_t col, double a)
{
  idx_t i;
  assert (m);
  MAT_END_COL_PARAM (m, col);
  assert (col >= 0 && col < mat_width (m));
  for (i = 0; i < mat_height (m); i++)
    m[i][col] += a;
}


void mat_col_decr (mat m, idx_t col, double a)
{
  idx_t i;
  assert (m);
  MAT_END_COL_PARAM (m, col);
  assert (col >= 0 && col < mat_width (m));
  for (i = 0; i < mat_height (m); i++)
    m[i][col] -= a;
}


void mat_col_mul_by (mat m, idx_t col, double a)
{
  idx_t i;
  assert (m);
  MAT_END_COL_PARAM (m, col);
  assert (col >= 0 && col < mat_width (m));
  for (i = 0; i < mat_height (m); i++)
    m[i][col] *= a;
}


void mat_col_div_by (mat m, idx_t col, double a)
{
  idx_t i;
  assert (m);
  MAT_END_COL_PARAM (m, col);
  assert (col >= 0 && col < mat_width (m));
  for (i = 0; i < mat_height (m); i++)
    m[i][col] /= a;
}


void imat_col_set (imat m, idx_t col, int a)
{
  idx_t i;
  assert (m);
  MAT_END_COL_PARAM (m, col);
  assert (col >= 0 && col < imat_width (m));
  for (i = 0; i < imat_height (m); i++)
    m[i][col] = a;
}


void imat_col_incr (imat m, idx_t col, int a)
{
  idx_t i;
  assert (m);
  MAT_END_COL_PARAM (m, col);
  assert (col >= 0 && col < imat_width (m));
  for (i = 0; i < imat_height (m); i++)
    m[i][col] += a;
}


void imat_col_decr (imat m, idx_t col, int a)
{
  idx_t i;
  assert (m);
  MAT_END_COL_PARAM (m, col);
  assert (col >= 0 && col < imat_width (m));
  for (i = 0; i < imat_height (m); i++)
    m[i][col] -= a;
}


void imat_col_mul_by (imat m, idx_t col, int a)
{
  idx_t i;
  assert (m);
  MAT_END_COL_PARAM (m, col);
  assert (col >= 0 && col < imat_width (m));
  for (i = 0; i < imat_height (m); i++)
    m[i][col] *= a;
}


void imat_col_div_by (imat m, idx_t col, int a)
{
  idx_t i;
  assert (m);
  MAT_END_COL_PARAM (m, col);
  assert (col >= 0 && col < imat_width (m));
  for (i = 0; i < imat_height (m); i++)
    m[i][col] /= a;
}


void mat_row_incr (mat m, idx_t row, double a)
{
  assert (m);
  MAT_END_ROW_PARAM (m, row);
  vec_incr (m[row], a);
}


void mat_row_decr (mat m, idx_t row, double a)
{
  assert (m);
  MAT_END_ROW_PARAM (m, row);
  vec_decr (m[row], a);
}


void mat_row_mul_by (mat m, idx_t row, double a)
{
  assert (m);
  MAT_END_ROW_PARAM (m, row);
  vec_mul_by (m[row], a);
}


void mat_row_div_by (mat m, idx_t row, double a)
{
  assert (m);
  MAT_END_ROW_PARAM (m, row);
  vec_div_by (m[row], a);
}


void imat_row_incr (imat m, idx_t row, int a)
{
  assert (m);
  MAT_END_ROW_PARAM (m, row);
  ivec_incr (m[row], a);
}


void imat_row_decr (imat m, idx_t row, int a)
{
  assert (m);
  MAT_END_ROW_PARAM (m, row);
  ivec_decr (m[row], a);
}


void imat_row_mul_by (imat m, idx_t row, int a)
{
  assert (m);
  MAT_END_ROW_PARAM (m, row);
  ivec_mul_by (m[row], a);
}


void imat_row_div_by (imat m, idx_t row, int a)
{
  assert (m);
  MAT_END_ROW_PARAM (m, row);
  ivec_div_by (m[row], a);
}



/*------------------------------------------------------------------------------*/
/** @name Component per component operations                                    */
/*------------------------------------------------------------------------------*/

void mat_elem_add (mat m1, mat m2)
{
  idx_t i, j;
  assert (m1 && m2);
  assert (mat_height (m1) == mat_height (m2));
  for (i = 0; i < mat_height (m1); i++)
    for (j = 0; j < mat_width (m1); j++)
      m1[i][j] += m2[i][j];
}


void mat_elem_sub (mat m1, mat m2)
{
  idx_t i, j;
  assert (m1 && m2);
  assert (mat_height (m1) == mat_height (m2));
  for (i = 0; i < mat_height (m1); i++)
    for (j = 0; j < mat_width (m1); j++)
      m1[i][j] -= m2[i][j];
}


void mat_elem_mul (mat m1, mat m2)
{
  idx_t i, j;
  assert (m1 && m2);
  assert (mat_height (m1) == mat_height (m2));
  for (i = 0; i < mat_height (m1); i++)
    for (j = 0; j < mat_width (m1); j++)
      m1[i][j] *= m2[i][j];
}


void mat_elem_div (mat m1, mat m2)
{
  idx_t i, j;
  assert (m1 && m2);
  assert (mat_height (m1) == mat_height (m2));
  for (i = 0; i < mat_height (m1); i++)
    for (j = 0; j < mat_width (m1); j++)
      m1[i][j] /= m2[i][j];
}

void mat_elem_abs (mat m)
{
  int i, j;
  assert(m); 
  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      m[i][j] = fabs (m[i][j]);
}

void mat_elem_exp (mat m)
{
  int i, j;
  assert(m); 
  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      m[i][j] = exp (m[i][j]);
}

void mat_elem_sqrt (mat m)
{
  int i, j;
  assert(m); 
  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      m[i][j] = sqrt (m[i][j]);
}

void mat_elem_sqr (mat m)
{
  int i, j;
  assert(m); 
  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      m[i][j] = m[i][j] * m[i][j];
}

void mat_elem_pow (mat m, double p)
{
  int i, j;
  assert(m);
  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      m[i][j] = pow (m[i][j], p);
}

void imat_elem_add (imat m1, imat m2)
{
  idx_t i, j;
  assert (m1 && m2);
  assert (imat_height (m1) == imat_height (m2));
  for (i = 0; i < imat_height (m1); i++)
    for (j = 0; j < imat_width (m1); j++)
      m1[i][j] += m2[i][j];
}


void imat_elem_sub (imat m1, imat m2)
{
  idx_t i, j;
  assert (m1 && m2);
  assert (imat_height (m1) == imat_height (m2));
  for (i = 0; i < imat_height (m1); i++)
    for (j = 0; j < imat_width (m1); j++)
      m1[i][j] -= m2[i][j];
}


void imat_elem_mul (imat m1, imat m2)
{
  idx_t i, j;
  assert (m1 && m2);
  assert (imat_height (m1) == imat_height (m2));
  for (i = 0; i < imat_height (m1); i++)
    for (j = 0; j < imat_width (m1); j++)
      m1[i][j] *= m2[i][j];
}


void imat_elem_div (imat m1, imat m2)
{
  idx_t i, j;
  assert (m1 && m2);
  assert (imat_height (m1) == imat_height (m2));
  for (i = 0; i < imat_height (m1); i++)
    for (j = 0; j < imat_width (m1); j++)
      m1[i][j] /= m2[i][j];
}


void imat_elem_sqr (imat m)
{
  int i, j;
  assert(m); 
  for (i = 0; i < imat_height (m); i++)
    for (j = 0; j < imat_width (m); j++)
      m[i][j] = m[i][j] * m[i][j];
}

void cmat_elem_add (cmat m1, cmat m2)
{
  idx_t i, j;
  assert (m1 && m2);
  assert (cmat_height (m1) == cmat_height (m2));
  for (i = 0; i < cmat_height (m1); i++)
    for (j = 0; j < cmat_width (m1); j++)
      m1[i][j] = cadd (m1[i][j], m2[i][j]);
}


void cmat_elem_sub (cmat m1, cmat m2)
{
  idx_t i, j;
  assert (m1 && m2);
  assert (cmat_height (m1) == cmat_height (m2));
  for (i = 0; i < cmat_height (m1); i++)
    for (j = 0; j < cmat_width (m1); j++)
      m1[i][j] = csub (m1[i][j], m2[i][j]);
}


void cmat_elem_mul (cmat m1, cmat m2)
{
  idx_t i, j;
  assert (m1 && m2);
  assert (cmat_height (m1) == cmat_height (m2));
  for (i = 0; i < cmat_height (m1); i++)
    for (j = 0; j < cmat_width (m1); j++)
      m1[i][j] = cmul (m1[i][j], m2[i][j]);
}


void cmat_elem_div (cmat m1, cmat m2)
{
  idx_t i, j;
  assert (m1 && m2);
  assert (cmat_height (m1) == cmat_height (m2));
  for (i = 0; i < cmat_height (m1); i++)
    for (j = 0; j < cmat_width (m1); j++)
      m1[i][j] = cdiv (m1[i][j], m2[i][j]);
}


/*------------------------------------------------------------------------------*/
/** @name Matrix arithmetic                                                     */
/*------------------------------------------------------------------------------*/

void mat_add (mat m1, mat m2)
{
  mat_elem_add (m1, m2);
}


void mat_sub (mat m1, mat m2)
{
  mat_elem_sub (m1, m2);
}


void imat_add (imat m1, imat m2)
{
  imat_elem_add (m1, m2);
}


void imat_sub (imat m1, imat m2)
{
  imat_elem_sub (m1, m2);
}


void cmat_add (cmat m1, cmat m2)
{
  cmat_elem_add (m1, m2);
}


void cmat_sub (cmat m1, cmat m2)
{
  cmat_elem_sub (m1, m2);
}


void mat_mul (mat out, mat m1, mat m2)
{
  int i, j, k;
  assert( m1 && m2 );
  it_assert (mat_width(m1)==mat_height(m2), "Size mismatch in mat_mul inputs");
  it_assert (mat_height(m1)==mat_height(out), "Size mismatch between input and output");
  it_assert (mat_width(m2)==mat_width(out), "Size mismatch between input and output");
  
  mat_zeros (out);

  for ( i = 0 ; i< mat_height (m1) ; i++ )
    for ( j = 0 ; j< mat_width (m2) ; j++ )
      for ( k = 0 ; k< mat_width (m1) ; k++ )
	out[i][j] += m1[i][k]*m2[k][j];
}


void imat_mul (imat out, imat m1, imat m2)
{
  int i, j, k;
  assert( m1 && m2 );
  it_assert (imat_width(m1)==imat_height(m2), "Size misimatch in imat_mul inputs");
  it_assert (imat_height(m1)==imat_height(out), "Size mismatch between input and output");
  it_assert (imat_width(m2)==imat_width(out), "Size mismatch between input and output");
  
  imat_zeros (out);

  for ( i = 0 ; i< imat_height (m1) ; i++ )
    for ( j = 0 ; j< imat_width (m2) ; j++ )
      for ( k = 0 ; k< imat_width (m1) ; k++ )
	out[i][j] += m1[i][k]*m2[k][j];
}


void mat_mul_transpose_left (mat out, mat m1, mat m2)
{
  int i, j, k;
  assert( m1 && m2 );
  it_assert (mat_height(m1)==mat_height(m2), "Size mismatch in mat_mul inputs");
  it_assert (mat_width(m1)==mat_height(out), "Size mismatch between input and output");
  it_assert (mat_width(m2)==mat_width(out), "Size mismatch between input and output");
  
  mat_zeros (out);
  
  for ( i= 0; i< mat_width(m1); i++ )
    for ( j= 0; j< mat_width(m2); j++ )
      for ( k= 0; k< mat_height(m1); k++ )
	out[i][j] += m1[k][i]*m2[k][j];
}


void mat_mul_transpose_right (mat out, mat m1, mat m2)
{
  int i, j, k;
  assert( m1 && m2 );
  it_assert (mat_height(m1)==mat_width(m2), "Size mismatch in mat_mul inputs");
  it_assert (mat_width(m1)==mat_height(out), "Size mismatch between input and output");
  it_assert (mat_height(m2)==mat_width(out), "Size mismatch between input and output");
  
  mat_zeros (out);
  
  for ( i= 0; i< mat_height(m1); i++ )
    for ( j= 0; j< mat_height(m2); j++ )
      for ( k= 0; k< mat_width(m1); k++ )
	out[i][j] += m1[i][k]*m2[j][k];
}


void mat_mul_transpose_leftright (mat out, mat m1, mat m2)
{
  int i, j, k;
  assert( m1 && m2 );
  it_assert (mat_width(m1)==mat_width(m2), "Size mismatch in mat_mul inputs");
  it_assert (mat_height(m1)==mat_height(out), "Size mismatch between input and output");
  it_assert (mat_height(m2)==mat_width(out), "Size mismatch between input and output");
  
  mat_zeros (out);

  for ( i= 0; i< mat_width(m1); i++ )
    for ( j= 0; j< mat_height(m2); j++ )
      for ( k= 0; k< mat_height(m1); k++ )
	out[i][j] += m1[k][i]*m2[j][k];
}


/*******************************/

mat mat_new_mul (mat m1, mat m2)
{
  idx_t i, j, k;
  mat  m;
  assert (m1 && m2);
  assert (mat_width (m1) == mat_height (m2));

  m = mat_new_zeros (mat_height (m1), mat_width (m2));

  for (i = 0; i < mat_height (m1); i++)
    for (j = 0; j < mat_width (m2); j++)
      for (k = 0; k < mat_width (m1); k++)
	m[i][j] += m1[i][k] * m2[k][j];
  return m;
}




mat mat_new_add (mat m1, mat m2)
{
  mat  m = mat_clone (m1);
  mat_add (m, m2);
  return m;
}


mat mat_new_sub (mat m1, mat m2)
{
  mat  m = mat_clone (m1);
  mat_add (m, m2);
  return m;
}


imat imat_new_add (imat m1, imat m2)
{
  imat m = imat_clone (m1);
  imat_add (m, m2);
  return m;
}


imat imat_new_sub (imat m1, imat m2)
{
  imat m = imat_clone (m1);
  imat_add (m, m2);
  return m;
}


imat imat_new_mul (imat m1, imat m2)
{
  idx_t i, j, k;
  imat m;
  assert (m1 && m2);
  assert (imat_width (m1) == imat_height (m2));

  m = imat_new_zeros (imat_height (m1), imat_width (m2));

  for (i = 0; i < imat_height (m1); i++)
    for (j = 0; j < imat_width (m2); j++)
      for (k = 0; k < imat_width (m1); k++)
	m[i][j] += m1[i][k] * m2[k][j];
  return m;
}


cmat cmat_new_add (cmat m1, cmat m2)
{
  cmat m = cmat_clone (m1);
  cmat_add (m, m2);
  return m;
}


cmat cmat_new_sub (cmat m1, cmat m2)
{
  cmat m = cmat_clone (m1);
  cmat_add (m, m2);
  return m;
}


cmat cmat_new_mul (cmat m1, cmat m2)
{
  idx_t i, j, k;
  cmat m;
  assert (m1 && m2);
  assert (cmat_width (m1) == cmat_height (m2));

  m = cmat_new_zeros (cmat_height (m1), cmat_width (m2));

  for (i = 0; i < cmat_height (m1); i++)
    for (j = 0; j < cmat_width (m2); j++)
      for (k = 0; k < cmat_width (m1); k++)
	m[i][j] = cadd (m[i][j], cmul (m1[i][k], m2[k][j]));
  return m;
}


vec mat_vec_mul (mat m, vec v)
{
  idx_t i, j;
  vec  r;
  assert (m && v);
  assert (mat_width (m) == vec_length (v));

  r = vec_new_zeros (mat_height (m));

  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      r[i] += m[i][j] * v[j];
  return r;
}


vec imat_vec_mul (imat m, vec v)
{
  idx_t i, j;
  vec  r;
  assert (m && v);
  assert (imat_width (m) == vec_length (v));

  r = vec_new_zeros (imat_height (m));

  for (i = 0; i < imat_height (m); i++)
    for (j = 0; j < imat_width (m); j++)
      r[i] += m[i][j] * v[j];
  return r;
}


ivec imat_bvec_mul (imat m, bvec v)
{
  idx_t i, j;
  ivec r;
  assert (m && v);
  assert (imat_width (m) == bvec_length (v));

  r = ivec_new_zeros (imat_height (m));

  for (i = 0; i < imat_height (m); i++)
    for (j = 0; j < imat_width (m); j++)
      r[i] += m[i][j] * v[j];
  return r;
}


vec mat_ivec_mul (mat m, ivec v)
{
  idx_t i, j;
  vec  r;
  assert (m && v);
  assert (mat_width (m) == ivec_length (v));

  r = vec_new_zeros (mat_height (m));

  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      r[i] += m[i][j] * v[j];
  return r;
}


cvec cmat_vec_mul (cmat m, vec v)
{
  idx_t i, j;
  cvec r;
  assert (m && v);
  assert (cmat_width (m) == vec_length (v));

  r = cvec_new_zeros (cmat_height (m));

  for (i = 0; i < cmat_height (m); i++)
    for (j = 0; j < cmat_width (m); j++)
      r[i] = cadd (r[i], cscale (m[i][j], v[j]));
  return r;
}


cvec cmat_cvec_mul (cmat m, cvec v)
{
  idx_t i, j;
  cvec r;
  assert (m && v);
  assert (cmat_width (m) == cvec_length (v));

  r = cvec_new_zeros (cmat_height (m));

  for (i = 0; i < cmat_height (m); i++)
    for (j = 0; j < cmat_width (m); j++)
      r[i] = cadd (r[i], cmul (m[i][j], v[j]));
  return r;
}


vec vec_mat_mul (vec v, mat m)
{
  idx_t i, j;
  vec  r;
  assert (m && v);
  assert (mat_height (m) == vec_length (v));

  r = vec_new_zeros (mat_width (m));

  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      r[j] += m[i][j] * v[i];
  return r;
}


ivec imat_ivec_mul (imat m, ivec v)
{
  idx_t i, j;
  ivec r;
  assert (m && v);
  assert (imat_width (m) == ivec_length (v));

  r = ivec_new_zeros (imat_height (m));

  for (i = 0; i < imat_height (m); i++)
    for (j = 0; j < imat_width (m); j++)
      r[i] += m[i][j] * v[j];
  return r;
}


ivec ivec_imat_mul (ivec v, imat m)
{
  idx_t i, j;
  ivec r;
  assert (m && v);
  assert (imat_height (m) == ivec_length (v));

  r = ivec_new_zeros (imat_width (m));

  for (i = 0; i < imat_height (m); i++)
    for (j = 0; j < imat_width (m); j++)
      r[j] += m[i][j] * v[i];
  return r;
}


/*----------------------------------------------------------------------------*/
/* Retrieve part or assemble matrices                                         */
mat mat_get_submatrix (mat m, idx_t r1, idx_t c1, idx_t r2, idx_t c2)
{
  mat  s;
  idx_t r;

  if (r2 < r1) {
    it_warning
      ("Upper row number %d is greater than bottom row %d. Returning an empty matrix.",
       r1, r2);
    return mat_new_zeros (0, 0);
  };

  if (c2 < c1) {
    it_warning
      ("Left column %d greater than right column %d. Returning an empty matrix.",
       c1, c2);
    return mat_new_zeros (0, 0);
  };

  MAT_END_ROW_PARAM (m, r1);
  MAT_END_ROW_PARAM (m, r2);
  MAT_END_COL_PARAM (m, c1);
  MAT_END_COL_PARAM (m, c2);

  it_assert (r1 >= 0 && c1 >= 0 && r2 >= 0 && c2 >= 0
	     && r1 < mat_height (m) && r2 < mat_height (m)
	     && c1 < mat_width (m) && c2 < mat_width (m),
	     "Invalid col or row number");

  s = mat_new (r2 - r1 + 1, c2 - c1 + 1);

  for (r = r1; r <= r2; r++)
    vec_init (s[r - r1], m[r] + c1, c2 - c1 + 1);

  return s;
}


imat imat_get_submatrix (imat m, idx_t r1, idx_t c1, idx_t r2, idx_t c2)
{
  imat s;
  idx_t r;

  if (r2 < r1) {
    it_warning
      ("Upper row number %d is greater than bottom row %d . Return void matrix",
       r1, r2);
    return imat_new_zeros (0, 0);
  };

  if (c2 < c1) {
    it_warning
      ("Left column %d greater than right row %d . Return void matrix", c1,
       c2);
    return imat_new_zeros (0, 0);
  };

  MAT_END_ROW_PARAM (m, r1);
  MAT_END_ROW_PARAM (m, r2);
  MAT_END_COL_PARAM (m, c1);
  MAT_END_COL_PARAM (m, c2);

  it_assert (r1 >= 0 && c1 >= 0 && r2 >= 0 && c2 >= 0
	     && r1 < imat_height (m) && r2 < imat_height (m)
	     && c1 < imat_width (m) && c2 < imat_width (m),
	     "Invalid col or row number");

  s = imat_new (r2 - r1 + 1, c2 - c1 + 1);

  for (r = r1; r <= r2; r++)
    ivec_init (s[r - r1], m[r] + c1, c2 - c1 + 1);

  return s;
}


bmat bmat_get_submatrix (bmat m, idx_t r1, idx_t c1, idx_t r2, idx_t c2)
{
  bmat s;
  idx_t r;

  if (r2 < r1) {
    it_warning
      ("Upper row number %d is greater than bottom row %d . Return void matrix",
       r1, r2);
    return bmat_new_zeros (0, 0);
  };

  if (c2 < c1) {
    it_warning
      ("Left column %d greater than right row %d . Return void matrix", c1,
       c2);
    return bmat_new_zeros (0, 0);
  };

  MAT_END_ROW_PARAM (m, r1);
  MAT_END_ROW_PARAM (m, r2);
  MAT_END_COL_PARAM (m, c1);
  MAT_END_COL_PARAM (m, c2);

  it_assert (r1 >= 0 && c1 >= 0 && r2 >= 0 && c2 >= 0
	     && r1 < bmat_height (m) && r2 < bmat_height (m)
	     && c1 < bmat_width (m) && c2 < bmat_width (m),
	     "Invalid col or row number");

  s = bmat_new (r2 - r1 + 1, c2 - c1 + 1);

  for (r = r1; r <= r2; r++)
    bvec_init (s[r - r1], m[r] + c1, c2 - c1 + 1);

  return s;
}


/*----------------------------------------------------------------------------*/
mat mat_set_submatrix (mat m, mat s, idx_t r, idx_t c)
{
  idx_t i, j;

  MAT_END_ROW_PARAM (m, r);
  MAT_END_COL_PARAM (m, c);

  it_assert (r >= 0 && c >= 0 && r < mat_height (m) && c < mat_width (m),
	     "Invalid col or row number");

  it_assert (r + mat_height (s) <= mat_height (m)
	     && c + mat_width (s) <= mat_width (m), "Submatrix is too big");


  for (i = 0; i < mat_height (s); i++)
    for (j = 0; j < mat_width (s); j++)
      m[r + i][c + j] = s[i][j];

  return m;
}


imat imat_set_submatrix (imat m, imat s, idx_t r, idx_t c)
{
  idx_t i, j;

  MAT_END_ROW_PARAM (m, r);
  MAT_END_COL_PARAM (m, c);

  it_assert (r >= 0 && c >= 0 && r < imat_height (m)
	     && c < imat_width (m), "Invalid col or row number");

  it_assert (r + imat_height (s) <= imat_height (m)
	     && c + imat_width (s) <= imat_width (m), "Submatrix is too big");


  for (i = 0; i < imat_height (s); i++)
    for (j = 0; j < imat_width (s); j++)
      m[r + i][c + j] = s[i][j];

  return m;
}


bmat bmat_set_submatrix (bmat m, bmat s, idx_t r, idx_t c)
{
  idx_t i, j;

  MAT_END_ROW_PARAM (m, r);
  MAT_END_COL_PARAM (m, c);

  it_assert (r >= 0 && c >= 0 && r < bmat_height (m)
	     && c < bmat_width (m), "Invalid col or row number");

  it_assert (r + bmat_height (s) <= bmat_height (m)
	     && c + bmat_width (s) <= bmat_width (m), "Submatrix is too big");


  for (i = 0; i < bmat_height (s); i++)
    for (j = 0; j < bmat_width (s); j++)
      m[r + i][c + j] = s[i][j];

  return m;
}


void mat_copy_row (vec v, mat m, idx_t r) 
{
  assert (m);
  MAT_END_ROW_PARAM (m, r);
  assert (r < mat_height (m));
  assert (vec_length (v) == mat_width (m));
  Vec_copy (v, m[r]);
}


void imat_copy_row (ivec v, imat m, idx_t r) {
  assert (m);
  MAT_END_ROW_PARAM (m, r);
  assert (r < imat_height (m));
  assert (ivec_length (v) == imat_width (m));
  Vec_copy (v, m[r]);
}


void bmat_copy_row (bvec v, bmat m, idx_t r) {
  assert (m);
  MAT_END_ROW_PARAM (m, r);
  assert (r < bmat_height (m));
  assert (bvec_length (v) == bmat_width (m));
  Vec_copy (v, m[r]);
}


void cmat_copy_row (cvec v, cmat m, idx_t r) {
  assert (m);
  MAT_END_ROW_PARAM (m, r);
  assert (r < cmat_height (m));
  assert (cvec_length (v) == cmat_width (m));
  Vec_copy (v, m[r]);
}


/*----------------------------------------------------------------------------*/
void mat_copy_col (vec v, mat m, idx_t c)
{
  idx_t i;

  assert (m);
  assert (v);
  MAT_END_COL_PARAM (m, c);
  assert (vec_length (v) == mat_height (m));
  assert (c < mat_width (m));

  for (i = 0; i < mat_height (m); i++)
    v[i] = m[i][c];
}


void imat_copy_col (ivec v, imat m, idx_t c)
{
  idx_t i;

  assert (m);
  assert (v);
  MAT_END_COL_PARAM (m, c);
  assert (ivec_length (v) == imat_height (m));
  assert (c < imat_width (m));

  for (i = 0; i < imat_height (m); i++)
    v[i] = m[i][c];
}


void bmat_copy_col (bvec v, bmat m, idx_t c)
{
  idx_t i;

  assert (m);
  assert (v);
  MAT_END_COL_PARAM (m, c);
  assert (bvec_length (v) == bmat_height (m));
  assert (c < bmat_width (m));

  for (i = 0; i < bmat_height (m); i++)
    v[i] = m[i][c];
}


void cmat_copy_col (cvec v, cmat m, idx_t c)
{
  idx_t i;

  assert (m);
  assert (v);
  MAT_END_COL_PARAM (m, c);
  assert (cvec_length (v) == cmat_height (m));
  assert (c < cmat_width (m));


  for (i = 0; i < cmat_height (m); i++)
    v[i] = m[i][c];
}


vec mat_get_row (mat m, idx_t r) {
  vec  v;
  assert (m);
  v = vec_new (mat_width (m));
  mat_copy_row (v, m, r);
  return v;
}


ivec imat_get_row (imat m, idx_t r) {
  ivec v;
  assert (m);
  v = ivec_new (imat_width (m));
  imat_copy_row (v, m, r);
  return v;
}


bvec bmat_get_row (bmat m, idx_t r) {
  bvec v;
  assert (m);
  v = bvec_new (bmat_width (m));
  bmat_copy_row (v, m, r);
  return v;
}


cvec cmat_get_row (cmat m, idx_t r) {
  cvec v;
  assert (m);
  v = cvec_new (cmat_width (m));
  cmat_copy_row (v, m, r);
  return v;
}

vec mat_get_col (mat m, idx_t c) {
  vec  v;
  assert (m);
  v = vec_new (mat_height (m));
  mat_copy_col (v, m, c);
  return v;
}


ivec imat_get_col (imat m, idx_t c) {
  ivec v;
  assert (m);
  v = ivec_new (imat_height (m));
  imat_copy_col (v, m, c);
  return v;
}


bvec bmat_get_col (bmat m, idx_t c) {
  bvec v;
  assert (m);
  v = bvec_new (bmat_height (m));
  bmat_copy_col (v, m, c);
  return v;
}


cvec cmat_get_col (cmat m, idx_t c) {
  cvec v;
  assert (m);
  v = cvec_new (cmat_height (m));
  cmat_copy_col (v, m, c);
  return v;
}

void Mat_set_col (Mat m, idx_t c, Vec v)
{
  idx_t i;
  size_t s;

  assert (m);
  assert (m[0]);
  assert (v);
  assert (Mat_height (m) == Vec_length (v));
  MAT_END_COL_PARAM (m, c);
  assert (c < Mat_width (m));
  s = Vec_element_size (m[0]);
  assert (s == Vec_element_size (v));

  for (i = 0; i < Mat_height (m); i++)
    memcpy ((char *) m[i] + c * s, (char *) v + s * i, s);
}


void mat_set_col (mat m, idx_t c, vec v)
{
  idx_t i;

  assert (m);
  assert (v);
  assert (mat_height (m) == vec_length (v));
  MAT_END_COL_PARAM (m, c);
  assert (c < mat_width (m));

  for (i = 0; i < mat_height (m); i++)
    m[i][c] = v[i];
}


void imat_set_col (imat m, idx_t c, ivec v)
{
  idx_t i;

  assert (m);
  assert (v);
  assert (imat_height (m) == ivec_length (v));
  MAT_END_COL_PARAM (m, c);
  assert (c < imat_width (m));

  for (i = 0; i < imat_height (m); i++)
    m[i][c] = v[i];
}


void bmat_set_col (bmat m, idx_t c, bvec v)
{
  idx_t i;

  assert (m);
  assert (v);
  assert (bmat_height (m) == bvec_length (v));
  MAT_END_COL_PARAM (m, c);
  assert (c < bmat_width (m));

  for (i = 0; i < bmat_height (m); i++)
    m[i][c] = v[i];
}


void cmat_set_col (cmat m, idx_t c, cvec v)
{
  idx_t i;

  assert (m);
  assert (v);
  assert (cmat_height (m) == cvec_length (v));
  MAT_END_COL_PARAM (m, c);
  assert (c < cmat_width (m));

  for (i = 0; i < cmat_height (m); i++)
    m[i][c] = v[i];
}


/*----------------------------------------------------------------------------*/

vec mat_get_diag (mat m)
{
  int i;  
  int n = __mat_c_my_min (mat_height (m), mat_width (m));
  vec diag = vec_new (n); 
  for (i = 0; i < n; i++)
    diag[i] = m[i][i];

  return diag;
}


ivec imat_get_diag (imat m)
{
  int i;  
  int n = __mat_c_my_min (imat_height (m), imat_width (m));
  ivec diag = ivec_new (n); 
  for (i = 0; i < n; i++)
    diag[i] = m[i][i];

  return diag;
}


bvec bmat_get_diag (bmat m)
{
  int i;  

  int n = __mat_c_my_min (bmat_height (m), bmat_width (m));
  bvec diag = bvec_new (n); 
  for (i = 0; i < n; i++)
    diag[i] = m[i][i];

  return diag;
}


cvec cmat_get_diag (cmat m)
{
  int i;  

  int n = __mat_c_my_min (cmat_height (m), cmat_width (m));
  cvec diag = cvec_new (n); 
  for (i = 0; i < n; i++) {
    diag[i] = m[i][i];
    diag[i] = m[i][i];
  }

  return diag;
}


void mat_diag (mat m, vec v)
{
  idx_t i, mi;
  mi = __mat_c_my_min (mat_height (m), mat_width (m));
  it_assert (mi == vec_length (v), "Incompatible Matrix and vector sizes");

  mat_zeros (m);
  for (i = 0; i < mi; i++)
    m[i][i] = v[i];
}


void imat_diag (imat m, ivec v)
{
  idx_t i, mi;
  mi = __mat_c_my_min (imat_height (m), imat_width (m));
  it_assert (mi == ivec_length (v), "Incompatible Matrix and vector sizes");

  imat_zeros (m);
  for (i = 0; i < mi; i++)
    m[i][i] = v[i];
}


void bmat_diag (bmat m, bvec v)
{
  idx_t i, mi;
  mi = __mat_c_my_min (bmat_height (m), bmat_width (m));
  it_assert (mi == bvec_length (v), "Incompatible Matrix and vector sizes");

  bmat_zeros (m);
  for (i = 0; i < mi; i++)
    m[i][i] = v[i];
}


void cmat_diag (cmat m, cvec v)
{
  idx_t i, mi;
  mi = __mat_c_my_min (cmat_height (m), cmat_width (m));
  it_assert (mi == cvec_length (v), "Incompatible Matrix and vector sizes");

  cmat_zeros (m);
  for (i = 0; i < mi; i++)
    m[i][i] = v[i];
}


void mat_set_diag (mat m, vec v)
{
  idx_t i, mi;
  mi = __mat_c_my_min (mat_height (m), mat_width (m));
  it_assert (mi == vec_length (v), "Incompatible Matrix and vector sizes");

  for (i = 0; i < mi; i++)
    m[i][i] = v[i];
}


void imat_set_diag (imat m, ivec v)
{
  idx_t i, mi;
  mi = __mat_c_my_min (imat_height (m), imat_width (m));
  it_assert (mi == ivec_length (v), "Incompatible Matrix and vector sizes");

  for (i = 0; i < mi; i++)
    m[i][i] = v[i];
}


void bmat_set_diag (bmat m, bvec v)
{
  idx_t i, mi;
  mi = __mat_c_my_min (bmat_height (m), bmat_width (m));
  it_assert (mi == bvec_length (v), "Incompatible Matrix and vector sizes");

  for (i = 0; i < mi; i++)
    m[i][i] = v[i];
}


void cmat_set_diag (cmat m, cvec v)
{
  idx_t i, mi;
  mi = __mat_c_my_min (cmat_height (m), cmat_width (m));
  it_assert (mi == cvec_length (v), "Incompatible Matrix and vector sizes");

  for (i = 0; i < mi; i++)
    m[i][i] = v[i];
}


/*----------------------------------------------------------------------------*/
/** @name Component per component arbitrar function                           */
/*----------------------------------------------------------------------------*/

void mat_apply_function (mat m, it_function_t function, it_args_t args)
{
  idx_t i;
  idx_t l;
  assert (m);
  l = mat_height (m);
  for (i = 0; i < l; i++)
    vec_apply_function (m[i], function, args);
}


mat mat_new_apply_function (mat m, it_function_t function, it_args_t args)
{
  mat  r;
  idx_t i;
  idx_t l;
  assert (m);
  l = mat_height (m);
  r = mat_clone (m);
  for (i = 0; i < l; i++)
    vec_apply_function (r[i], function, args);

  return (r);
}


void imat_apply_function (imat m, it_ifunction_t function, it_args_t args)
{
  idx_t i;
  idx_t l;
  assert (m);
  l = imat_height (m);
  for (i = 0; i < l; i++)
    ivec_apply_function (m[i], function, args);
}


imat imat_new_apply_function (imat m, it_ifunction_t function, it_args_t args)
{
  imat r;
  idx_t i;
  idx_t l;
  assert (m);
  l = imat_height (m);
  r = imat_clone (m);
  for (i = 0; i < l; i++)
    ivec_apply_function (r[i], function, args);

  return (r);
}



/*----------------------------------------------------------------------------*/
/*                Special Matrices                                            */
/*----------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
void mat_rand (mat m)
{
  int  i, j;
  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      m[i][j] = it_rand ();
}


void mat_randn (mat m)
{
  int  i, j;
  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      m[i][j] = it_randn ();
}

/*---------------------------------------------------------------------------*/
mat mat_new_rand (idx_t h, idx_t w)
{
  int  i, j;
  mat  m = mat_new (h, w);

  for (i = 0; i < h; i++)
    for (j = 0; j < w; j++)
      m[i][j] = it_rand ();

  return m;
}


mat mat_new_randn (idx_t h, idx_t w)
{
  int  i, j;
  mat  m = mat_new (h, w);

  for (i = 0; i < h; i++)
    for (j = 0; j < w; j++)
      m[i][j] = it_randn ();

  return m;
}


mat mat_cov (mat m)
{
  int i, j;
  int n; 
  vec vec_temp_1, vec_temp_2;
  mat m_cov; 

  assert(m); 
  assert(mat_height(m)==mat_width(m)); 

  n = mat_height (m);
  m_cov = mat_new(n, n);

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	vec_temp_1 = mat_get_row (m, i);
	vec_temp_2 = mat_get_row (m, j);
	m_cov[i][j] = vec_cov (vec_temp_1, vec_temp_2);
	
	vec_delete (vec_temp_1);
	vec_delete (vec_temp_2);
      }

  return m_cov;
}


