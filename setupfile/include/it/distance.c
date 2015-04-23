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
  Distance measures
  Copyright (C) 2005 Vivien Chappelier, Herve Jegou
*/

#include <it/types.h>
#include <it/distance.h>
#include <it/source_func.h>
#include <math.h>

/*---------------------------------------------------------------------------*/
int vec_distance_hamming (vec v1, vec v2)
{
  int  lmin, lmax, d = 0;
  idx_t i;

  if (vec_length (v1) > vec_length (v2)) {
    lmin = vec_length (v2);
    lmax = vec_length (v1);
  }
  else {
    lmin = vec_length (v1);
    lmax = vec_length (v2);
  }
  for (i = 0; i < lmin; i++)
    if (v1[i] != v2[i])
      d++;
  return d + lmax - lmin;
}


int ivec_distance_hamming (ivec v1, ivec v2)
{
  int  lmin, lmax, d = 0;
  idx_t i;

  if (ivec_length (v1) > ivec_length (v2)) {
    lmin = ivec_length (v2);
    lmax = ivec_length (v1);
  }
  else {
    lmin = ivec_length (v1);
    lmax = ivec_length (v2);
  }
  for (i = 0; i < lmin; i++)
    if (v1[i] != v2[i])
      d++;
  return d + lmax - lmin;
}


int bvec_distance_hamming (bvec v1, bvec v2)
{
  int  lmin, lmax, d = 0;
  idx_t i;

  if (bvec_length (v1) > bvec_length (v2)) {
    lmin = bvec_length (v2);
    lmax = bvec_length (v1);
  }
  else {
    lmin = bvec_length (v1);
    lmax = bvec_length (v2);
  }
  for (i = 0; i < lmin; i++)
    if (v1[i] != v2[i])
      d++;
  return d + lmax - lmin;
}


int cvec_distance_hamming (cvec v1, cvec v2)
{
  int  lmin, lmax, d = 0;
  idx_t i;

  if (cvec_length (v1) > cvec_length (v2)) {
    lmin = cvec_length (v2);
    lmax = cvec_length (v1);
  }
  else {
    lmin = cvec_length (v1);
    lmax = cvec_length (v2);
  }
  for (i = 0; i < lmin; i++)
    if (!ceq (v1[i], v2[i]))
      d++;
  return d + lmax - lmin;
}


/*-----------------------------------------------------------------------*/
double vec_ser (vec v1, vec v2)
{
  int  lmin, lmax, d = 0;
  idx_t i;

  if (vec_length (v1) > vec_length (v2)) {
    lmin = vec_length (v2);
    lmax = vec_length (v1);
  }
  else {
    lmin = vec_length (v1);
    lmax = vec_length (v2);
  }
  for (i = 0; i < lmin; i++)
    if (v1[i] != v2[i])
      d++;

  return (d + lmax - vec_length (v2)) / (double) vec_length (v1);
}


double ivec_ser (ivec v1, ivec v2)
{
  int  lmin, lmax, d = 0;
  idx_t i;

  if (ivec_length (v1) > ivec_length (v2)) {
    lmin = ivec_length (v2);
    lmax = ivec_length (v1);
  }
  else {
    lmin = ivec_length (v1);
    lmax = ivec_length (v2);
  }
  for (i = 0; i < lmin; i++)
    if (v1[i] != v2[i])
      d++;

  return (d + lmax - ivec_length (v2)) / (double) ivec_length (v1);
}


double bvec_ber (bvec v1, bvec v2)
{
  int  lmin, lmax, d = 0;
  idx_t i;

  if (bvec_length (v1) > bvec_length (v2)) {
    lmin = bvec_length (v2);
    lmax = bvec_length (v1);
  }
  else {
    lmin = bvec_length (v1);
    lmax = bvec_length (v2);
  }
  for (i = 0; i < lmin; i++)
    if (v1[i] != v2[i])
      d++;

  return (d + lmax - bvec_length (v2)) / (double) bvec_length (v1);
}



/*-----------------------------------------------------------------------*/
int ivec_distance_levenshtein (ivec v1, ivec v2, int cost_ins, int cost_del,
			       int cost_sub)
{
  /* i1 and i2 respectively index the positions in vectors v1 and v2 */
  idx_t i1, i2;
  int  d;

  imat dislev = imat_new (ivec_length (v1) + 1, ivec_length (v2) + 1);
  imat_set (dislev, 2 * ivec_length (v1));

  /* Initializations */
  for (i1 = 0; i1 < ivec_length (v1) + 1; i1++)
    dislev[i1][0] = i1;

  for (i2 = 0; i2 < ivec_length (v2) + 1; i2++)
    dislev[0][i2] = i2;

  for (i1 = 0; i1 < ivec_length (v1); i1++)
    for (i2 = 0; i2 < ivec_length (v2); i2++) {
      /* First try the substitution */
      if (v1[i1] == v2[i2])
	dislev[i1 + 1][i2 + 1] = dislev[i1][i2];
      else
	dislev[i1 + 1][i2 + 1] = dislev[i1][i2] + cost_sub;

      /* Insertion and deletion costs */
      if (dislev[i1 + 1][i2] + cost_del < dislev[i1 + 1][i2 + 1])
	dislev[i1 + 1][i2 + 1] = dislev[i1 + 1][i2] + cost_del;

      if (dislev[i1][i2 + 1] + cost_ins < dislev[i1 + 1][i2 + 1])
	dislev[i1 + 1][i2 + 1] = dislev[i1][i2 + 1] + cost_ins;
    }

  d = dislev[ivec_length (v1)][ivec_length (v2)];
  imat_delete (dislev);
  return d;
}

/*-----------------------------------------------------------------------*/
double vec_distance_norm (vec v1, vec v2, double norm)
{
  double s = 0;
  idx_t i;
  assert (vec_length (v1) == vec_length (v2));

  for (i = 0; i < vec_length (v1); i++)
    s += pow (fabs (v1[i] - v2[i]), norm);

  return pow (s, 1. / (double) norm);
}


double mat_distance_norm (mat m1, mat m2, double norm)
{
  double s = 0;
  idx_t i, j;
  assert (mat_width (m1) == mat_width (m2));
  assert (mat_height (m1) == mat_height (m2));

  for (i = 0; i < mat_height (m1); i++)
    for (j = 0; j < mat_width (m1); j++)
      s += pow (fabs (m1[i][j] - m2[i][j]), norm);

  return pow (s, 1. / (double) norm);
}


int ivec_distance_norm1 (ivec v1, ivec v2)
{
  int s = 0;
  idx_t i;
  assert (ivec_length (v1) == ivec_length (v2));

  for (i = 0; i < ivec_length (v1); i++)
    s += abs (v1[i] - v2[i]);

  return s;
}


/*-----------------------------------------------------------------------*/
double vec_distance_mse (vec v1, vec v2, double rec_value)
{
  double s = 0, d;
  idx_t i;
  vec  vmin, vmax;		/* vmin is the shortest vector and vmax the longer one */

  if (vec_length (v1) > vec_length (v2)) {
    vmin = v2;
    vmax = v1;
  }
  else {
    vmin = v1;
    vmax = v2;
  }

  for (i = 0; i < vec_length (vmin); i++) {
    d = vmax[i] - vmin[i];
    s += d * d;
  }

  for (; i < vec_length (vmax); i++) {
    d = vmax[i] - rec_value;
    s += d * d;
  }
  return s / vec_length (vmax);
}


/*-----------------------------------------------------------------------*/
double mat_distance_mse (mat m1, mat m2, double rec_value)
{
  double s = 0;
  idx_t i;
  mat  mmin, mmax;		/* mmin is the shortest matrix and mmax the tallest one */

  if (mat_height (m1) > mat_height (m2)) {
    mmin = m2;
    mmax = m1;
  }
  else {
    mmin = m1;
    mmax = m2;
  }

  for (i = 0; i < mat_height (mmin); i++)
    s += vec_distance_mse (mmax[i], mmin[i], rec_value);

  for (; i < mat_height (mmax); i++)
    s += vec_distance_mse (mmax[i], vec_null, rec_value);

  return s / mat_height (mmax);
}

double ivec_distance_mse (ivec v1, ivec v2, double rec_value)
{
  double s = 0, d;
  idx_t i;
  ivec vmin, vmax;		/* vmin is the shortest vector and vmax the longer one */

  if (ivec_length (v1) > ivec_length (v2)) {
    vmin = v2;
    vmax = v1;
  }
  else {
    vmin = v1;
    vmax = v2;
  }

  for (i = 0; i < ivec_length (vmin); i++) {
    d = vmax[i] - vmin[i];
    s += d * d;
  }

  for (; i < ivec_length (vmax); i++) {
    d = vmax[i] - rec_value;
    s += d * d;
  }
  return s / (double) ivec_length (vmax);
}

double imat_distance_mse (imat m1, imat m2, double rec_value)
{
  double s = 0;
  idx_t i;
  imat mmin, mmax;		/* mmin is the shortest matrix and mmax the tallest one */

  if (imat_height (m1) > imat_height (m2)) {
    mmin = m2;
    mmax = m1;
  }
  else {
    mmin = m1;
    mmax = m2;
  }

  for (i = 0; i < imat_height (mmin); i++)
    s += ivec_distance_mse (mmax[i], mmin[i], rec_value);

  for (; i < imat_height (mmax); i++)
    s += ivec_distance_mse (mmax[i], ivec_null, rec_value);

  return s / imat_height (mmax);
}


long ivec_distance_sqr (ivec v1, ivec v2)
{
  assert (ivec_length (v1) == ivec_length (v2));

  long s = 0, d;
  idx_t i;

  for (i = 0; i < ivec_length (v1); i++) {
    d = v2[i] - v1[i];
    s += d * (long) d;
  }

  return s;
}


/*-----------------------------------------------------------------------*/
/* Return the Kullback-Leibler pseudo-distance between distribution pdf1 and pdf2. */
double vec_distance_kullback_leibler (vec pdf1, vec pdf2)
{
  idx_t i;
  double d = 0;
  assert (vec_length (pdf1) == vec_length (pdf2));
  assert (is_valid_pdf (pdf1, 1e-10) && is_valid_pdf (pdf2, 1e-10));


  for (i = 0; i < vec_length (pdf1); i++) {
    if (pdf1[i] != 0) {
      if (pdf2[i] == 0)
	return INT_MAX;
      else
	d += pdf1[i] * log (pdf1[i] / pdf2[i]);
    }
  }
  return d / log (2);
}




/*-----------------------------------------------------------------------*/
/* compute the matrix of distance associated with a given vector         */
mat compute_distance_matrix (mat v, double nr)
{
  int n = mat_height (v), i, j;

  mat dis = mat_new (n, n);

  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < n ; j++) {
      dis[i][j] = vec_distance_norm (v[i], v[j], nr);
    }

  return dis;
}


/*----------------------------------------------------------------*/
/*             Distances between sparse vectors                   */

/* compute the L1 distance between two sparse vectors */
double spvec_distance_norm1 (ivec svi1, vec sv1, ivec svi2, vec sv2)
{
  double s = 0;
  idx_t i1 = 0, i2 = 0;

  assert (ivec_length (svi1) == vec_length (sv1));
  assert (ivec_length (svi2) == vec_length (sv2));

  while (1) {

    if (i1 == ivec_length (svi1)) {
      while (i2 < ivec_length (svi2)) 
	s += fabs (sv2[i2++]);
      break;
    }

    if (i2 == ivec_length (svi2)) {
      while (i1 < ivec_length (svi1)) 
	s += fabs (sv1[i1++]);
      break;
    }

    if (svi1[i1] == svi2[i2]) {
      s += fabs (sv1[i1] - sv2[i2]);
      i1++;
      i2++;
    }

    else {
      if (svi1[i1] < svi2[i2]) 
	s += fabs (sv1[i1++]);
      else
	s += fabs (sv2[i2++]);
    }
  }

  return s;
}


int spivec_distance_norm1 (ivec svi1, ivec sv1, ivec svi2, ivec sv2)
{
  int s = 0;
  idx_t i1 = 0, i2 = 0;

  assert (ivec_length (svi1) == ivec_length (sv1));
  assert (ivec_length (svi2) == ivec_length (sv2));

  while (1) {

    if (i1 == ivec_length (svi1)) {
      while (i2 < ivec_length (svi2)) 
	s += abs (sv2[i2++]);
      break;
    }

    if (i2 == ivec_length (svi2)) {
      while (i1 < ivec_length (svi1)) 
	s += abs (sv1[i1++]);
      break;
    }

    if (svi1[i1] == svi2[i2]) {
      s += abs (sv1[i1] - sv2[i2]);
      i1++;
      i2++;
    }

    else {
      if (svi1[i1] < svi2[i2]) 
	s += abs (sv1[i1++]);
      else
	s += abs (sv2[i2++]);
    }
  }

  return s;
}


/* compute the square of the euclidean distance between two sparse vectors */
double spvec_distance_sqr (ivec svi1, vec sv1, ivec svi2, vec sv2)
{
  double s = 0;
  idx_t i1 = 0, i2 = 0;

  assert (ivec_length (svi1) == vec_length (sv1));
  assert (ivec_length (svi2) == vec_length (sv2));

  while (1) {

    if (i1 == ivec_length (svi1)) {
      while (i2 < ivec_length (svi2)) {
	s += sv2[i2] * sv2[i2];
	i2++;
      }
      break;
    }

    if (i2 == ivec_length (svi2)) {
      while (i1 < ivec_length (svi1)) {
	s += sv1[i1] * sv1[i1];
	i1++;
      }
      break;
    }

    if (svi1[i1] == svi2[i2]) {
      s += (sv1[i1] - sv2[i2]) * (sv1[i1] - sv2[i2]);
      i1++;
      i2++;
    }

    else {
      if (svi1[i1] < svi2[i2]) {
	s += sv1[i1] * sv1[i1];
	i1++;
      }
      else {
	s += sv2[i2] * sv2[i2];
	i2++;
      }
    }
  }

  return s;
}


int spivec_distance_sqr (ivec svi1, ivec sv1, ivec svi2, ivec sv2)
{
  double s = 0;
  idx_t i1 = 0, i2 = 0;

  assert (ivec_length (svi1) == ivec_length (sv1));
  assert (ivec_length (svi2) == ivec_length (sv2));

  while (1) {

    if (i1 == ivec_length (svi1)) {
      while (i2 < ivec_length (svi2)) {
	s += sv2[i2] * sv2[i2];
	i2++;
      }
      break;
    }

    if (i2 == ivec_length (svi2)) {
      while (i1 < ivec_length (svi1)) {
	s += sv1[i1] * sv1[i1];
	i1++;
      }
      break;
    }

    if (svi1[i1] == svi2[i2]) {
      s += (sv1[i1] - sv2[i2]) * (sv1[i1] - sv2[i2]);
      i1++;
      i2++;
    }

    else {
      if (svi1[i1] < svi2[i2]) {
	s += sv1[i1] * sv1[i1];
	i1++;
      }
      else {
	s += sv2[i2] * sv2[i2];
	i2++;
      }
    }
  }

  return s;
}


double spvec_distance_norm2 (ivec svi1, vec sv1, ivec svi2, vec sv2)
{
  return sqrt (spvec_distance_sqr (svi1, sv1, svi2, sv2));
}
