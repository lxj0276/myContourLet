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
  Linear algebra functions
  Copyright (C) 2006 Vivien Chappelier, Herve Jegou, François Cayre
*/

#include <it/linalg.h>
#include <it/mat.h>
#include <it/vec.h>

/* 
  Static helper functions: 
  - ccdiv( ... ): I know it is ugly. Should be performed 
    with standard libit cdiv. Too lazy to do it right 
    now. 
  - maxx( double, double ): ugly too. 

-> Francois, tu seras flagellé sur la place publique pour cette insulte 
   au dogme que tu es censé défendre
*/

static void ccdiv (double xr, double xi, double yr, double yi, double *cdivr,
		   double *cdivi)
{
  double r, d;

  if (fabs (yr) > fabs (yi)) {
    r = yi / yr;
    d = yr + r * yi;
    *cdivr = (xr + r * xi) / d;
    *cdivi = (xi - r * xr) / d;
  }
  else {
    r = yr / yi;
    d = yi + r * yr;
    *cdivr = (r * xr + xi) / d;
    *cdivi = (r * xi - xr) / d;
  }
}


static double maxx (double a, double b)
{
  if (a > b)
    return (a);
  else
    return (b);
}

int mat_is_symmetric( mat a ) 
{

  idx_t i, j; 

  if( mat_width(a)!=mat_height(a) ) 
    return 0; 

  for ( i= 0; i< mat_width(a); i++ )
    for ( j= 0; j< i; j++ ) 
      if ( a[i][j]!=a[j][i] ) return 0; 

  return 1; 
}

/* 
   SVD decomposition 

   Set either U or V (or both) to  NULL to  *not* compute 
   singular subspaces. 

   If one wants to compute singular subspaces, make  sure 
   allocation of U and V matrices has been performed in a 
   proper way (check size). 

   In particular: 
   - U = mat_new( m, min( m, n ) ); 
   - V = mat_new( n, n ); 

   NOTE: This implementation is known to work with m>=n 
   matrices and *may* fail with m<n matrices. We chose 
   the conservative approach. 
   
*/

vec mat_svd (mat M, mat U, mat V)
{
  /* Initialize */

  int  i, j, k, kase, ks, m, n, p, pp, iter, nu, nct, nrt, wantu = 0, wantv =
    0;
  int  mmax;

  double t, f, g, cs, sn, eps, scale, sp, spm1, epm1, sk, ek, b, c, shift;
  mat A; 
  vec  s, e, work;

  m = mat_height (M);
  n = mat_width (M);

  /* Comment this next line if you want a non conservative SVD */
  it_assert (m >= n, "This SVD is meant for m > n matrices");

  nu = m > n ? n : m;

  s = vec_new_zeros (m + 1 > n ? n : m + 1);

  if (U) {
    it_assert (mat_height (U) == m
	       && mat_width (U) == nu,
	       "U matrix should be a m by min(m+1,n) real matrix");
    mat_zeros (U);
    wantu = 1;
  }

  if (V) {
    it_assert (mat_height (V) == mat_width (V), "Matrix V must be square");
    it_assert (mat_height (V) == n, "Size mismatch for matrix V");
    mat_zeros (V);
    wantv = 1;
  }

  A = mat_clone( M );

  e = vec_new_zeros (n);
  work = vec_new_zeros (m);

  /* 
     Reduce A to bidiagonal form, storing the diagonal elements
     in s and the super-diagonal elements in e.
   */

  nct = m - 1 > n ? n : m - 1;
  nrt = 0 > (n - 2 > m ? m : n - 2) ? 0 : (n - 2 > m ? m : n - 2);

  mmax = nct > nrt ? nct : nrt;

  for (k = 0; k < mmax; k++) {
    if (k < nct) {
      /* 
         Compute the transformation for the k-th column and
         place the k-th diagonal in s[k].
         Compute 2-norm of k-th column without under/overflow.
       */

      s[k] = 0;
      for (i = k; i < m; i++)
	s[k] = sqrt (s[k] * s[k] + A[i][k] * A[i][k]);

      if (s[k]) {
	if (A[k][k] < 0.0)
	  s[k] = -s[k];

	for (i = k; i < m; i++)
	  A[i][k] /= s[k];

	A[k][k] += 1.0;
      }
      s[k] = -s[k];
    }
    for (j = k + 1; j < n; j++) {
      if ((k < nct) && (s[k] != 0.0)) {
	/* Apply the transformation */

	t = 0.;
	for (i = k; i < m; i++)
	  t += A[i][k] * A[i][j];

	t = -t / A[k][k];
	for (i = k; i < m; i++)
	  A[i][j] += t * A[i][k];

      }

      /* 
         Place the k-th row of A into e for the
         subsequent calculation of the row transformation
       */

      e[j] = A[k][j];
    }
    if (wantu && (k < nct)) {
      /* 
         Place the transformation in U for subsequent back
         multiplication.
       */

      for (i = k; i < m; i++)
	U[i][k] = A[i][k];

    }
    if (k < nrt) {
      /* 
         Compute the k-th row transformation and place the
         k-th super-diagonal in e[k].
         Compute 2-norm without under/overflow.
       */

      e[k] = 0;
      for (i = k + 1; i < n; i++)
	e[k] = sqrt (e[k] * e[k] + e[i] * e[i]);

      if (e[k]) {
	if (e[k + 1] < 0.0)
	  e[k] = -e[k];

	for (i = k + 1; i < n; i++)
	  e[i] /= e[k];

	e[k + 1] += 1.0;
      }
      e[k] = -e[k];
      if ((k + 1 < m) && e[k]) {

	/* Apply the transformation */

	for (i = k + 1; i < m; i++)
	  work[i] = 0.0;

	for (j = k + 1; j < n; j++)
	  for (i = k + 1; i < m; i++)
	    work[i] += e[j] * A[i][j];

	for (j = k + 1; j < n; j++) {
	  t = -e[j] / e[k + 1];
	  for (i = k + 1; i < m; i++)
	    A[i][j] += t * work[i];
	}
      }
      if (wantv) {
	/* 
	   Place the transformation in V for subsequent
	   back multiplication
	 */

	for (i = k + 1; i < n; i++)
	  V[i][k] = e[i];

      }
    }
  }

  /* Set up the final bidiagonal matrix or order p */

  p = n > m + 1 ? m + 1 : n;

  if (nct < n)
    s[nct] = A[nct][nct];
  if (m < p)
    s[p - 1] = 0.0;
  if (nrt + 1 < p)
    e[nrt] = A[nrt][p - 1];

  e[p - 1] = 0.0;

  /* If required, generate U */

  if (wantu) {
    for (j = nct; j < nu; j++) {
      for (i = 0; i < m; i++)
	U[i][j] = 0.0;

      U[j][j] = 1.0;
    }
    for (k = nct - 1; k >= 0; k--) {
      if (s[k]) {
	for (j = k + 1; j < nu; j++) {
	  t = 0.;
	  for (i = k; i < m; i++)
	    t += U[i][k] * U[i][j];

	  t = -t / U[k][k];
	  for (i = k; i < m; i++)
	    U[i][j] += t * U[i][k];

	}
	for (i = k; i < m; i++)
	  U[i][k] = -U[i][k];

	U[k][k] = 1.0 + U[k][k];
	for (i = 0; i < k - 1; i++)
	  U[i][k] = 0.0;

      }
      else {
	for (i = 0; i < m; i++)
	  U[i][k] = 0.0;

	U[k][k] = 1.0;
      }
    }
  }

  /* If required, generate V */

  if (wantv) {
    for (k = n - 1; k >= 0; k--) {
      if ((k < nrt) && e[k]) {
	for (j = k + 1; j < nu; j++) {
	  t = 0.;
	  for (i = k + 1; i < n; i++)
	    t += V[i][k] * V[i][j];

	  t = -t / V[k + 1][k];
	  for (i = k + 1; i < n; i++)
	    V[i][j] += t * V[i][k];

	}
      }
      for (i = 0; i < n; i++)
	V[i][k] = 0.0;

      V[k][k] = 1.0;
    }
  }


  /* Main iteration loop for the singular values */

  pp = p - 1;
  iter = 0;
  eps = pow (2.0, -52.0);

  while (p > 0) {
    /* Here is where a test for too many iterations would go */

    /* 
       This section of the program inspects for
       negligible elements in the s and e arrays.  On
       completion the variables kase and k are set as follows.

       kase = 1     if s(p) and e[k-1] are negligible and k<p
       kase = 2     if s(k) is negligible and k<p
       kase = 3     if e[k-1] is negligible, k<p, and
       s(k), ..., s(p) are not negligible (qr step).
       kase = 4     if e(p-1) is negligible (convergence).
     */

    for (k = p - 2; k >= -1; k--) {
      if (k == -1)
	break;

      if (fabs (e[k]) <= eps * (fabs (s[k]) + fabs (s[k + 1]))) {
	e[k] = 0.0;
	break;
      }
    }

    if (k == p - 2)
      kase = 4;

    else {
      for (ks = p - 1; ks >= k; ks--) {
	if (ks == k)
	  break;

	t =
	  (ks != p ? fabs (e[ks]) : 0.) + (ks !=
					   k + 1 ? fabs (e[ks - 1]) : 0.);
	if (fabs (s[ks]) <= eps * t) {
	  s[ks] = 0.0;
	  break;
	}
      }
      if (ks == k)
	kase = 3;
      else if (ks == p - 1)
	kase = 1;
      else {
	kase = 2;
	k = ks;
      }
    }
    k++;

    /* Perform the task indicated by kase */

    switch (kase) {
      /* Deflate negligible s(p) */

    case 1:{
	f = e[p - 2];
	e[p - 2] = 0.;
	for (j = p - 2; j >= k; j--) {
	  t = sqrt (s[j] * s[j] + f * f);
	  cs = s[j] / t;
	  sn = f / t;
	  s[j] = t;
	  if (j != k) {
	    f = -sn * e[j - 1];
	    e[j - 1] = cs * e[j - 1];
	  }
	  if (wantv) {
	    for (i = 0; i < n; i++) {
	      t = cs * V[i][j] + sn * V[i][p - 1];
	      V[i][p - 1] = -sn * V[i][j] + cs * V[i][p - 1];
	      V[i][j] = t;
	    }
	  }
	}
      }
      break;

      /* Split at negligible s(k) */

    case 2:{
	f = e[k - 1];
	e[k - 1] = 0.0;
	for (j = k; j < p; j++) {
	  t = sqrt (s[j] * s[j] + f * f);
	  cs = s[j] / t;
	  sn = f / t;
	  s[j] = t;
	  f = -sn * e[j];
	  e[j] = cs * e[j];
	  if (wantu) {
	    for (i = 0; i < m; i++) {
	      t = cs * U[i][j] + sn * U[i][k - 1];
	      U[i][k - 1] = -sn * U[i][j] + cs * U[i][k - 1];
	      U[i][j] = t;
	    }
	  }
	}
      }
      break;

      /* Perform one qr step */

    case 3:{

	/* Calculate the shift */
	scale =
	  maxx (maxx
		(maxx
		 (maxx (fabs (s[p - 1]), fabs (s[p - 2])), fabs (e[p - 2])),
		 fabs (s[k])), fabs (e[k]));

	sp = s[p - 1] / scale;
	spm1 = s[p - 2] / scale;
	epm1 = e[p - 2] / scale;
	sk = s[k] / scale;
	ek = e[k] / scale;
	b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
	c = (sp * epm1) * (sp * epm1);
	shift = 0.;
	if (b || c) {
	  shift = sqrt (b * b + c);
	  if (b < 0.0)
	    shift = -shift;

	  shift = c / (b + shift);
	}
	f = (sk + sp) * (sk - sp) + shift;
	g = sk * ek;

	/* Chase zeros */

	for (j = k; j < p - 1; j++) {
	  t = sqrt (f * f + g * g);
	  cs = f / t;
	  sn = g / t;
	  if (j != k)
	    e[j - 1] = t;

	  f = cs * s[j] + sn * e[j];
	  e[j] = cs * e[j] - sn * s[j];
	  g = sn * s[j + 1];
	  s[j + 1] = cs * s[j + 1];
	  if (wantv) {
	    for (i = 0; i < n; i++) {
	      t = cs * V[i][j] + sn * V[i][j + 1];
	      V[i][j + 1] = -sn * V[i][j] + cs * V[i][j + 1];
	      V[i][j] = t;
	    }
	  }
	  t = sqrt (f * f + g * g);
	  cs = f / t;
	  sn = g / t;
	  s[j] = t;
	  f = cs * e[j] + sn * s[j + 1];
	  s[j + 1] = -sn * e[j] + cs * s[j + 1];
	  g = sn * e[j + 1];
	  e[j + 1] = cs * e[j + 1];
	  if (wantu && (j < m - 1)) {
	    for (i = 0; i < m; i++) {
	      t = cs * U[i][j] + sn * U[i][j + 1];
	      U[i][j + 1] = -sn * U[i][j] + cs * U[i][j + 1];
	      U[i][j] = t;
	    }
	  }
	}
	e[p - 2] = f;
	iter = iter + 1;
      }

      break;

      /* Convergence */

    case 4:{

	/* Make the singular values positive */

	if (s[k] <= 0.0) {
	  s[k] = (s[k] < 0.0 ? -s[k] : 0.0);
	  if (wantv)
	    for (i = 0; i <= pp; i++)
	      V[i][k] = -V[i][k];
	}

	/* Order the singular values */

	while (k < pp) {
	  if (s[k] >= s[k + 1])
	    break;

	  t = s[k];
	  s[k] = s[k + 1];
	  s[k + 1] = t;

	  if (wantv && (k < n - 1)) {
	    for (i = 0; i < n; i++)
	      {
		t = V[i][k + 1];
		V[i][k + 1] = V[i][k];
		V[i][k] = t;
	      }
	  }

	  if (wantu && (k < m - 1)) {
	    for (i = 0; i < m; i++)
	      {
		t = U[i][k + 1];
		U[i][k + 1] = U[i][k];
		U[i][k] = t;
	      }
	  }

	  k++;
	}
	iter = 0;
	p--;
      }
      break;
    }
  }

  vec_delete (e);
  vec_delete (work);

  return (s);
}


/* 
   This is the recommended way of computing the rank of a matrix. 
 */

int mat_rank (mat a)
{
  int  i, r = 0;
  double eps = pow (2., -52.);
  double tol;
  vec  s = mat_svd (a, NULL, NULL);

  tol = (double) (mat_height (a) >
		  mat_width (a) ? mat_height (a) : mat_width (a)) * s[0] * eps;

  for (i = 0; i < vec_length (s); i++)
    if (s[i] > tol)
      r++;

  vec_delete (s);

  return (r);
}


/*
  Condition number of a matrix is the ratio of the two extremal singular values. 
 */
double mat_cond (mat a)
{
  double cond;
  vec  s = mat_svd (a, NULL, NULL);

  cond =
    s[0] / s[mat_height (a) >
	     mat_width (a) ? mat_width (a) - 1 : mat_height (a) - 1];

  vec_delete (s);

  return (cond);

}

/*
  2-norm of a matrix is the largest singular value. 
 */

double mat_norm2 (mat a)
{

  double norm2;
  vec  s = mat_svd (a, NULL, NULL);

  norm2 = s[0];

  vec_delete (s);

  return (norm2);

}

static mat mat_tridiag (vec e, vec d, mat evec)
{

  int  i, j, k, n;
  double scale, f, g, h, hh;

  it_assert (mat_width (evec) == mat_height (evec), "Matrix must be square");

  n = mat_height (evec);

  vec_zeros (e);

  for (j = 0; j < n; j++)
    d[j] = evec[n - 1][j];

  /* Householder reduction to tridiagonal form */

  for (i = n - 1; i > 0; i--) {
    /* Scale to avoid under/overflow */

    scale = 0.0;
    h = 0.0;

    for (k = 0; k < i; k++)
      scale += fabs (d[k]);

    if (!scale) {
      e[i] = d[i - 1];
      for (j = 0; j < i; j++) {
	d[j] = evec[i - 1][j];
	evec[i][j] = 0.0;
	evec[j][i] = 0.0;
      }
    }
    else {
      /* Generate Householder vector */

      for (k = 0; k < i; k++) {
	d[k] /= scale;
	h += d[k] * d[k];
      }
      f = d[i - 1];
      g = sqrt (h);
      if (f > 0)
	g = -g;

      e[i] = scale * g;
      h = h - f * g;
      d[i - 1] = f - g;
      for (j = 0; j < i; j++)
	e[j] = 0.0;

      /* Apply similarity transformation to remaining columns */

      for (j = 0; j < i; j++) {
	f = d[j];
	evec[j][i] = f;
	g = e[j] + evec[j][j] * f;
	for (k = j + 1; k <= i - 1; k++) {
	  g += evec[k][j] * d[k];
	  e[k] += evec[k][j] * f;
	}
	e[j] = g;
      }
      f = 0.0;
      for (j = 0; j < i; j++) {
	e[j] /= h;
	f += e[j] * d[j];
      }
      hh = f / (h + h);
      for (j = 0; j < i; j++)
	e[j] -= hh * d[j];

      for (j = 0; j < i; j++) {
	f = d[j];
	g = e[j];
	for (k = j; k <= i - 1; k++)
	  evec[k][j] -= (f * e[k] + g * d[k]);

	d[j] = evec[i - 1][j];
	evec[i][j] = 0.0;
      }
    }
    d[i] = h;
  }

  /* Accumulate transformations */

  for (i = 0; i < n - 1; i++) {
    evec[n - 1][i] = evec[i][i];
    evec[i][i] = 1.0;
    h = d[i + 1];
    if (h != 0.0) {
      for (k = 0; k <= i; k++)
	d[k] = evec[k][i + 1] / h;

      for (j = 0; j <= i; j++) {
	g = 0.0;
	for (k = 0; k <= i; k++)
	  g += evec[k][i + 1] * evec[k][j];

	for (k = 0; k <= i; k++)
	  evec[k][j] -= g * d[k];

      }
    }

    for (k = 0; k <= i; k++)
      evec[k][i + 1] = 0.0;

  }
  for (j = 0; j < n; j++) {
    d[j] = evec[n - 1][j];
    evec[n - 1][j] = 0.0;
  }
  evec[n - 1][n - 1] = 1.0;
  e[0] = 0.0;

  return (evec);

}


/* Symmetric tridiagonal QL algorithm */
static vec mat_tridiag_ql (vec e, vec d, mat evec)
{
  int  i, j, k, l, m, n, iter;
  double f, tst1, eps, p, g, r, el1, dl1, h, c, c2, c3, s, s2;

  n = mat_height (evec);

  for (i = 1; i < n; i++)
    e[i - 1] = e[i];

  e[n - 1] = 0.0;

  f = 0.0;
  tst1 = 0.0;
  eps = pow (2.0, -52.0);

  for (l = 0; l < n; l++) {
    /* Find small subdiagonal element */

    tst1 =
      tst1 > fabs (d[l]) + fabs (e[l]) ? tst1 : fabs (d[l]) + fabs (e[l]);
    m = l;

    while (m < n) {
      if (fabs (e[m]) <= eps * tst1)
	break;

      m++;
    }

    /* 
       If m == l, d[l] is an eigenvalue,
       otherwise, iterate.
     */

    if (m > l) {
      iter = 0;

      do {
	iter++;			/* (Could check iteration count here.) */

	/* Compute implicit shift */

	g = d[l];
	p = (d[l + 1] - g) / (2.0 * e[l]);
	r = sqrt (p * p + 1.);
	if (p < 0)
	  r = -r;

	d[l] = e[l] / (p + r);
	d[l + 1] = e[l] * (p + r);
	dl1 = d[l + 1];
	h = g - d[l];
	for (i = l + 2; i < n; i++)
	  d[i] -= h;

	f = f + h;

	/* Implicit QL transformation */

	p = d[m];
	c = 1.0;
	c2 = c;
	c3 = c;
	el1 = e[l + 1];
	s = 0.0;
	s2 = 0.0;

	for (i = m - 1; i >= l; i--) {
	  c3 = c2;
	  c2 = c;
	  s2 = s;
	  g = c * e[i];
	  h = c * p;
	  r = sqrt (p * p + e[i] * e[i]);
	  e[i + 1] = s * r;
	  s = e[i] / r;
	  c = p / r;
	  p = c * d[i] - s * g;
	  d[i + 1] = h + s * (c * g + s * d[i]);

	  /* Accumulate transformation */

	  for (k = 0; k < n; k++) {
	    h = evec[k][i + 1];
	    evec[k][i + 1] = s * evec[k][i] + c * h;
	    evec[k][i] = c * evec[k][i] - s * h;
	  }
	}
	p = -s * s2 * c3 * el1 * e[l] / dl1;
	e[l] = s * p;
	d[l] = c * p;

	/* Check for convergence */

      } while (fabs (e[l]) > eps * tst1);
    }
    d[l] = d[l] + f;
    e[l] = 0.0;
  }

  /* Sort eigenvalues and corresponding vectors */

  for (i = 0; i < n - 1; i++) {
    k = i;
    p = d[i];

    for (j = i + 1; j < n; j++)
      if (d[j] < p) {
	k = j;
	p = d[j];
      }

    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (j = 0; j < n; j++) {
	p = evec[j][i];
	evec[j][i] = evec[j][k];
	evec[j][k] = p;
      }
    }
  }

  return (d);
}


/*  
    I doubt there's any use for this function. Could be 
    declared static. 
*/

mat mat_hessenberg (mat a, mat V)
{

  int  i, j, m, n;
  int  low = 0, high;
  double scale, f, g, h;
  mat  H;
  vec  ort;

  it_assert (mat_height (a) == mat_width (a), "Matrix must be square");

  it_assert (mat_height (V) == mat_width (V),
	     "Transformation matrix must be square");

  it_assert (mat_width (a) == mat_width (V),
	     "Matrices must have same dimensions");

  H = mat_clone (a);

  n = mat_height (a);

  mat_zeros (V);
  ort = vec_new_zeros (n);

  high = n - 1;

  for (m = low + 1; m <= high - 1; m++) {

    /* Scale column */

    scale = 0.0;
    for (i = m; i <= high; i++)
      scale = scale + fabs (H[i][m - 1]);

    if (scale) {
      /* Compute Householder transformation */

      h = 0.0;
      for (i = high; i >= m; i--) {
	ort[i] = H[i][m - 1] / scale;
	h += ort[i] * ort[i];
      }
      g = sqrt (h);

      if (ort[m] > 0)
	g = -g;

      h = h - ort[m] * g;
      ort[m] = ort[m] - g;

      /* 
         Apply Householder similarity transformation
         H = (I-u*u'/h)*H*(I-u*u')/h)
       */

      for (j = m; j < n; j++) {
	f = 0.0;
	for (i = high; i >= m; i--)
	  f += ort[i] * H[i][j];

	f /= h;

	for (i = m; i <= high; i++)
	  H[i][j] -= f * ort[i];
      }

      for (i = 0; i <= high; i++) {
	f = 0.0;
	for (j = high; j >= m; j--)
	  f += ort[j] * H[i][j];

	f /= h;

	for (j = m; j <= high; j++)
	  H[i][j] -= f * ort[j];

      }
      ort[m] = scale * ort[m];
      H[m][m - 1] = scale * g;

    }				/* IF scale */
  }

  /* Accumulate transformations (Algol's ortran) */

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      V[i][j] = (i == j ? 1.0 : 0.0);

  for (m = high - 1; m >= low + 1; m--) {
    if (H[m][m - 1]) {
      for (i = m + 1; i <= high; i++)
	ort[i] = H[i][m - 1];

      for (j = m; j <= high; j++) {
	g = 0.0;
	for (i = m; i <= high; i++)
	  g += ort[i] * V[i][j];

	/* Double division avoids possible underflow */
	g = (g / ort[m]) / H[m][m - 1];
	for (i = m; i <= high; i++)
	  V[i][j] += g * ort[i];
      }
    }
  }

  vec_delete (ort);

  return (H);

}

static cvec mat_h2rs (mat H, mat evec)
{

  int  nn, n, low, high;
  int  i, j, k, l, m, mmax, mmin, iter, notlast;
  double eps;
  double exshift = 0.0;
  double p = 0, q = 0, r = 0, s = 0, z = 0, t, w, x, y;
  double norm = 0.0;
  double ra, sa, vr, vi;
  double cdivr, cdivi;
  cvec eval;
  vec  e, d;

  it_assert (mat_width (H) == mat_height (H), "Matrix must be square");

  nn = mat_height (H);
  n = nn - 1;
  low = 0;
  high = nn - 1;
  eps = pow (2.0, -52.0);

  e = vec_new_zeros (nn);
  d = vec_new_zeros (nn);

  eval = cvec_new_zeros (nn);

  /* Store roots isolated by balanc and compute matrix norm */

  /* Initialize */

  for (i = 0; i < nn; i++) {
    if (i < low || i > high) {
      d[i] = H[i][i];
      e[i] = 0.0;
    }

    mmax = i - 1 > 0 ? i - 1 : 0;

    for (j = mmax; j < nn; j++)
      norm += fabs (H[i][j]);

  }

  /* Outer loop over eigenvalue index */

  iter = 0;

  while (n >= low) {

    /* Look for single small sub-diagonal element */

    l = n;

    while (l > low) {
      s = fabs (H[l - 1][l - 1]) + fabs (H[l][l]);

      if (!s)
	s = norm;

      if (fabs (H[l][l - 1]) < eps * s)
	break;

      l--;
    }

    /* Check for convergence */
    /* One root found */

    if (l == n) {
      H[n][n] += exshift;
      d[n] = H[n][n];
      e[n] = 0.0;
      n--;
      iter = 0;

      /* Two roots found */

    }
    else if (l == n - 1) {
      w = H[n][n - 1] * H[n - 1][n];
      p = (H[n - 1][n - 1] - H[n][n]) / 2.;
      q = p * p + w;
      z = sqrt (fabs (q));
      H[n][n] += exshift;
      H[n - 1][n - 1] += exshift;
      x = H[n][n];

      /* Real pair */

      if (q >= 0) {
	if (p >= 0)
	  z = p + z;
	else
	  z = p - z;

	d[n - 1] = x + z;
	d[n] = d[n - 1];
	if (z)
	  d[n] = x - w / z;

	e[n - 1] = 0.0;
	e[n] = 0.0;
	x = H[n][n - 1];
	s = fabs (x) + fabs (z);
	p = x / s;
	q = z / s;
	r = sqrt (p * p + q * q);
	p /= r;
	q /= r;

	/* Row modification */

	for (j = n - 1; j < nn; j++) {
	  z = H[n - 1][j];
	  H[n - 1][j] = q * z + p * H[n][j];
	  H[n][j] = q * H[n][j] - p * z;
	}

	/* Column modification */

	for (i = 0; i <= n; i++) {
	  z = H[i][n - 1];
	  H[i][n - 1] = q * z + p * H[i][n];
	  H[i][n] = q * H[i][n] - p * z;
	}

	/* Accumulate transformations */

	for (i = low; i <= high; i++) {
	  z = evec[i][n - 1];
	  evec[i][n - 1] = q * z + p * evec[i][n];
	  evec[i][n] = q * evec[i][n] - p * z;
	}

	/* Complex pair */

      }
      else {
	d[n - 1] = x + p;
	d[n] = x + p;
	e[n - 1] = z;
	e[n] = -z;
      }

      n -= 2;
      iter = 0;

      /* No convergence yet */

    }
    else {
      /* Form shift */

      x = H[n][n];
      y = 0.0;
      w = 0.0;
      if (l < n) {
	y = H[n - 1][n - 1];
	w = H[n][n - 1] * H[n - 1][n];
      }

      /* Wilkinson's original ad hoc shift */

      if (iter == 10) {
	exshift += x;
	for (i = low; i <= n; i++)
	  H[i][i] -= x;

	s = fabs (H[n][n - 1]) + fabs (H[n - 1][n - 2]);
	x = y = 0.75 * s;
	w = -0.4375 * s * s;
      }

      /* MATLAB's new ad hoc shift */

      if (iter == 30) {
	s = (y - x) / 2.0;
	s = s * s + w;
	if (s > 0) {
	  s = sqrt (s);
	  if (y < x)
	    s = -s;
	  s = x - w / ((y - x) / 2. + s);
	  for (i = low; i <= n; i++)
	    H[i][i] -= s;
	  exshift += s;
	  x = y = w = 0.964;
	}
      }

      iter++;			/* (Could check iteration count here.) */

      /* Look for two consecutive small sub-diagonal elements */

      m = n - 2;
      while (m >= l) {
	z = H[m][m];
	r = x - z;
	s = y - z;
	p = (r * s - w) / H[m + 1][m] + H[m][m + 1];
	q = H[m + 1][m + 1] - z - r - s;
	r = H[m + 2][m + 1];
	s = fabs (p) + fabs (q) + fabs (r);
	p /= s;
	q /= s;
	r /= s;
	if (m == l)
	  break;

	if (fabs (H[m][m - 1]) * (fabs (q) + fabs (r)) <
	    eps * (fabs (p) *
		   (fabs (H[m - 1][m - 1]) + fabs (z) +
		    fabs (H[m + 1][m + 1]))))
	  break;

	m--;
      }

      for (i = m + 2; i <= n; i++) {
	H[i][i - 2] = 0.;
	if (i > m + 2)
	  H[i][i - 3] = 0.;
      }

      /* Double QR step involving rows l:n and columns m:n */

      for (k = m; k <= n - 1; k++) {
	notlast = (k != n - 1);
	if (k != m) {
	  p = H[k][k - 1];
	  q = H[k + 1][k - 1];
	  r = (notlast ? H[k + 2][k - 1] : 0.);
	  x = fabs (p) + fabs (q) + fabs (r);
	  if (x) {
	    p /= x;
	    q /= x;
	    r /= x;
	  }
	}

	if (!x)
	  break;

	s = sqrt (p * p + q * q + r * r);
	if (p < 0)
	  s = -s;

	if (s) {
	  if (k != m)
	    H[k][k - 1] = -s * x;
	  else if (l != m)
	    H[k][k - 1] = -H[k][k - 1];

	  p += s;
	  x = p / s;
	  y = q / s;
	  z = r / s;
	  q /= p;
	  r /= p;

	  /* Row modification */

	  for (j = k; j < nn; j++) {
	    p = H[k][j] + q * H[k + 1][j];
	    if (notlast) {
	      p = p + r * H[k + 2][j];
	      H[k + 2][j] = H[k + 2][j] - p * z;
	    }
	    H[k][j] = H[k][j] - p * x;
	    H[k + 1][j] = H[k + 1][j] - p * y;
	  }

	  /* Column modification */

	  mmin = n < k + 3 ? n : k + 3;

	  for (i = 0; i <= mmin; i++) {
	    p = x * H[i][k] + y * H[i][k + 1];
	    if (notlast) {
	      p += z * H[i][k + 2];
	      H[i][k + 2] -= p * r;
	    }
	    H[i][k] -= p;
	    H[i][k + 1] -= p * q;
	  }

	  /* Accumulate transformations */

	  for (i = low; i <= high; i++) {
	    p = x * evec[i][k] + y * evec[i][k + 1];
	    if (notlast) {
	      p += z * evec[i][k + 2];
	      evec[i][k + 2] -= p * r;
	    }
	    evec[i][k] -= p;
	    evec[i][k + 1] -= p * q;
	  }

	}			/* (s != 0) */
      }				/* k loop */
    }				/* check convergence */
  }				/* while (n >= low) */

  /* Backsubstitute to find vectors of upper triangular form */

  /* We return null vector if there is a problem */
  if (!norm)
    return (eval);

  for (n = nn - 1; n >= 0; n--) {
    p = d[n];
    q = e[n];

    /* Real vector */

    if (!q) {
      l = n;
      H[n][n] = 1.;
      for (i = n - 1; i >= 0; i--) {
	w = H[i][i] - p;
	r = 0.;
	for (j = l; j <= n; j++)
	  r += H[i][j] * H[j][n];

	if (e[i] < 0.0) {
	  z = w;
	  s = r;
	}
	else {
	  l = i;
	  if (!e[i]) {
	    if (w != 0.0)
	      H[i][n] = -r / w;
	    else
	      H[i][n] = -r / (eps * norm);

	    /* Solve real equations */

	  }
	  else {
	    x = H[i][i + 1];
	    y = H[i + 1][i];
	    q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
	    t = (x * s - z * r) / q;
	    H[i][n] = t;
	    if (fabs (x) > fabs (z))
	      H[i + 1][n] = (-r - w * t) / x;
	    else
	      H[i + 1][n] = (-s - y * t) / z;
	  }

	  /* Overflow control */

	  t = fabs (H[i][n]);
	  if ((eps * t) * t > 1)
	    for (j = i; j <= n; j++)
	      H[j][n] /= t;
	}
      }

      /* Complex vector */

    }
    else if (q < 0) {
      l = n - 1;

      /* Last vector component imaginary so matrix is triangular */

      if (fabs (H[n][n - 1]) > fabs (H[n - 1][n])) {
	H[n - 1][n - 1] = q / H[n][n - 1];
	H[n - 1][n] = -(H[n][n] - p) / H[n][n - 1];
      }
      else {
	ccdiv (0.0, -H[n - 1][n], H[n - 1][n - 1] - p, q, &cdivr, &cdivi);
	H[n - 1][n - 1] = cdivr;
	H[n - 1][n] = cdivi;
      }
      H[n][n - 1] = 0.;
      H[n][n] = 1.;
      for (i = n - 2; i >= 0; i--) {
	ra = 0.0;
	sa = 0.0;
	for (j = l; j <= n; j++) {
	  ra += H[i][j] * H[j][n - 1];
	  sa += H[i][j] * H[j][n];
	}
	w = H[i][i] - p;

	if (e[i] < 0.0) {
	  z = w;
	  r = ra;
	  s = sa;
	}
	else {
	  l = i;
	  if (e[i]) {
	    ccdiv (-ra, -sa, w, q, &cdivr, &cdivi);
	    H[i][n - 1] = cdivr;
	    H[i][n] = cdivi;
	  }
	  else {

	    /* Solve complex equations */

	    x = H[i][i + 1];
	    y = H[i + 1][i];
	    vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
	    vi = (d[i] - p) * 2.0 * q;
	    if (!vr && !vi)
	      vr =
		eps * norm * (fabs (w) + fabs (q) + fabs (x) + fabs (y) +
			      fabs (z));

	    ccdiv (x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi,
		   &cdivr, &cdivi);
	    H[i][n - 1] = cdivr;
	    H[i][n] = cdivi;
	    if (fabs (x) > (fabs (z) + fabs (q))) {
	      H[i + 1][n - 1] = (-ra - w * H[i][n - 1] + q * H[i][n]) / x;
	      H[i + 1][n] = (-sa - w * H[i][n] - q * H[i][n - 1]) / x;
	    }
	    else {
	      ccdiv (-r - y * H[i][n - 1], -s - y * H[i][n], z, q, &cdivr,
		     &cdivi);
	      H[i + 1][n - 1] = cdivr;
	      H[i + 1][n] = cdivi;
	    }
	  }

	  /* Overflow control */

	  t =
	    fabs (H[i][n - 1]) >
	    fabs (H[i][n]) ? fabs (H[i][n - 1]) : fabs (H[i][n]);

	  if ((eps * t) * t > 1.)
	    for (j = i; j <= n; j++) {
	      H[j][n - 1] /= t;
	      H[j][n] /= t;
	    }
	}
      }
    }
  }

  /* Vectors of isolated roots */

  for (i = 0; i < nn; i++)
    if (i < low || i > high)
      for (j = i; j < nn; j++)
	evec[i][j] = H[i][j];

  /* Back transformation to get eigenvectors of original matrix */

  for (j = nn - 1; j >= low; j--) {
    for (i = low; i <= high; i++) {
      z = 0.;

      mmin = j < high ? j : high;

      for (k = low; k <= mmin; k++)
	z += evec[i][k] * H[k][j];
      evec[i][j] = z;
    }
  }

  for (i = 0; i < vec_length (d); i++) {
    creal (eval[i]) = d[i];
    cimag (eval[i]) = e[i];
  }

  vec_delete (e);
  vec_delete (d);

  return (eval);

}

/*
  
  Eigenvalues/eigenvectors decomposition of a real symmetric 
  matrix. The output is real in every cases (Hermitian). 

 */

vec mat_eig_sym (mat a, mat evec)
{

  int  i, j;
  vec  e, d;

  it_assert (mat_width (a) == mat_height (a), "Matrix must be square");
  it_assert (evec != NULL,
	     "Memory for matrix of eigenvectors must be reserved");
  it_assert (mat_width (evec) == mat_height (evec),
	     "Matrix of eigenvectors must be square");
  it_assert (mat_width (evec) == mat_width (a),
	     "Matrices dimensions must agree");

  e = vec_new_zeros (mat_width (a));
  d = vec_new_zeros (mat_width (a));

  for (i = 0; i < mat_height (a); i++)
    for (j = 0; j < mat_width (a); j++)
      evec[i][j] = a[i][j];

  mat_tridiag (e, d, evec);
  mat_tridiag_ql (e, d, evec);

  vec_delete (e);

  return (d);
}

/*

  Eigenvalues/eigenvectors decomposition of a real *NON* symmetric 
  matrix. Well, it works with symmetric matrices too, but it is 
  slower. The output is generally complex. 

 */

cvec mat_eig (mat a, cmat evec)
{

  int  i, j, n;
  cvec eval;
  mat  h, V;
  double z;

  it_assert (mat_width (a) == mat_height (a), "Matrix must be square");
  it_assert (evec != NULL,
	     "Memory for matrix of eigenvectors must be reserved");
  it_assert (cmat_width (evec) == cmat_height (evec),
	     "Matrix of eigenvectors must be square");
  it_assert (cmat_width (evec) == mat_width (a),
	     "Matrices dimensions must agree");

  V = mat_clone (a);

  n = mat_width (a);

  h = mat_hessenberg (a, V);
  eval = mat_h2rs (h, V);
  mat_delete (h);

  /* Outputs (possibly) complex normalized eigenvectors from real Schur form */

  for (j = 0; j < n; j++) {
    if (cimag (eval[j])) {

      z = 0.;

      for (i = 0; i < n; i++) {
	creal (evec[i][j]) = creal (evec[i][j + 1]) = V[i][j];
	cimag (evec[i][j]) = V[i][j + 1];
	cimag (evec[i][j + 1]) = -V[i][j + 1];
	z +=
	  creal (evec[i][j]) * creal (evec[i][j]) +
	  cimag (evec[i][j]) * cimag (evec[i][j]);
      }

      z = 1. / sqrt (z);

      for (i = 0; i < n; i++) {
	creal (evec[i][j]) *= z;
	cimag (evec[i][j]) *= z;
	creal (evec[i][j + 1]) *= z;
	cimag (evec[i][j + 1]) *= z;
      }

      j++;			/* Skip conj */
    }
    else {

      z = 0.;

      for (i = 0; i < n; i++) {
	cimag (evec[i][j]) = 0.;
	creal (evec[i][j]) = V[i][j];
	z += V[i][j] * V[i][j];
      }

      z = 1. / sqrt (z);

      /* Mimic Octave behaviour */
      if (V[0][j] < 0.)
	z = -z;

      for (i = 0; i < n; i++)
	creal (evec[i][j]) *= z;

    }

  }

  mat_delete (V);

  return (eval);

}


/*
  Like the Hessenberg function, it could be declared static. 
 */
mat mat_real_schur (mat a)
{
  mat  h, e;
  cvec eval;

  it_assert (mat_height (a) == mat_width (a), "Matrix must be square");

  e = mat_new_zeros (mat_height (a), mat_width (a));

  h = mat_hessenberg (a, e);
  eval = mat_h2rs (h, e);

  cvec_delete (eval);
  mat_delete (e);

  return (h);
}

double  mat_lu( mat a, ivec piv ) 
{

  idx_t i = 0; 
  idx_t j = 0; 
  idx_t k = 0; 
  idx_t pi = 0; 
  idx_t xp, yp; 
  idx_t t = 0; 

  int pivsign = 1; 

  size_t n, m, p; 

  double max = 0.0; 
  double eps = IT_EPSILON*IT_EPSILON; 

  it_assert( a, "Please use an existing matrix" ); 
  it_assert( pc, "Please use an existing pc vector" ); 
  it_assert( pr, "Please use an existing pr vector" ); 
  it_assert( ivec_length( piv ) == mat_height( a ), "Vector piv must have matrix height size" ); 

  m = mat_height( a );
  n = mat_width( a );

  for ( i= 0; i< ivec_length( piv ); i++ ) 
    piv[i] = i; 

  p = (n<=m)?n:m; 

  /* Iterate on columns */
  for ( k= 0; k< p; k++ ) 
    {
      /* Find pivot */ 
      max = a[k][k]; 
      xp = yp = k; 

      for ( i= k; i< m; i++ ) 
	if ( fabs(max) < fabs( a[i][k] ) ) 
	  {
	    max = a[i][k]; 
	    yp = i; 
	  }
            
      if ( fabs(max) > eps ) /* We have found a non-zero pivot */
	{
	  /* Swap rows if necessary */ 
	  if ( yp!=k ) 
	    {
	      mat_swap_rows( a, yp, k );
	      
	      /* Update pivot vector */ 
	      t = piv[yp]; 
	      piv[yp] = piv[k]; 
	      piv[k] = t; 

	      pivsign = -pivsign; 
	      
	    }      

	  /* Compute U coeffs */
	  /* Iterate on rows for U */
	  for ( j= k; j< n; j++ )
	    for ( pi= 0; pi< k; pi++ )
	      a[k][j] -= a[k][pi]*a[pi][j];
	  
	  /* Compute L coeffs */
	  for ( i= k+1; i< m; i++ ) 
	    {
	      for ( pi= 0; pi< k; pi++ ) 
		a[i][k] -= a[i][pi]*a[pi][k]; 	  
	      a[i][k] /= a[k][k]; 
	    } /* L */	  

	}
      else /* matrix is singular */
	{
	  it_warning( "Matrix is singular" );
	  return( 0. ); 
	}

    }

  /* Computing determinant : priceless... */
  max = (double)pivsign; 

  for ( i= 0; i< n; i++ ) 
    max *= a[i][i];

  return( max );

}

/* Compute determinant of a matrix. Returns zero if the matrix is singular. */

double mat_det (mat a)
{
  double det = 0.;
  ivec piv; 
  mat  lu = mat_clone (a);

  piv = ivec_new_zeros (mat_height (a));

  /* Use total pivoting strategy, useless? */
  det = mat_lu (lu, piv ); 

  ivec_delete (piv);
  mat_delete (lu);

  return (det);
}


/* Use LU for solving A x = b. */
vec mat_solve_vec (mat A, vec b)
{
  size_t m, n; 
  idx_t i, k; 
  ivec piv; 
  mat LU; 
  vec X; 

  m  = mat_height( A ); 
  n  = mat_width( A ); 

  it_assert( m==vec_length(b), "Matrices size mismatch" );

  piv = ivec_new( m ); 
  LU = mat_clone( A ); 

  if ( !mat_lu( LU, piv ) ) 
    {
      it_warning( "Matrix is singular, aborting" ); 
      ivec_delete( piv ); 
      mat_delete( LU ); 
      return( NULL );
    }
  
  X = vec_new( ivec_length(piv) );

  for ( i= 0; i< ivec_length(piv); i++ ) 
    X[i] = b[piv[i]]; 

  for ( k= 0; k< n; k++ ) 
    for ( i= k+1; i< n; i++ ) 
      X[i] -= X[k]*LU[i][k]; 

  for ( k= n-1; k >= 0; k-- ) 
    {
      X[k] /= LU[k][k]; 

      for ( i= 0; i< k; i++ )
	X[i] -= X[k]*LU[i][k];
    }

  ivec_delete( piv ); 
  mat_delete( LU );

  return( X ); 
}


/*  Use LU for solving several equations A x_i = b_i.
    The various b_i are to be stored column-wise in B.
    Return column-wise solutions.                          */

mat mat_solve_mat( mat A, mat B ) 
{
  size_t m, n, nx; 
  ivec piv; 
  idx_t i, j, k; 
  mat LU, X; 

  m  = mat_height( A ); 
  n  = mat_width( A ); 
  nx = mat_width( B ); 

  it_assert( m==mat_height(B), "Matrices size mismatch" );

  piv = ivec_new( m ); 
  LU = mat_clone( A ); 

  if ( !mat_lu( LU, piv ) ) 
    {
      it_warning( "Matrix is singular, aborting" ); 
      ivec_delete( piv ); 
      mat_delete( LU ); 
      return( NULL );
    }
  
  X = mat_new( ivec_length(piv), nx ); 

  for ( i= 0; i< ivec_length(piv); i++ ) 
    for ( j= 0; j< nx; j++ ) 
      X[i][j] = B[piv[i]][j]; 

  for ( k= 0; k< n; k++ ) 
    for ( i= k+1; i< n; i++ ) 
      for ( j= 0; j< nx; j++ ) 
	X[i][j] -= X[k][j]*LU[i][k]; 

  for ( k= n-1; k >= 0; k-- ) 
    {
      for ( j= 0; j< nx; j++ ) 
	X[k][j] /= LU[k][k]; 

      for ( i= 0; i< k; i++ ) 
	for ( j= 0; j< nx; j++ ) 
	  X[i][j] -= X[k][j]*LU[i][k];
    }

  ivec_delete( piv );
  mat_delete( LU );

  return( X ); 
}

static mat mat_inv_direct (mat A)
{
  mat  I, Id;

  it_assert (mat_width (A) == mat_height (A), "Matrix must be square");

  Id = mat_new_zeros (mat_height (A), mat_width (A));
  mat_eye (Id);

  I = mat_solve_mat (A, Id);

  mat_delete (Id);
  return (I);
}

/* 
   QR decomposition
*/
vec mat_qr( mat A ) 
{

  size_t m, n; 
  idx_t i, j, k; 
  vec Rdiag; 
  double nrm, s; 

  m = mat_height( A );
  n = mat_width( A );

  Rdiag = vec_new( n ); 

  for ( k= 0; k< n; k++ )
    {
      nrm = 0.; 
      for ( i= k; i< m; i++ ) 
	nrm = hypot( nrm, A[i][k] ); 

      if ( nrm ) 
	{
	  if ( A[k][k] < 0 ) 
	    nrm = -nrm; 

	  for ( i= k; i< m; i++ ) 
	    A[i][k] /= nrm; 

	  A[k][k] += 1.; 

	  for ( j= k+1; j< n; j++ ) 
	    {
	      s = 0.; 
	      for ( i= k; i< m; i++ ) 
		s+= A[i][k]*A[i][j]; 

	      s = -s / A[k][k]; 

	      for ( i= k; i< m; i++ ) 
		A[i][j] += s*A[i][k];
	    }

	}

      Rdiag[k] = -nrm; 
    }

  return Rdiag; 

}

mat mat_ls( mat A, mat B ) 
{

  size_t m, n, nx; 
  idx_t i, j, k; 
  double s; 
  mat QR, X; 
  vec Rdiag; 

  m =  mat_height( A );
  n =  mat_width( A );
  nx = mat_width( B );

  it_assert( m==mat_height(B), "Matrices size mismatch" ); 

  QR = mat_clone( A ); 

  Rdiag = mat_qr( QR ); 

  for ( j= 0; j< n; j++ ) 
    if ( !Rdiag[j] ) 
      {
	it_warning( "Matrix is rank deficient, aborting" ); 
	vec_delete( Rdiag ); 
	mat_delete( QR ); 
	return NULL; 
      }

  X = mat_clone( B ); 

  for ( k= 0; k< n; k++ ) 
    for ( j= 0; j< nx; j++ ) 
      {
	s = 0.; 
	for ( i= k; i< m; i++ ) 
	  s+= QR[i][k]*X[i][j]; 
	
	s = -s / QR[k][k]; 
	
	for ( i= k; i< m; i++ ) 
	  X[i][j] += s*QR[i][k];
	
      }
  
  for ( k= n-1; k >= 0; k-- ) 
    {
      for ( j= 0; j< nx; j++ ) 
	X[k][j] /= Rdiag[k]; 

      for ( i= 0; i< k; i++ ) 
	for ( j= 0; j< nx; j++ ) 
	  X[i][j] -= X[k][j]*QR[i][k]; 
    }

  vec_delete( Rdiag );
  mat_delete( QR );

  return X; 

}

/*----------------------------------------------------------------------------*/
/* Orthonormalisation in-place function.                                      */
/* Implement Gram-Schmidt procedure                                           */
void mat_gs (mat A) 
{
  int  i, j, k;
  vec  v = vec_new (mat_height (A)), proj;
  double innerp;		/* for the inner product */

  for (i = 1; i < mat_width (A); i++) {
    vec_zeros (v);

    for (j = 0; j < i; j++) {
      innerp = 0.;
      for (k = 0; k < mat_height (A); k++)
	innerp += A[k][i] * A[k][j];

      proj = mat_get_col (A, j);
      vec_mul_by (proj, innerp / pow( vec_norm( proj, 2.), 2.));

      vec_add (v, proj);
      vec_delete (proj);
    }

    /* deduce v from column vector i */
    for (j = 0; j < mat_height (A); j++)
      A[j][i] -= v[j];
  }

  /* Proceed to normalization when possible
     Priority is given to orthogonalisation
   */
  if ( mat_width( A ) < mat_height( A ) )
    mat_cols_normalize( A, 2. );

  vec_delete (v);
}

void mat_cholesky( mat a ) 
{

  idx_t i = 0; 
  idx_t j = 0; 
  idx_t k = 0; 
  size_t n = 0; 

  it_assert( a, "Please use an existing matrix" ); 
  it_assert( mat_width( a )==mat_height( a ), "Matrix must be square" );

  n = mat_width( a );

  for ( k= 0; k< n; k++ ) 
    {
      a[k][k] = sqrt( a[k][k] );
      for ( i= k+1; i< n; i++ )
	a[i][k] /= a[k][k];
      for ( j= k+1; j< n; j++ ) 
	for ( i= j; i< n; i++ ) 
	  a[i][j] -= a[i][k]*a[j][k];
    }

  for ( j= 0; j< n; j++ ) 
    for ( i= 0; i< j; i++ ) 
      a[i][j] = 0.0;

  return;
}

/*

  Inversion of a matrix using LU. 

 */

mat mat_inv (mat m)
{

  mat  inv;
  idx_t i, j;

  inv = mat_inv_direct (m);

  if (inv) {
    for (i = 0; i < mat_height (m); i++)
      for (j = 0; j < mat_width (m); j++)
	m[i][j] = inv[i][j];

    mat_delete (inv);

    return (m);

  }
  else
    return (NULL);

}


/*----------------------------------------------------------------------------*/
mat mat_new_inv (mat m)
{
  return mat_inv_direct (m);
}


/*----------------------------------------------------------------------------*/
cmat cmat_inv (cmat m)
{
  idx_t i, j, k;
  cplx coef;
  cvec f = cvec_new (cmat_height (m));

  it_assert (cmat_height (m) == cmat_width (m), "Matrix must be square");

  /* Gauss method */
  for (i = 0; i < cmat_height (m); i++) {
    if (cnorm (m[i][i]) < IT_EPSILON * 1e-5) {	/* Singular matrix ? */
      cmat_delete (m);
      cvec_delete (f);
      it_warning
	("Matrix is singular or you are unlucky. Try another method for inversion.\n");
      m = cmat_new_void ();
      return m;
    }

    f[i] = m[i][i];
    cvec_div_by (m[i], f[i]);

    for (j = 0; j < cmat_height (m); j++) {
      if (i != j) {
	coef = m[j][i];
	m[j][i] = cplx_0;

	for (k = 0; k < cmat_width (m); k++)
	  m[j][k] = csub (m[j][k], cmul (coef, m[i][k]));
      }
    }
  }

  /* Normalization  */
  for (i = 0; i < cmat_height (m); i++)
    cvec_div_by (m[i], f[i]);

  cvec_delete (f);
  return m;
}


/*----------------------------------------------------------------------------*/
cmat cmat_new_inv (cmat m)
{
  cmat r = cmat_clone (m);
  return cmat_inv (r);
}
