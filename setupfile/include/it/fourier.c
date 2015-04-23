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
  Discrete Fourier transform
  Copyright (C) 2005 Vivien Chappelier.

  The current algorithm is O(nlogn) for any sizes. Power-of-two sizes use the
  FFT, while other sizes use DFT computed as a special case of z-transform.
  The z-transform itself done by convolution, which is done by FFT.

  This is certainly not the fastest FFT in the world, possible optimizations
  include using higher-radix FFTs and special cases for FFT of real data.

  Many of the tricks used here are well-documented in the following books,
  please refer to them for details on the algorithms:
  - "DSP Algorithms for Programmers"
    by Jorg Arndt
    [http://www.jjj.de/fxt/fxtbook.pdf] 
  - "Numerical Recipies in C"
    by W.H. Press, S.A. Teukolsky, W.T. Wetterling, B.P. Flannery
    [http://www.library.cornell.edu/nr/bookcpdf.html]

*/

#include <it/types.h>
#include <it/fourier.h>
#include <it/io.h>
#include <it/math.h>

/* buffer size for which memory access latency is not limiting performance */
#define CACHE_SIZE 128*1024	/* sane default, to be automatically tuned later */
#define DFT_LAST   16		/* DFT of smaller size are done directly */

static unsigned int cache_size = CACHE_SIZE;	/* tunable */

/* FFT context, read-only once initialized */
static int fft_init_done = 0;
static int fft_log2size;
static cvec dft_roots[DFT_LAST];
static cvec fft_roots;
static ivec fft_perm;

/* return i such that 1 << i = n if n is a non-zero power of 2 */
static inline int intlog2 (int n)
{
  int  i;

  for (i = 0; n; i++, n >>= 1);

  return (i - 1);
}

/* initialize the FFT context */
void fft_init (void)
{
  int  size;
  int  i, k, m;

  if (fft_init_done)
    return;

  /* FFT data and precomputed roots must both fit in the cache */
  size = cache_size / sizeof (cplx) / 2;
  fft_log2size = intlog2 (size);
  fft_roots = cvec_new_unit_roots (size);
  cvec_push (fft_roots, cplx_1);	/* XXX: should not be needed */

  /* precompute the bit inversion permutation */
  fft_perm = ivec_new (size);
  fft_perm[0] = 0;
  fft_perm[size - 1] = size - 1;
  for (i = 1, k = 0; i < size - 1; i++) {
    m = size;
    do {
      m >>= 1;
      k ^= m;
    }
    while ((k & m) == 0);
    fft_perm[i] = k;
  }

  /* precompute the first few roots of unity for direct DFT */
  for (i = 2; i < DFT_LAST; i++)
    dft_roots[i] = cvec_new_unit_roots (i);

  fft_init_done = 1;
}

/* permute the elements of an array of complex numbers by inverting the bit */
/* representation of the indexes */
static inline void cplx_array_bitrev_permute (cplx * v, int n)
{
  int  i, m, r;

  if (n <= 2)
    return;

  r = 0;
  for (i = 1; i < n - 1; i++) {
    m = n;
    do {
      m >>= 1;
      r ^= m;
    }
    while ((r & m) == 0);
    if (r > i)
      cplx_swap (v[i], v[r]);
  }
}

/* The Core FFT algorithm.
 * It is quite restrictive as the input must have a length which is a power of
 * two and less than half the cache size divided by the size of the complex
 * type. All other DFT algorithm are derived from this one.
 * No checks are made at this point, they are the responsibility of the
 * wrappers. The index permutation is done outside this function as well.
 */

/* forward FFT by decimation in frequency */
static void __fft_dif (cplx * v, int log2n)
{
  cvec roots = fft_roots;
  int  log2f = fft_log2size;
  idx_t n, m, mh, s;
  idx_t k, l;
  idx_t log2m;
  cplx a, b;
  cplx *r;

  n = 1 << log2n;

  for (log2m = log2n; log2m > 0; log2m--) {
    m = 1 << log2m;
    mh = m >> 1;
    s = 1 << (log2f - log2m);

    r = &roots[1 << log2f];
    for (k = 0; k < mh; k++) {
      for (l = 0; l < n; l += m) {
	a = v[l + k];
	b = v[l + k + mh];
	v[l + k] = cadd (a, b);
	v[l + k + mh] = cmul (csub (a, b), *r);
      }
      r -= s;
    }
  }
}

/* backward FFT by decimation in time */
static void __ifft_dit (cplx * v, int log2n)
{
  cvec roots = fft_roots;
  int  log2f = fft_log2size;
  idx_t n, m, mh, s;
  idx_t k, l;
  idx_t log2m;
  cplx a, b;
  cplx *r;

  n = 1 << log2n;

  for (log2m = 1; log2m <= log2n; log2m++) {
    m = 1 << log2m;
    mh = m >> 1;
    s = 1 << (log2f - log2m);

    r = &roots[0];
    for (k = 0; k < mh; k++) {
      for (l = 0; l < n; l += m) {
	a = v[l + k];
	b = cmul (v[l + k + mh], *r);
	v[l + k] = cadd (a, b);
	v[l + k + mh] = csub (a, b);
      }
      r += s;
    }
  }
}

/* The Recursive FFT algorithm.
 * If the input vector fits in the cache, the Core FFT algorithm is called.
 * Otherwise, the MFA algorithm is applied on the input vector, considered as
 * a matrix with __fft_size rows. Each column is transformed efficiently
 * by the core FFT algorithm, and multiplied by the appropriate roots of unity.
 * The rows are then transformed by the Recursive FFT algorithm. Putting all
 * the coefficient back in order is done elsewhere as some FFT-based algorithms
 * don't need it.
 */

static void _fft (cplx * v, int log2n)
{
  int  log2f = fft_log2size;
  int  f = 1 << log2f;
  int  n = 1 << log2n;
  cvec tmp;
  int  i, j, k;
  int  r, log2r;
  double c;
  cplx *p;
  cplx e;

  if (log2n <= log2f) {
    /* fits in the cache */
    __fft_dif (v, log2n);
    return;
  }

  /* use the MFA algorithm considering an size x r matrix */
  log2r = log2n - log2f;
  r = (1 << log2r);

  /* transform the columns */
  tmp = cvec_new (f);

  for (j = 0; j < r; j++) {

    /* fetch a column */
    for (i = 0; i < f; i++)
      tmp[i] = v[j + r * i];

    /* transform it */
    __fft_dif (tmp, log2f);

    /* store, permute and multiply terms by e^(- 2*cplx_I*M_PI*i*j/n) */
    p = &v[j];
    *p = tmp[0];
    p += r;

    c = -2 * M_PI * j / n;
    for (i = 1; i < f - 1; i++) {

      k = fft_perm[i];

      e.r = cos (c * i);
      e.i = sin (c * i);

      *p = cmul (tmp[k], e);
      p += r;
    }

    e.r = cos (c * (f - 1));
    e.i = sin (c * (f - 1));
    *p = cmul (tmp[f - 1], e);
  }

  cvec_delete (tmp);

  /* transform the rows */
  p = v;
  for (j = 0; j < f; j++) {
    _fft (p, log2r);
    p += r;
  }
}

static void _ifft (cplx * v, int log2n)
{
  int  log2f = fft_log2size;
  int  f = 1 << log2f;
  int  n = 1 << log2n;
  cvec tmp;
  int  i, j, k;
  int  r, log2r;
  double c;
  cplx *p;
  cplx e;

  if (log2n <= log2f) {
    /* fits in the cache */
    __ifft_dit (v, log2n);
    return;
  }

  /* use the MFA algorithm considering an size x r matrix */
  log2r = log2n - log2f;
  r = (1 << log2r);

  /* transform the rows */
  p = v;
  for (j = 0; j < f; j++) {
    _ifft (p, log2r);
    p += r;
  }

  /* transform the columns */
  tmp = cvec_new (f);

  for (j = 0; j < r; j++) {

    /* fetch, permute and multiply terms by e^(2*cplx_I*M_PI*i*j/n) */
    p = &v[j];
    tmp[0] = *p;
    p += r;

    c = 2 * M_PI * j / n;
    for (i = 1; i < f - 1; i++) {

      k = fft_perm[i];

      e.r = cos (c * i);
      e.i = sin (c * i);

      tmp[k] = cmul (*p, e);
      p += r;
    }

    e.r = cos (c * (f - 1));
    e.i = sin (c * (f - 1));
    tmp[f - 1] = cmul (*p, e);

    /* transform it */
    __ifft_dit (tmp, log2f);

    /* store a column */
    for (i = 0; i < f; i++)
      v[j + r * i] = tmp[i];
  }

  cvec_delete (tmp);
}

/* Put back the FFT coefficient in correct order */
/* TODO: compute the permutation explicitly to get rid of recursivity */
/* the transform on an index [ a b c d ] of f+f+f+<f bits is          */
/*   f   f   f   <f          <f  f   f   f                            */
/* [ a   b   c   d  ]  ->  [ d'  c   b   a  ]                         */
/* with d' the bit reversed binary representation of d                */

static void _fft_demangle (cplx * v, int log2n)
{
  int  log2f = fft_log2size;
  int  f = 1 << log2f;
  int  n = 1 << log2n;
  int  i, j;
  int  r, log2r;
  cplx *p;

  if (log2n <= log2f) {
    cplx_array_bitrev_permute (v, n);
    return;
  }

  log2r = log2n - log2f;
  r = (1 << log2r);

  /* transpose the rows */
  if (log2r > log2f) {
    p = v;
    for (j = 0; j < f; j++) {
      _fft_demangle (p, log2r);
      cplx_array_bitrev_permute (p, r);
      p += r;
    }
  }

  cplx_array_bitrev_permute (v, n);

  p = v;
  for (i = 0; i < r; i++) {
    cplx_array_bitrev_permute (p, f);
    p += f;
  }
}

static void _ifft_mangle (cplx * v, int log2n)
{
  int  log2f = fft_log2size;
  int  f = 1 << log2f;
  int  n = 1 << log2n;
  int  i, j;
  int  r, log2r;
  cplx *p;

  if (log2n <= log2f) {
    cplx_array_bitrev_permute (v, n);
    return;
  }

  log2r = log2n - log2f;
  r = (1 << log2r);

  /* transpose */
  p = v;
  for (i = 0; i < r; i++) {
    cplx_array_bitrev_permute (p, f);
    p += f;
  }

  cplx_array_bitrev_permute (v, n);

  /* transform the rows */
  if (log2r > log2f) {
    p = v;
    for (j = 0; j < f; j++) {
      cplx_array_bitrev_permute (p, r);
      _ifft_mangle (p, log2r);
      p += r;
    }
  }
}

/* Fast Fourier Transform */
cvec it_fft (cvec v)
{
  int  n, log2n;
  cvec cv;

  assert (cvec_length (v));

  n = cvec_length (v);
  log2n = intlog2 (n);
  assert ((1 << log2n) == n);
  fft_init ();
  cv = cvec_clone (v);
  _fft (cv, log2n);
  _fft_demangle (cv, log2n);
  return (cv);
}

/* Inverse Fast Fourier Transform */
cvec it_ifft (cvec v)
{
  int  n, log2n;
  cvec cv;

  assert (cvec_length (v));

  n = cvec_length (v);
  log2n = intlog2 (n);
  assert ((1 << log2n) == n);
  fft_init ();
  cv = cvec_clone (v);
  _ifft_mangle (cv, log2n);
  _ifft (cv, log2n);
  cvec_div_by_real (cv, n);
  return (cv);
}

/* convolution by FFT */
cvec cvec_fft_conv (cvec a, cvec b)
{
  cvec ta, tb;
  int  i, l, n;
  int  log2n;

  /* copy and zero pad input to 2^n */
  ta = cvec_clone (a);
  tb = cvec_clone (b);
  l = cvec_length (a);
  if (cvec_length (b) > l)
    l = cvec_length (b);
  log2n = intlog2 (l);
  if ((1 << log2n) != l)
    log2n++;
  log2n++;
  n = 1 << log2n;

  cvec_set_length (ta, n);
  cvec_set_length (tb, n);

  for (i = cvec_length (a); i < n; i++)
    ta[i] = cplx_0;
  for (i = cvec_length (b); i < n; i++)
    tb[i] = cplx_0;

  /* convolve by product in FFT domain */
  fft_init ();
  _fft (ta, log2n);
  _fft (tb, log2n);
  cvec_mul (ta, tb);
  cvec_delete (tb);
  _ifft (ta, log2n);
  cvec_div_by_real (ta, n);
  cvec_set_length (ta, cvec_length (a) + cvec_length (b) - 1);
  return ta;
}

/* correlation by FFT */
cvec cvec_fft_corr (cvec a, cvec b)
{
  cvec ta, tb;
  int  i, l, n;
  int  log2n;

  /* copy and zero pad input to 2^n */
  ta = cvec_clone (a);
  tb = cvec_clone (b);
  l = cvec_length (a);
  if (cvec_length (b) > l)
    l = cvec_length (b);
  log2n = intlog2 (l);
  if ((1 << log2n) != l)
    log2n++;
  log2n++;
  n = 1 << log2n;

  cvec_set_length (ta, n);
  cvec_set_length (tb, n);

  for (i = cvec_length (a); i < n; i++)
    ta[i] = cplx_0;
  for (i = cvec_length (b); i < n; i++)
    tb[i] = cplx_0;

  /* correlate by conjugate product in FFT domain */
  fft_init ();
  _fft (ta, log2n);
  _fft (tb, log2n);
  cvec_conj_mul (ta, tb);
  cvec_delete (tb);
  _ifft (ta, log2n);
  cvec_div_by_real (ta, n);
  cvec_set_length (ta, l);
  return ta;
}

/* autocorrelation by FFT */
cvec cvec_fft_autocorr (cvec a)
{
  cvec ta;
  int  i, l, n;
  int  log2n;

  /* copy and zero pad input to 2^n */
  ta = cvec_clone (a);
  l = cvec_length (a);
  log2n = intlog2 (l);
  if ((1 << log2n) != l)
    log2n++;
  log2n++;
  n = 1 << log2n;

  cvec_set_length (ta, n);

  for (i = cvec_length (a); i < n; i++)
    ta[i] = cplx_0;

  /* correlate by conjugate product in FFT domain */
  fft_init ();
  _fft (ta, log2n);
  cvec_abssqr (ta);
  _ifft (ta, log2n);
  cvec_div_by_real (ta, n);
  cvec_set_length (ta, l);
  return ta;
}

/* z-transform through FFT */
cvec it_fzt (cvec v, cplx z)
{
  int  i, l, n;
  int  log2n;
  cplx c, d;
  cvec av, bv, cv;
  double p;
  double rho;
  double theta;
  double crho;
  double s;

  assert (cvec_length (v));
  assert (cnorm (z) != 0);

  l = cvec_length (v);

  log2n = intlog2 (l);
  if ((1 << log2n) != l)
    log2n++;
  log2n++;
  n = 1 << log2n;

  av = cvec_new (l);
  bv = cvec_new_zeros (n);
  cv = cvec_new (n);

  /* copy and zero pad */
  for (i = 0; i < cvec_length (v); i++)
    cv[i] = v[i];

  for (i = cvec_length (v); i < n; i++)
    cv[i] = cplx_0;

  /* polar representation of z */
  rho = cnorm (z);
  theta = cang (z);

  /* elementwise multiply by z^{i^2/2} */
  bv[0] = cplx_1;
  s = 1.0 / n;

  for (i = 1; i < l; i++) {
    p = i * i / 2.0;
    crho = pow (rho, p);
    creal (c) = crho * cos (theta * p);
    cimag (c) = crho * sin (theta * p);
    d = cinv (c);
    bv[i] = d;
    bv[n - i] = d;		/* mirror needed for the convolution */
    av[i] = cscale (c, s);
    cv[i] = cmul (cv[i], c);
  }

  /* convolve in FFT domain by z^{-i^2/2} */
  fft_init ();
  _fft (cv, log2n);
  _fft (bv, log2n);
  cvec_mul (cv, bv);
  _ifft (cv, log2n);
  cvec_delete (bv);

  /* elementwise multiply by z^{i^2/2} */
  cv[0].r *= s;
  cv[0].i *= s;
  for (i = 1; i < l; i++)
    cv[i] = cmul (cv[i], av[i]);

  cvec_delete (av);
  cvec_set_length (cv, l);

  return (cv);
}

/* DFT through FFT */
#define FDFT_FORW 0
#define FDFT_BACK 1
static cvec it_fdft (cvec v, int inv)
{
  int  i, l, n;
  int  log2n;
  cplx c, d;
  cvec av, bv, cv;
  double p;
  double s;

  assert (cvec_length (v));

  l = cvec_length (v);

  log2n = intlog2 (l);
  if ((1 << log2n) != l)
    log2n++;
  log2n++;
  n = 1 << log2n;

  av = cvec_new (l);
  bv = cvec_new_zeros (n);
  cv = cvec_new (n);

  /* copy and zero pad */
  for (i = 0; i < cvec_length (v); i++)
    cv[i] = v[i];

  for (i = cvec_length (v); i < n; i++)
    cv[i] = cplx_0;

  /* elementwise multiply by z^{i^2/2} */
  bv[0] = cplx_1;
  if (inv)
    s = 1.0 / n / l;
  else
    s = 1.0 / n;

  for (i = 1; i < l; i++) {
    p = i * i / 2.0;
    creal (c) = creal (d) = cos (M_PI * i * i / l);
    if (inv) {
      cimag (c) = sin (M_PI * i * i / l);
      cimag (d) = -cimag (c);
    }
    else {
      cimag (c) = sin (-M_PI * i * i / l);
      cimag (d) = -cimag (c);
    }
    bv[i] = d;
    bv[n - i] = d;		/* mirror needed for the convolution */
    av[i] = cscale (c, s);
    cv[i] = cmul (cv[i], c);
  }

  /* convolve in FFT domain by z^{-i^2/2} */
  fft_init ();
  _fft (cv, log2n);
  _fft (bv, log2n);
  cvec_mul (cv, bv);
  _ifft (cv, log2n);
  cvec_delete (bv);

  /* elementwise multiply by z^{i^2/2} */
  cv[0].r *= s;
  cv[0].i *= s;
  for (i = 1; i < l; i++)
    cv[i] = cmul (cv[i], av[i]);

  cvec_delete (av);
  cvec_set_length (cv, l);

  return (cv);
}

/* Slow Discrete Fourier Transform */
static cvec it_sdft (cvec v)
{
  idx_t n, k, m;
  idx_t N;
  cvec cv;
  cvec roots;

  N = cvec_length (v);
  assert (N < DFT_LAST);
  cv = cvec_new (N);

  fft_init ();
  roots = dft_roots[N];
  for (k = 0; k < N; k++) {
    cv[k] = cplx_0;
    m = 0;
    for (n = 0; n < N; n++) {
      cv[k] = cadd (cv[k], cmul (v[n], roots[m]));
      m += N - k;
      m %= N;
    }
  }
  return cv;
}

/* Slow Inverse Discrete Fourier Transform */
static cvec it_sidft (cvec v)
{
  idx_t n, k, m;
  idx_t N;
  cvec cv;
  cvec roots;

  N = cvec_length (v);
  assert (N < DFT_LAST);
  cv = cvec_new (N);

  fft_init ();
  roots = dft_roots[N];
  for (n = 0; n < N; n++) {
    cv[n] = cplx_0;
    m = 0;
    for (k = 0; k < N; k++) {
      cv[n] = cadd (cv[n], cmul (v[k], roots[m]));
      m += n;
      m %= N;
    }
  }
  cvec_div_by_real (cv, N);
  return cv;
}

/* Discrete Fourier Transform */
cvec it_dft (cvec v)
{
  int  l;
  int  log2n;

  l = cvec_length (v);

  if (l == 0)
    return cvec_new (0);

  if (l == 1)
    return cvec_clone (v);

  log2n = intlog2 (l);

  if (l == (1 << log2n))
    return it_fft (v);

  if (l < DFT_LAST)
    return it_sdft (v);

  return it_fdft (v, FDFT_FORW);
}

/* Inverse Discrete Fourier Transform */
cvec it_idft (cvec v)
{
  int  l;
  int  log2n;

  l = cvec_length (v);

  if (l == 0)
    return cvec_new (0);

  if (l == 1)
    return cvec_clone (v);

  log2n = intlog2 (l);

  if (l == (1 << log2n))
    return it_ifft (v);

  if (l < DFT_LAST)
    return it_sidft (v);

  return it_fdft (v, FDFT_BACK);
}

/* Discrete Fourier Transform for real input data */
cvec it_dft_real (vec v)
{
  /* TODO: this is an ugly wrapper, optimized functions should be written */
  cvec cv = vec_to_cvec (v);
  cvec r = it_dft (cv);
  cvec_delete (cv);
  return r;
}

/* Inverse Discrete Fourier Transform for real input data */
vec it_idft_real (cvec cv)
{
  /* TODO: this is an ugly wrapper, optimized functions should be written */
  int  i, l;
  cvec r = it_idft (cv);
  vec  v;

  l = cvec_length (cv);
  v = vec_new (l);

  for (i = 0; i < l; i++)
    v[i] = creal (r[i]);

  cvec_delete (r);
  return v;
}

/*-------------- everything below will be deprecated in later releases ------*/
#if 0
/* do the real DFT in a clever way using the symmetries of the transform.  */
/* the real vector is packed in a complex vector of half the length of the */
/* real vector and the result of the transform is decoupled using a        */
/* mathematical relation involving the roots of unity.                     */
cvec __it_dft_real_fast (vec v)
{
  cvec cv, cvt;
  cvec roots;
  idx_t k, l;

  l = vec_length (v);
  assert (l % 2 == 0);		/* require an even-sized vector */

  cv = cvec_new (l / 2);

  /* transform the vec into a cvec as cv[k] = v[2k] + i v[2k+1] */
  /* Note: this could be done by casting if we ensure the structure is packed */
  /*       there are gcc options for this, I don't know about VC     [Vivien] */
  for (k = 0; k < l / 2; k++) {
    creal (cv[k]) = v[2 * k];
    cimag (cv[k]) = v[2 * k + 1];
  }

  /* call the mighty DFT */
  cvt = it_dft (cv);
  cvec_delete (cv);

  /* compute the conjugate of roots of unity times i */
  roots = cvec_new_unit_roots (l);
  for (k = 0; k < l; k++) {
    double r, i;
    r = creal (roots[k]);
    i = cimag (roots[k]);
    creal (roots[k]) = i;
    cimag (roots[k]) = r;
  }

  /* separate the result, we use the fact that [c.f. Numerical Recipies in C] */
  /* vt = 1/2 * [ (cvt[k] + cvt[N-k]^*) - i (cvt[k] - cvt[N-k]^*) W_N^k^* ]   */
  cvec_set_length (cvt, l + 1);	/* add one to put the low freq at the end */
  cvt[l / 2] = cvt[0];		/* copy the first element so that cvt[l/2] is defined */
  for (k = 0; k <= l / 2; k++) {
    cplx s = cconj (cvt[l / 2 - k]);	/* TODO insert cvt[0] at position N */
    cplx c = cvt[k];

    c = csub (cadd (c, s), cmul (csub (c, s), roots[k]));

    /* store in place */
    cvt[l - k] = cconj (c);
  }
  /* roots, bloody roots */
  cvec_delete (roots);

  /* now mirror the beginning */
  for (k = 0; k < l / 2; k++)
    cvt[k] = cconj (cvt[l - k]);
  cvec_set_length (cvt, l);	/* forget about the last coeff */

  /* the 1/2 factor in the expression of vt */
  cvec_div_by_real (cvt, 2);

  return (cvt);
}
#endif
