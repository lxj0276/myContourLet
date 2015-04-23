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
  Quantizers
  Copyright (C) 2005 Vivien Chappelier
*/

#include <stdarg.h>
#include <it/types.h>
#include <it/vec.h>
#include <it/math.h>
#include <it/quantizer.h>
#include <it/io.h>

/* function prototypes */
static void scalar_quantizer_destructor (it_object_t * it_this);
static vec scalar_quantizer_get_codebook_range (it_scalar_quantizer_t * q,
						int start, int end);
static void scalar_quantizer_get_index_range (it_scalar_quantizer_t * q,
					      int *_imin, int *_imax);
static void scalar_quantizer_set_index_range (it_scalar_quantizer_t * q,
					      int _imin, int _imax);
static int scalar_quantizer_scalar_quantize (it_scalar_quantizer_t * q,
					     double v);
static double scalar_quantizer_scalar_dequantize (it_scalar_quantizer_t *
						  it_this, int q);
static void scalar_quantizer_set_codebook (it_scalar_quantizer_t * it_this,
					   vec codebook, int first);
static ivec scalar_quantizer_quantize (it_quantizer_t * it_this, vec v);
static vec scalar_quantizer_dequantize (it_quantizer_t * it_this, ivec q);
static unsigned int scalar_quantizer_get_cardinal (it_quantizer_t * it_this);

static int uniform_quantizer_scalar_quantize (it_scalar_quantizer_t *
					      it_this, double v);
static double uniform_quantizer_scalar_dequantize (it_scalar_quantizer_t *
						   it_this, int q);
static vec uniform_quantizer_get_codebook_range (it_scalar_quantizer_t *
						 it_this, int s, int e);
static void uniform_quantizer_set_index_range (it_scalar_quantizer_t * q,
					       int _imin, int _imax);
static void uniform_quantizer_set_codebook (it_scalar_quantizer_t * it_this,
					    vec codebook, int first);

static void trellis_coded_quantizer_destructor (it_object_t * it_this);
static ivec trellis_coded_quantizer_quantize (it_quantizer_t * it_this,
					      vec v);
static vec trellis_coded_quantizer_dequantize (it_quantizer_t * it_this,
					       ivec q);
static unsigned int trellis_coded_quantizer_get_cardinal (it_quantizer_t *
							  it_this);


/* instanciation of a scalar quantizer */
it_instanciate (it_scalar_quantizer_t)
{
  vec  codebook;

  /* start variable argument list */
  it_new_args_start ();

  /* call the parent constructor */
  it_construct (it_quantizer_t);

  /* set the appropriate magic */
  it_set_magic (it_this, it_scalar_quantizer_t);

  /* overload the virtual destructor */
  it_overload (it_this, it_object_t, destructor, scalar_quantizer_destructor);

  /* assign methods */
  IT_QUANTIZER (it_this)->quantize = scalar_quantizer_quantize;
  IT_QUANTIZER (it_this)->dequantize = scalar_quantizer_dequantize;
  IT_QUANTIZER (it_this)->get_cardinal = scalar_quantizer_get_cardinal;
  it_this->set_codebook = scalar_quantizer_set_codebook;
  it_this->get_codebook_range = scalar_quantizer_get_codebook_range;
  it_this->get_index_range = scalar_quantizer_get_index_range;
  it_this->set_index_range = scalar_quantizer_set_index_range;
  it_this->scalar_quantize = scalar_quantizer_scalar_quantize;
  it_this->scalar_dequantize = scalar_quantizer_scalar_dequantize;

  /* construction */
  codebook = it_new_args_next (vec);
  it_this->first = it_new_args_next (int);

  if (codebook)
    scalar_quantizer_set_codebook (it_this, codebook, it_this->first);
  else {
    it_this->codebook = NULL;
    it_scalar_quantizer_set_index_range (it_this, 0, -1);
  }

  it_new_args_stop ();
  return (it_this);
}

static void scalar_quantizer_destructor (it_object_t * it_this)
{
  it_scalar_quantizer_t *scalar_quantizer = IT_SCALAR_QUANTIZER (it_this);

  if (scalar_quantizer->codebook)
    vec_delete (scalar_quantizer->codebook);

  /* call the parent destructor */
  scalar_quantizer->it_overloaded (destructor) (it_this);
}

static vec scalar_quantizer_get_codebook_range (it_scalar_quantizer_t * q,
						int start, int end)
{
  if (start < q->imin)
    start = q->imin;
  if (end > q->imax)
    end = q->imax;
  return (vec_get_subvector (q->codebook, start, end));
}

static void scalar_quantizer_get_index_range (it_scalar_quantizer_t * q,
					      int *_imin, int *_imax)
{
  *_imax = q->imax;
  *_imin = q->imin;
}

static void scalar_quantizer_set_index_range (it_scalar_quantizer_t * q,
					      int _imin, int _imax)
{
  if (_imin < _imax) {
    q->imax = _imax;
    q->imin = _imin;
  }
  else {
    q->imax = _imin;
    q->imin = _imax;
  }

  if (q->imin < q->first)
    q->imin = q->first;
  if (q->codebook) {
    if (q->imax > q->first + (int) vec_length (q->codebook) - 1)
      q->imax = q->first + (int) vec_length (q->codebook) - 1;
  }
  else {
    q->imax = q->imin;
  }
}

/* scalar quantization */
static int scalar_quantizer_scalar_quantize (it_scalar_quantizer_t * q,
					     double v)
{
  int  s, e, h;
  double c;
  int  first = q->first;
  vec  codebook = q->codebook;

  /* logarithmic search for the index of the closest vector */
  s = q->imin - first;
  e = q->imax - first;

  while (e > s + 1) {
    /* compare with the median of the range */
    h = (s + e) / 2;
    c = codebook[h];

    /* narrow down the range */
    if (v < c)
      e = h;
    else
      s = h;
  }

  /* find the closest vector */
  if (v - codebook[s] < codebook[e] - v)
    return (s + first);
  else
    return (e + first);
}

/* scalar dequantization */
static double scalar_quantizer_scalar_dequantize (it_scalar_quantizer_t *
						  it_this, int q)
{
  int  first = it_this->first;
  return (it_this->codebook[q - first]);
}

/* change the codebook */
static void scalar_quantizer_set_codebook (it_scalar_quantizer_t * it_this,
					   vec codebook, int first)
{
  if (!codebook) {
    fprintf (stderr,
	     "scalar_quantizer_set_codebook called with a NULL codebook\n");
    abort ();
  }

  it_this->codebook = vec_clone (codebook);
  it_this->first = first;
  it_scalar_quantizer_set_index_range (it_this, first,
				       first + vec_length (codebook) - 1);
  /* sort the codebook in increasing order */
  vec_sort (it_this->codebook);
}

/* quantize the vector v */
static ivec scalar_quantizer_quantize (it_quantizer_t * it_this, vec v)
{
  it_scalar_quantizer_t *scalar_quantizer = IT_SCALAR_QUANTIZER (it_this);
  int  i;
  int  l = vec_length (v);
  ivec ret = ivec_new (l);

  for (i = 0; i < l; i++)
    ret[i] = it_scalar_quantizer_scalar_quantize (scalar_quantizer, v[i]);

  return (ret);
}

/* dequantize the indices vector q */
static vec scalar_quantizer_dequantize (it_quantizer_t * it_this, ivec q)
{
  it_scalar_quantizer_t *scalar_quantizer = IT_SCALAR_QUANTIZER (it_this);
  int  i;
  int  l = ivec_length (q);
  vec  ret = vec_new (l);

  for (i = 0; i < l; i++)
    ret[i] = it_scalar_quantizer_scalar_dequantize (scalar_quantizer, q[i]);

  return (ret);
}

/* get the number of codewords for the uniform scalar quantizer */
static unsigned int scalar_quantizer_get_cardinal (it_quantizer_t * it_this)
{
  int  imin = IT_SCALAR_QUANTIZER (it_this)->imin;
  int  imax = IT_SCALAR_QUANTIZER (it_this)->imax;

  return (imax - imin + 1);	/* may overflow, that's ok, 0 means 'infinite' */
}

/*---------------------- Uniform Quantizer -----------------------------*/

it_instanciate (it_uniform_quantizer_t)
{
  it_new_args_start ();

  /* construct a scalar quantizer with empty codebook */
  it_construct_va (it_scalar_quantizer_t) (it_va, NULL, 0);
  it_set_magic (it_this, it_uniform_quantizer_t);

  /* assign/overload methods */
  IT_SCALAR_QUANTIZER (it_this)->scalar_quantize =
    uniform_quantizer_scalar_quantize;
  IT_SCALAR_QUANTIZER (it_this)->scalar_dequantize =
    uniform_quantizer_scalar_dequantize;
  IT_SCALAR_QUANTIZER (it_this)->get_codebook_range =
    uniform_quantizer_get_codebook_range;
  IT_SCALAR_QUANTIZER (it_this)->set_index_range =
    uniform_quantizer_set_index_range;
  IT_SCALAR_QUANTIZER (it_this)->set_codebook =
    uniform_quantizer_set_codebook;

  /* construction */
  it_this->center = it_new_args_next (double);
  it_this->step = it_new_args_next (double);
  it_this->factor = it_new_args_next (double);
  it_this->deadzone = it_this->step * it_this->factor;
  it_scalar_quantizer_clear_index_range (it_this);

  it_new_args_stop ();
  return (it_this);
}

static int uniform_quantizer_scalar_quantize (it_scalar_quantizer_t *
					      it_this, double v)
{
  int  imin = it_this->imin;
  int  imax = it_this->imax;
  double min = IT_UNIFORM_QUANTIZER (it_this)->min;
  double max = IT_UNIFORM_QUANTIZER (it_this)->max;
  double center = it_quantizer_get_center (it_this);
  double step = it_quantizer_get_step (it_this);
  double dead_step = it_quantizer_get_deadzone_step (it_this);

  if (v <= min)
    return (imin);
  if (v >= max)
    return (imax);

  v -= center;			/* center the value */

  if (v > 0)			/* divide and round */
    return ((int) ((v + step - dead_step / 2) / step));
  else
    return ((int) ((v - step + dead_step / 2) / step));
}

/* uniform scalar dequantization */
static double uniform_quantizer_scalar_dequantize (it_scalar_quantizer_t *
						   it_this, int q)
{
  double center = it_quantizer_get_center (it_this);
  double step = it_quantizer_get_step (it_this);
  double dead_step = it_quantizer_get_deadzone_step (it_this);

  if (q == 0)
    return (center);

  if (q > 0)
    return (center + q * step + (dead_step - step) / 2);
  else
    return (center + q * step - (dead_step - step) / 2);
}

/* dynamically generate the codebook... use with caution */
static vec uniform_quantizer_get_codebook_range (it_scalar_quantizer_t *
						 it_this, int s, int e)
{
  int  imin = it_this->imin;
  int  imax = it_this->imax;
  vec  book;
  int  i;

  if (it_this->codebook)
    return (it_scalar_quantizer_get_codebook_range (it_this, s, e));

  if (s < imin)
    s = imin;
  if (e > imax)
    e = imax;
  book = vec_new (e - s + 1);

  for (i = s; i <= e; i++)
    book[i - s] = it_scalar_quantizer_scalar_dequantize (it_this, i);

  return (book);
}

static void uniform_quantizer_set_index_range (it_scalar_quantizer_t *
					       it_this, int _imin, int _imax)
{
  if (_imin < _imax) {
    it_this->imax = _imax;
    it_this->imin = _imin;
  }
  else {
    it_this->imax = _imin;
    it_this->imin = _imax;
  }

  IT_UNIFORM_QUANTIZER (it_this)->min = it_dequantize (it_this, _imin);
  IT_UNIFORM_QUANTIZER (it_this)->max = it_dequantize (it_this, _imax);
}

/* change the codebook */
static void uniform_quantizer_set_codebook (it_scalar_quantizer_t * it_this,
					    vec codebook, int first)
{
  fprintf (stderr, "cannot change the codebook of a uniform quantizer\n");
  abort ();
}

/*----------------- Trellis Coded Quantizer ----------------------------*/

/* gray code to binary code */
#define gray_to_bin(x) ((x) ^ (x >> 1))

it_instanciate (it_trellis_coded_quantizer_t)
{
  int  i, j, k;
  int  n, r, m;

  vec  codebook;
  vec  subcodebook;
  it_convolutional_code_t *code;
  it_scalar_quantizer_t *quantizer;
  it_scalar_quantizer_t **coset_quantizers;

  it_new_args_start ();

  it_construct (it_quantizer_t);
  it_set_magic (it_this, it_trellis_coded_quantizer_t);

  /* assign/overload methods */
  it_overload (it_this, it_object_t, destructor,
	       trellis_coded_quantizer_destructor);
  IT_QUANTIZER (it_this)->quantize = trellis_coded_quantizer_quantize;
  IT_QUANTIZER (it_this)->dequantize = trellis_coded_quantizer_dequantize;
  IT_QUANTIZER (it_this)->get_cardinal = trellis_coded_quantizer_get_cardinal;

  /* construction */
  code = it_new_args_next (it_convolutional_code_t *);
  quantizer = it_new_args_next (it_scalar_quantizer_t *);
  it_this->code = code;		/* TODO: deep copy? */

  /* retrieve the codebook */
  codebook = it_quantizer_get_codebook (quantizer);
  n = vec_length (codebook);
  r = code->n;
  m = (1 << r);

  /* set partitioning of the codebook */
  it_this->coset_quantizers = coset_quantizers =
    (it_scalar_quantizer_t **) malloc (m * sizeof (it_scalar_quantizer_t *));

  for (j = 0; j < m; j++) {
    subcodebook = vec_new ((n + m - j - 1) / m);

    for (i = j, k = 0; i < n; i += m, k++)
      subcodebook[k] = codebook[i];

    coset_quantizers[gray_to_bin (j)] = it_scalar_quantizer_new (subcodebook);
    vec_delete (subcodebook);
  }

  it_new_args_stop ();
  return (it_this);
}

/* TCQ destructor */
static void trellis_coded_quantizer_destructor (it_object_t * it_this)
{
  it_trellis_coded_quantizer_t *tcq = IT_TRELLIS_CODED_QUANTIZER (it_this);
  it_scalar_quantizer_t **coset_quantizers = tcq->coset_quantizers;
  int  m = 1 << tcq->code->n;	/* number of cosets */
  int  i;

  for (i = 0; i < m; i++)
    it_delete (coset_quantizers[i]);

  /* TODO: destruct the code? (if deep copy) */

  /* call the parent destructor */
  tcq->it_overloaded (destructor) (it_this);
}

/* quantize the vector v */
static ivec trellis_coded_quantizer_quantize (it_quantizer_t * it_this, vec v)
{
  it_trellis_coded_quantizer_t *tcq = IT_TRELLIS_CODED_QUANTIZER (it_this);
  it_scalar_quantizer_t **coset_quantizers = tcq->coset_quantizers;
  it_convolutional_code_t *cc = tcq->code;
  mat  metrics;
  imat bestword;
  int  i, j, q;
  double r;
  int  l = vec_length (v);
  int  n_labels = cc->n_labels;
  ivec ret = ivec_new (l);
  ivec path;
  ivec books;

  /* fill the metrics and best words */
  metrics = mat_new (n_labels, l);
  bestword = imat_new (n_labels, l);

  for (i = 0; i < l; i++)
    for (j = 0; j < n_labels; j++) {
      /* quantize in the coset associated with the branch label */
      q = it_quantize (coset_quantizers[j], v[i]);
      r = it_dequantize (coset_quantizers[j], q);
      /* compute the distortion associated to that label */
      metrics[j][i] = -(r - v[i]) * (r - v[i]);
      bestword[j][i] = q;
    }

  /* search for the sequence of minimal distortion */
  path = it_viterbi_decode_symbolic (cc, metrics);

  /* get the sequence of codebooks */
  books = it_cc_encode_symbolic (cc, path);

  /* merge the index and path bits to codeword bits */
  for (i = 0; i < l; i++)
    ret[i] = (bestword[books[i]][i] << 1) | path[i];

  mat_delete (metrics);
  imat_delete (bestword);

  ivec_delete (path);
  ivec_delete (books);

  return (ret);
}

/* dequantize the indices vector q */
static vec trellis_coded_quantizer_dequantize (it_quantizer_t * it_this,
					       ivec q)
{
  it_trellis_coded_quantizer_t *tcq = IT_TRELLIS_CODED_QUANTIZER (it_this);
  it_scalar_quantizer_t **coset_quantizers = tcq->coset_quantizers;
  it_convolutional_code_t *cc = tcq->code;
  int  i;
  int  l = ivec_length (q);
  vec  ret = vec_new (l);
  ivec path;			/* path bits */
  ivec books;			/* sequence of codebooks */
  int  word;

  /* split the index into path bits and codeword bits */
  path = ivec_new (l);
  for (i = 0; i < l; i++)
    path[i] = q[i] & 1;

  /* get the sequence of codebooks */
  books = it_cc_encode_symbolic (cc, path);

  /* dequantize using the appropriate codebook */
  for (i = 0; i < l; i++) {
    word = q[i] >> 1;
    ret[i] = it_dequantize (coset_quantizers[books[i]], word);
  }

  ivec_delete (path);
  ivec_delete (books);
  return (ret);
}


/* get the number of codewords for the TCQ quantizer */
static unsigned int trellis_coded_quantizer_get_cardinal (it_quantizer_t *
							  it_this)
{
  it_trellis_coded_quantizer_t *tcq = IT_TRELLIS_CODED_QUANTIZER (it_this);
  it_convolutional_code_t *cc = tcq->code;
  int  k = cc->k;		/* number of path bits */
  int  c = it_quantizer_get_cardinal (tcq->coset_quantizers[0]);	/* number of codewords */

  return ((1 << k) * c);
}


/*----------------- Helper functions ----------------------------*/

/* used to compute the distortion in Lloyd-Max */
it_function_args (lloyd_max_distortion)
{
  double expectation;		/* reconstruction value */
  it_function_t function;	/* function to integrate */
  it_args_t args;		/* and its arguments */
};

it_function (lloyd_max_distortion)
{
  double E = it_this->expectation;
  it_function_t f = it_this->function;
  it_args_t args = it_this->args;
  double d;

  d = x - E;

  return (d * d * f (x, args));
}


/* find the optimal codebook for a stationnary i.i.d. source */
/* the algorithm stops when the relative difference in distortion */
/* is smaller than IT_EPSILON. */
vec lloyd_max (it_function_t function, it_args_t args,
	       double a, double b, int N)
{
  /* for the initial guess, use a scalar quantizer of step (b-a)/N */
  double step = (b - a) / N;
  vec  threshold = vec_new_arithm (a, step, N + 1);
  vec  codebook = vec_new_arithm (a, step, N);
  int  i;
  double density;
  double expectation;
  double distortion, old_distortion;
  it_function_args (itf_integrate) density_args;
  it_function_args (itf_expectation) expectation_args;
  it_function_args (lloyd_max_distortion) distortion_args;
  it_function_args (itf_integrate) integrate_args;

  /* for integrating f(x) */
  density_args.function = function;
  density_args.args = args;
  /* for integrating x * f(x) */
  expectation_args.function = function;
  expectation_args.args = args;
  /* for integrating (x - E[X])^2 * f(x) */
  distortion_args.function = function;
  distortion_args.args = args;
  integrate_args.function = lloyd_max_distortion;
  integrate_args.args = &distortion_args;

  distortion = HUGE_VAL;

  do {

    old_distortion = distortion;
    distortion = 0;

    /* compute centroids */
    for (i = 0; i < N; i++) {
      /* integrate f(x) between threshold[i] and threshold[i+1] */
      density_args.a = threshold[i];
      density = itf_integrate (threshold[i + 1], &density_args);

      /* integrate xf(x) between threshold[i] and threshold[i+1] */
      expectation_args.a = threshold[i];
      expectation = itf_expectation (threshold[i + 1], &expectation_args);

      codebook[i] = expectation / density;

      /* integrate (x - E)^2 f(x) */
      integrate_args.a = threshold[i];
      distortion_args.expectation = codebook[i];
      distortion += itf_integrate (threshold[i + 1], &integrate_args);
    }

    /* compute thresholds */
    for (i = 0; i < N - 1; i++)
      threshold[i + 1] = (codebook[i] + codebook[i + 1]) / 2.0;

  }
  while ((old_distortion - distortion) / distortion > IT_EPSILON);

  return (codebook);
}

/* Matrix quantization: each line is quantized independently */
imat __it_quantize_mat (it_quantizer_t * q, mat m)
{
  idx_t i, h;
  imat r;
  h = mat_height (m);
  r = (ivec *) Vec_new (ivec, h);
  for (i = 0; i < h; i++)
    r[i] = it_quantize_vec (q, m[i]);
  return (r);
}

/* Matrix dequantization: each line is dequantized independently */
mat __it_dequantize_mat (it_quantizer_t * q, imat m)
{
  idx_t i, h;
  mat  r;
  h = imat_height (m);
  r = (vec *) Vec_new (vec, h);
  for (i = 0; i < h; i++)
    r[i] = it_dequantize_vec (q, m[i]);
  return (r);
}
