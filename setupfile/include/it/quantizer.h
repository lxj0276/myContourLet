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
  Base class for quantization.
  Vivien Chappelier <vivien.chappelier@irisa.fr>
*/

#ifndef __it_quantizer_h
#define __it_quantizer_h

#include <it/types.h>
#include <it/vec.h>
#include <it/mat.h>
#include <it/convcode.h> /* needed for TCQ */

#ifdef __cplusplus
extern "C" {
#endif

/* Generic vector quantizer class */
typedef struct _quantizer_ {
  it_extends(it_object_t);

  /* Quantize a vector v to a vector of indexes */
  ivec (* quantize)(struct _quantizer_ *quantizer, vec v);

  /* Dequantize the indices to the reconstruction vector */
  vec  (* dequantize)(struct _quantizer_ *quantizer, ivec Q);

  /* Return the cardinal of the codebook (0 means infinite or too large) */
  unsigned int  (* get_cardinal)(struct _quantizer_ *quantizer);

} it_quantizer_t;

#define IT_QUANTIZER(q) IT_CAST(it_quantizer_t, q)

static inline it_instanciate(it_quantizer_t)
{
  it_construct(it_object_t);
  it_set_magic(it_this, it_quantizer_t);
  return(it_this);
}

#define it_quantizer_quantize(it_this, v) __it_quantizer_quantize(IT_QUANTIZER(it_this), v)
static inline ivec __it_quantizer_quantize(it_quantizer_t *it_this, vec v)
{
  return(it_this->quantize(it_this, v));
}

#define it_quantizer_dequantize(it_this, i) __it_quantizer_dequantize(IT_QUANTIZER(it_this), i)
static inline vec __it_quantizer_dequantize(it_quantizer_t *it_this, ivec i)
{
  return(it_this->dequantize(it_this, i));
}

#define it_quantizer_get_cardinal(it_this) __it_quantizer_get_cardinal(IT_QUANTIZER(it_this))
static inline unsigned int __it_quantizer_get_cardinal(it_quantizer_t *it_this)
{
  return(it_this->get_cardinal(it_this));
}

/* Scalar quantizer */
typedef struct _scalar_quantizer_ {
  it_extends(it_quantizer_t);
  
  /* overload the virtual destructor */
  void (* it_overloaded(destructor))(it_object_t *it_this);

  /* set a new codebook for the scalar quantizer */
  void (* use_codebook)(struct _scalar_quantizer_ *quantizer, vec codebook, int first);

  /* return the quantizer codebook from index start to index end */
  vec (* get_codebook_range)(struct _scalar_quantizer_ *quantizer, int start, int end);

  /* get/set/clear the quantizer index range (inclusive) */
  void (* get_index_range)(struct _scalar_quantizer_ *quantizer, int * _imin, int * _imax);
  void (* set_index_range)(struct _scalar_quantizer_ *quantizer, int _imin, int _imax);

  /* quantize a scalar */
  int (* scalar_quantize)(struct _scalar_quantizer_ *quantizer, double v);

  /* dequantize an index */
  double (* scalar_dequantize)(struct _scalar_quantizer_ *quantizer, int i);

  /* quantizer codebook */
  vec codebook;

  /* valid index range */
  int imin;
  int imax;
  int first;

} it_scalar_quantizer_t;

#define IT_SCALAR_QUANTIZER(q) IT_CAST(it_scalar_quantizer_t, q)

it_instanciate(it_scalar_quantizer_t);

/* create a new quantizer from a codebook and the index corresponding to the first value in the codebook */
static inline it_scalar_quantizer_t *it_scalar_quantizer_new_start_index(vec codebook, int first)
{
  return(it_new_va(it_scalar_quantizer_t)(it_va, codebook, first));
}

#define it_scalar_quantizer_new(codebook) it_scalar_quantizer_new_start_index(codebook, 0)

/* set a new codebook for the scalar quantizer */
#define it_scalar_quantizer_use_codebook(it_this, codebook, first) \
          __it_scalar_quantizer_use_codebook(IT_SCALAR_QUANTIZER(it_this), codebook, first)
static inline void __it_scalar_quantizer_use_codebook(it_scalar_quantizer_t *it_this, vec codebook, int first) {
  it_this->use_codebook(it_this, codebook, first);
}

/* get the codebook of the quantizer from index 'start' to index 'end' */
#define it_scalar_quantizer_get_codebook_range(it_this, start, end) \
          __it_scalar_quantizer_get_codebook_range(IT_SCALAR_QUANTIZER(it_this), start, end)
static inline vec __it_scalar_quantizer_get_codebook_range(it_scalar_quantizer_t *it_this, int start, int end) {
  return(it_this->get_codebook_range(it_this, start, end));
}

#define it_scalar_quantizer_get_codebook(it_this) \
          it_scalar_quantizer_get_codebook_range(it_this, INT_MIN, INT_MAX)

#define it_scalar_quantizer_get_index_range(it_this, imin, imax) \
          __it_scalar_quantizer_get_index_range(IT_SCALAR_QUANTIZER(it_this), imin, imax)
static inline void __it_scalar_quantizer_get_index_range(it_scalar_quantizer_t *it_this, int *imin, int *imax) {
  it_this->get_index_range(it_this, imin, imax);
}

#define it_scalar_quantizer_set_index_range(it_this, imin, imax) \
          __it_scalar_quantizer_set_index_range(IT_SCALAR_QUANTIZER(it_this), imin, imax)
static inline void __it_scalar_quantizer_set_index_range(it_scalar_quantizer_t *it_this, int imin, int imax) {
  it_this->set_index_range(it_this, imin, imax);
}

#define it_scalar_quantizer_clear_index_range(it_this) \
          it_scalar_quantizer_set_index_range(it_this, INT_MIN, INT_MAX)

#define it_scalar_quantizer_scalar_quantize(it_this, v) \
          __it_scalar_quantizer_scalar_quantize(IT_SCALAR_QUANTIZER(it_this), v)
static inline int __it_scalar_quantizer_scalar_quantize(it_scalar_quantizer_t *it_this, double v)
{
  return(it_this->scalar_quantize(it_this, v));
}

#define it_scalar_quantizer_scalar_dequantize(it_this, i) \
          __it_scalar_quantizer_scalar_dequantize(IT_SCALAR_QUANTIZER(it_this), i)
static inline double __it_scalar_quantizer_scalar_dequantize(it_scalar_quantizer_t *it_this, int i)
{
  return(it_this->scalar_dequantize(it_this, i));
}

/* Uniform scalar quantizer */
typedef struct _uniform_quantizer_ {
  it_extends(it_scalar_quantizer_t);

  /* quantizer step size */
  double step;

  /* quantizer center */
  double center;

  /* range */
  double min;
  double max;

  /* dead zone scale factor */
  double factor;

  /* dead zone step (=factor*step) */
  double deadzone; 

} it_uniform_quantizer_t;

it_instanciate(it_uniform_quantizer_t);

#define IT_UNIFORM_QUANTIZER(q) IT_CAST(it_uniform_quantizer_t, q)

/* constructor from center and step */
static inline it_uniform_quantizer_t *it_uniform_quantizer_new_from_center(double center, double step, double factor)
{
  return(it_new_va(it_uniform_quantizer_t)(it_va, center, step, factor));
}

/* contructor from number of codewords and range */
static inline it_uniform_quantizer_t *it_uniform_quantizer_new_from_range(int n, double min, double max)
{
  double center;
  double step;
  it_uniform_quantizer_t *quantizer;

  step = (max - min) / n;
  center = min + step / 2;
  quantizer = it_uniform_quantizer_new_from_center(center, step, 1.0);
  it_scalar_quantizer_set_index_range(IT_SCALAR_QUANTIZER(quantizer), 0, n - 1);

  return(quantizer);
}

/* get the step size of the quantizer */
#define it_uniform_quantizer_get_step(q) (IT_UNIFORM_QUANTIZER(q)->step)

/* get the center of the quantizer */
#define it_uniform_quantizer_get_center(q) (IT_UNIFORM_QUANTIZER(q)->center)

/* get the dead zone step factor */
#define it_uniform_quantizer_get_deadzone_factor(q) (IT_UNIFORM_QUANTIZER(q)->factor)

/* get the dead zone step */
#define it_uniform_quantizer_get_deadzone_step(q)   (IT_UNIFORM_QUANTIZER(q)->deadzone)

/* TCQ class */
typedef struct _trellis_coded_quantizer_ {
  it_extends(it_quantizer_t);

  /* convolutional code used for set partitioning */
  it_convolutional_code_t *code;

  /* scalar quantizers for each coset */
  it_scalar_quantizer_t **coset_quantizers;

  /* overload the virtual destructor */
  void (* it_overloaded(destructor))(it_object_t *it_this);

} it_trellis_coded_quantizer_t;

it_instanciate(it_trellis_coded_quantizer_t);

#define IT_TRELLIS_CODED_QUANTIZER(q) IT_CAST(it_trellis_coded_quantizer_t, q)

/* Initialize a TCQ from a quantizer by partitioning it */
#define it_trellis_coded_quantizer_new_partition(code, quantizer) \
          __it_trellis_coded_quantizer_new_partition(IT_CONVOLUTIONAL_CODE(code), IT_SCALAR_QUANTIZER(quantizer))
static inline it_trellis_coded_quantizer_t *__it_trellis_coded_quantizer_new_partition(it_convolutional_code_t *code,
										       it_scalar_quantizer_t *quantizer)
{
  return(it_new_va(it_trellis_coded_quantizer_t)(it_va, code, quantizer));  
}

/* Matrix quantization: each line is quantized independently */
#define it_quantize_mat(q, m) __it_quantize_mat(IT_QUANTIZER(q), m)
imat __it_quantize_mat(it_quantizer_t *q, mat m);
#define it_dequantize_mat(q, m) __it_dequantize_mat(IT_QUANTIZER(q), m)
mat __it_dequantize_mat(it_quantizer_t *q, imat m);

/* shortnames for some functions */
#define it_quantize(q, v) it_scalar_quantizer_scalar_quantize(q, v)
#define it_dequantize(q, i) it_scalar_quantizer_scalar_dequantize(q, i)
#define it_quantize_vec(q, v) it_quantizer_quantize(q, v)
#define it_dequantize_vec(q, v) it_quantizer_dequantize(q, v)

#define it_quantizer_use_codebook(q, codebook, start) \
          it_scalar_quantizer_use_codebook(q, codebook, start)
#define it_quantizer_get_codebook_range(q, start, end) \
          it_scalar_quantizer_get_codebook_range(q, start, end)
#define it_quantizer_get_codebook(q) \
          it_scalar_quantizer_get_codebook(q)
#define it_quantizer_get_index_range(q, imin, imax) \
          it_scalar_quantizer_get_index_range(q, imin, imax)
#define it_quantizer_set_index_range(q, imin, imax) \
          it_scalar_quantizer_set_index_range(q, imin, imax)
#define it_quantizer_clear_index_range(q) \
          it_scalar_quantizer_clear_index_range(q)
#define it_quantizer_scalar_quantize(q, v) \
          it_scalar_quantizer_scalar_quantize(q, v)
#define it_quantizer_scalar_dequantize(q, i) \
          it_scalar_quantizer_scalar_dequantize(q, i)
#define it_quantizer_get_step(q) it_uniform_quantizer_get_step(q)
#define it_quantizer_get_center(q) it_uniform_quantizer_get_center(q)
#define it_quantizer_get_deadzone_factor(q) \
          it_uniform_quantizer_get_deadzone_factor(q) 
#define it_quantizer_get_deadzone_step(q) \
          it_uniform_quantizer_get_deadzone_step(q) 

/* helper functions (to be moved to quantizer_func.h */

/* generate the codebook for the Lloyd-Max quantizer on N levels  */
/* This is the optimal codebook for a stationnary i.i.d. source.  */
/* [a,b] is the interval where to expect the pdf to be non-zero   */
/* a good guess is [ mean - std_dev, mean + std_dev ].            */
/* Note that returned centroids may lie outside this interval     */
/* The algorithm stops when the relative difference in distortion */
/* is smaller than IT_EPSILON.                                    */

vec lloyd_max(it_function_t function, it_args_t args,
	      double a, double b, int N);

#ifdef __cplusplus
}
#endif /* extern "C" */
#endif
