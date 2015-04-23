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
  Base class for convolutional codes.
  Vivien Chappelier <vivien.chappelier@irisa.fr>
*/

#ifndef __it_convcode_h
#define __it_convcode_h

#include <it/types.h>
#include <it/mat.h>
#include <it/vec.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Generic convolutional code class */
typedef struct _convolutional_code_ {
  it_extends(it_object_t);
  
  /* overload the virtual destructor */
  void (* it_overloaded(destructor))(it_object_t *it_this);

  /* Encode a binary vector with the convolutional code */
  bvec (* encode)(struct _convolutional_code_ *convolutional_code, bvec v);

  /* Encode a 2^k-valued vector to a 2^n-valued vector */
  ivec (* encode_symbolic)(struct _convolutional_code_ *cc, ivec b);

  unsigned int k;                 /* input bits */
  unsigned int n;                 /* output bits */
  unsigned int memory;            /* memory of the code */
  unsigned int n_labels;
  unsigned int n_states;
  unsigned int n_branches;
  imat next;                      /* next state   [n_states][n_branches] */
  imat label;                     /* branch label [n_states][n_branches] */

  /*  unsigned long Q;*/    /* feedback polynomial Q (0 is non-systematic) */

  imat generators;    /* generator polynomials */ 

} it_convolutional_code_t;

#define IT_CONVOLUTIONAL_CODE(q) IT_CAST(it_convolutional_code_t, q)

it_instanciate(it_convolutional_code_t);

/* create a new convolutional code from its generator polynomials */
static inline it_convolutional_code_t *it_convolutional_code_new(imat generators)
{
  return(it_new_va(it_convolutional_code_t)(it_va, generators));
}

#define it_convolutional_code_encode(it_this, b) __it_convolutional_code_encode(IT_CONVOLUTIONAL_CODE(it_this), b)
static inline bvec __it_convolutional_code_encode(it_convolutional_code_t *it_this, bvec b)
{
  return(it_this->encode(it_this, b));
}

#define it_convolutional_code_encode_symbolic(it_this, b) __it_convolutional_code_encode_symbolic(IT_CONVOLUTIONAL_CODE(it_this), b)
static inline ivec __it_convolutional_code_encode_symbolic(it_convolutional_code_t *it_this, ivec b)
{
  return(it_this->encode_symbolic(it_this, b));
}

#define it_convolutional_code_trellis_next(cc, s, b) ((cc)->next[s][b])
#define it_convolutional_code_trellis_label(cc, s, b) ((cc)->label[s][b])

/* shortnames for some functions */
#define it_cc_encode(cc, b) it_convolutional_code_encode(cc, b)
#define it_cc_encode_symbolic(cc, b) it_convolutional_code_encode_symbolic(cc, b)
#define it_cc_decode(cc, v) it_convolutional_code_decode(cc, v)
#define it_cc_next(cc, s, b)  it_convolutional_code_trellis_next(cc, s, b)
#define it_cc_label(cc, s, b) it_convolutional_code_trellis_label(cc, s, b)

/* helper functions */
/* TODO: move to convcode_func.h */
ivec it_viterbi_decode_symbolic(it_convolutional_code_t *cc, mat metrics);
bvec it_viterbi_decode(it_convolutional_code_t *cc, mat metrics);

#ifdef __cplusplus
}
#endif /* extern "C" */
#endif


