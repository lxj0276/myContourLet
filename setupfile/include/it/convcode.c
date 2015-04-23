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
  Convolutional codes
  Copyright (C) 2005 Vivien Chappelier
*/

#include <it/types.h>
#include <it/vec.h>
#include <it/convcode.h>
#include <it/io.h>
#include <it/hmmalgo.h>

#define bit(x, n) (((x) >> (n)) & 1)

#define DECLARE_CONV_CODE(x, k, n) \
  int it_cc_generators_##x[k*n]

/* some convolutional and recursive convolutional codes */
DECLARE_CONV_CODE (conv_code_1, 1, 2) = {
02, 03};

DECLARE_CONV_CODE (conv_code_2, 1, 2) = {
05, 07};

DECLARE_CONV_CODE (conv_code_3, 1, 2) = {
015, 017};

DECLARE_CONV_CODE (conv_code_4, 1, 2) = {
023, 033};

DECLARE_CONV_CODE (conv_code_6, 1, 2) = {
0133, 0171};

/* local helper functions */
static int fls (int c);
static int encode_conv (it_convolutional_code_t * cc, int *state, int input,
			imat G, int Q, int k, int n);

/* function prototypes */
static void convolutional_code_destructor (it_object_t * it_this);
static bvec convolutional_code_encode (it_convolutional_code_t * cc, bvec b);
static ivec convolutional_code_encode_symbolic (it_convolutional_code_t * cc,
						ivec b);

/* instanciation of a convolutional code */
it_instanciate (it_convolutional_code_t)
{
  imat generators;
  int  Q;
  int  i, j, c;
  int  s, n, k, l;
  int  state, symbol;
  int  n_states, n_symbols, n_labels;

  it_new_args_start ();

  it_construct (it_object_t);
  it_set_magic (it_this, it_convolutional_code_t);

  it_overload (it_this, it_object_t, destructor,
	       convolutional_code_destructor);

  /* assign methods */
  it_this->encode = convolutional_code_encode;
  it_this->encode_symbolic = convolutional_code_encode_symbolic;

  /* construction */
  it_this->generators = generators = imat_clone (it_new_args_next (imat));
  it_this->Q = Q = it_new_args_next (int);
  it_this->k = k = imat_height (generators);
  it_this->n = n = imat_width (generators);

  /* find the total memory of the code */
  l = 0;
  for (j = 0; j < k; j++) {
    c = Q;
    for (i = 0; i < n; i++)
      c |= generators[j][i];
    l += fls (c) - 1;
  }

  /* allocate automaton */
  it_this->memory = l;
  it_this->n_states = n_states = 1 << l;
  it_this->n_symbols = n_symbols = 1 << k;
  it_this->n_labels = n_labels = 1 << n;

  it_this->automaton = imat_new_set (DUMM_NODE, 1 << n, n_states);
  it_this->output = imat_new (1 << k, n_states);
  it_this->input = imat_new_set (NULL_INDEX, 1 << n, n_states);

  /* start from the zero state */
  it_this->logalpha_0 = vec_new_set (-HUGE_VAL, n_states);
  it_this->logalpha_0[0] = 0.0;

  /* end in whatever state */
  it_this->logbeta_end = vec_new_zeros (n_states);

  /* each branch has equal a priori probability of being taken */
  it_this->next_state_logpt = mat_new_zeros (n_labels, n_states);

  /* compute label and next state for each state and symbol */
  for (state = 0; state < n_states; state++) {
    for (symbol = 0; symbol < n_symbols; symbol++) {
      s = state;
      c = encode_conv (it_this, &s, symbol, generators, Q, k, n);

      it_this->automaton[c][state] = s;
      it_this->output[symbol][state] = c;
      it_this->input[c][state] = symbol;
    }
  }

  it_new_args_stop ();
  return (it_this);
}


static void convolutional_code_destructor (it_object_t * it_this)
{
  it_convolutional_code_t *convolutional_code =
    IT_CONVOLUTIONAL_CODE (it_this);

  if (convolutional_code->generators)
    imat_delete (convolutional_code->generators);

  vec_delete (convolutional_code->logalpha_0);
  vec_delete (convolutional_code->logbeta_end);
  mat_delete (convolutional_code->next_state_logpt);
  imat_delete (convolutional_code->automaton);
  imat_delete (convolutional_code->input);
  imat_delete (convolutional_code->output);

  /* call the parent destructor */
  convolutional_code->it_overloaded (destructor) (it_this);
}


static bvec convolutional_code_encode (it_convolutional_code_t * cc, bvec b)
{
  int  i, j;
  int  k = cc->k;
  int  n = cc->n;
  idx_t L = bvec_length (b) / k;
  bvec r = bvec_new (L * n);
  int  state = 0;
  int  symbol;
  int  label;

  assert (L * k == bvec_length (b));

  for (i = 0; i < L; i++) {
    symbol = 0;
    /* read k bits */
    for (j = 0; j < k; j++) {
      assert (!(b[i * k + j] >> 1));
      symbol |= b[i * k + j] << (k - 1 - j);
    }

    /* get the label and next state for that symbol */
    label = it_cc_label (cc, state, symbol);
    state = it_cc_next (cc, state, symbol);

    /* write n bits */
    for (j = 0; j < n; j++)
      r[i * n + j] = bit (label, n - 1 - j);
  }

  return (r);
}


static ivec convolutional_code_encode_symbolic (it_convolutional_code_t * cc,
						ivec s)
{
  int  i;
  int  L = ivec_length (s);
  ivec r = ivec_new (L);
  int  state = 0;

  for (i = 0; i < L; i++) {
    /* get the label and next state for that symbol */
    r[i] = it_cc_label (cc, state, s[i]);
    state = it_cc_next (cc, state, s[i]);
  }
  return (r);
}


/* ------------ helper functions  ------------------*/
static int encode_conv (it_convolutional_code_t * cc, int *state,	/* state of the coder */
			int input,	/* k bits */
			imat G,	/* generator matrix [k][n] */
			int Q,	/* feedback polynomial */
			int k,	/* number of input bits */
			int n)
{				/* number of output bits */
  int  i, j, b, r, t;
  int  c, s, l;
  int  m;
  int  feedback;

  s = *state;
  c = m = r = 0;

  for (j = 0; j < k; j++) {

    /* get the memory of this register */
    t = Q;
    for (i = 0; i < n; i++)
      t |= G[j][i];
    l = fls (t) - 1;

    /* compute feedback value */
    feedback = bit (input, j);
    for (i = 0; i < l; i++)
      feedback ^= bit (Q, i) & bit (s, i);

    /* generate redundancy bits */
    for (i = 0; i < n; i++) {
      c ^= (bit (G[j][i], l) & feedback) << i;
      for (b = 0; b < l; b++)
	c ^= (bit (G[j][i], b) & bit (s, b)) << i;
    }

    /* update coder state */
    r |= ((s & ((1 << l) - 1)) >> 1) << m;
    if (l)
      r |= feedback << (m + l - 1);
    s >>= l;
    m += l;
  }

  *state = r;

  return (c);
}


/* find last set bit */
static int fls (int c)
{
  int  l = 0;

  while (c) {
    l++;
    c >>= 1;
  }
  return (l);
}


ivec it_viterbi_decode_symbolic (it_convolutional_code_t * cc, mat metrics)
{
  ivec branches;
  imat next_state = cc->automaton;
  mat  next_state_logpt = cc->next_state_logpt;
  vec  logalpha_0 = cc->logalpha_0;
  vec  logbeta_end = cc->logbeta_end;
  int  j, l, k;

  /* decode the most probable sequence */
  branches =
    logviterbi (metrics, next_state, next_state_logpt, logalpha_0,
		logbeta_end);

  /* find the path leading to that sequence */
  j = 0;
  for (k = 0; k < ivec_length (branches); k++) {
    l = branches[k];
    branches[k] = cc->input[l][j];
    j = next_state[l][j];
  }

  return (branches);
}


bvec it_viterbi_decode (it_convolutional_code_t * cc, mat metrics)
{
  int  K, L;
  ivec branches;
  bvec r;
  int  i, j;

  K = cc->k;
  L = mat_width (metrics);
  branches = it_viterbi_decode_symbolic (cc, metrics);
  r = bvec_new (L * K);

  for (i = 0; i < L; i++)
    for (j = 0; j < K; j++)
      r[i * K + j] = bit (branches[i], K - 1 - j);

  ivec_delete (branches);

  return (r);
}
