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
  1D wavelet transform using lifting.
  Copyright (C) 2005 Vivien Chappelier, Herve Jegou
*/

#include <it/types.h>
#include <it/wavelet.h>
#include <it/io.h>

/* This computes the wavelet decomposition using the lifting implementation.
*/

/* some predefined filters */
static const double __it_wavelet_lifting_97_steps[] =
  { -1.586134342, -0.052980118, 0.882911075, 0.443506852 };
static const it_wavelet_lifting_t __it_wavelet_lifting_97 = {
  4,
  1.149604398,
  __it_wavelet_lifting_97_steps
};
it_wavelet_lifting_t const *it_wavelet_lifting_97 = &__it_wavelet_lifting_97;

static const double __it_wavelet_lifting_53_steps[] = { -0.5, 0.25 };
static const it_wavelet_lifting_t __it_wavelet_lifting_53 = {
  2,
  1.41421,
  __it_wavelet_lifting_53_steps
};
it_wavelet_lifting_t const *it_wavelet_lifting_53 = &__it_wavelet_lifting_53;

/*
  X[0] = X[0] + alpha * (X[-1] + X[1]); odd
  X[0] = X[0] + beta  * (X[-1] + X[1]);  even
  X[0] = X[0] + gamma * (X[-1] + X[1]); odd
  X[0] = X[0] + delta * (X[-1] + X[1]);  even
  X[0] = -scale * X[0];                 odd
  X[0] = X[0] / scale;                   even
*/

#define shift_up(x, l) (((x) + ~(-1 << (l))) >> (l))

/* compute the next level decomposition using the lifting method. */
static int __wavelet_split (it_wavelet_t * wavelet)
{
  int  i, j, l, h;
  int  length;
  int  level;
  int  levels;
  int  count;
  double s;
  double const *step;
  double scale;
  vec  v;
  vec  r;

  assert (wavelet);

  if (wavelet->level == wavelet->levels)
    return (-IT_EINVAL);

  v = wavelet->current;
  r = wavelet->next;
  length = wavelet->length;
  level = wavelet->level;
  levels = wavelet->levels;
  count = wavelet->lifting->count / 2;
  step = wavelet->lifting->step;
  scale = wavelet->lifting->scale;

  l = shift_up (length, level);
  h = shift_up (length, level + 1);

  /* lifting */
  for (j = 0; j < count; j++) {
    /* stage 1 : odd samples lifting */
    s = step[2 * j];
    if (l & 1) {
      for (i = 1; i < l - 1; i += 2)
	v[i] += s * (v[i - 1] + v[i + 1]);
    }
    else {
      for (i = 1; i < l - 1; i += 2)
	v[i] += s * (v[i - 1] + v[i + 1]);
      v[i] += s * 2 * v[i - 1];
    }

    /* stage 2 : even samples lifting */
    s = step[2 * j + 1];
    if (l & 1) {
      v[0] += s * 2 * v[1];
      for (i = 2; i < l - 1; i += 2)
	v[i] += s * (v[i - 1] + v[i + 1]);
      if (l > 2)
	v[i] += s * 2 * v[i - 1];
    }
    else {
      v[0] += s * 2 * v[1];
      for (i = 2; i < l - 1; i += 2)
	v[i] += s * (v[i - 1] + v[i + 1]);
    }
  }

  /* scale and demultiplex */
  for (i = 0; i < l - 1; i += 2) {
    r[(i >> 1)] = v[i] * scale;
    r[h + (i >> 1)] = v[i + 1] / scale;
  }
  if (l & 1)
    r[(i >> 1)] = v[i] * scale;

  /* swap the vectors */
  wavelet->current = r;
  wavelet->next = v;

  wavelet->level++;

  return (0);
}

/* reconstruct the previous level of decomposition
   by inverting the lifting steps. */
static int __wavelet_merge (it_wavelet_t * wavelet)
{
  int  i, j, l, h;
  int  length;
  int  level;
  int  levels;
  int  count;
  double const *step;
  double scale;
  double s;
  vec  v;
  vec  r;

  assert (wavelet);

  if (wavelet->level == 0)
    return (-IT_EINVAL);

  wavelet->level--;

  r = wavelet->current;
  v = wavelet->next;
  length = wavelet->length;
  level = wavelet->level;
  levels = wavelet->levels;
  count = wavelet->lifting->count / 2;
  step = wavelet->lifting->step;
  scale = wavelet->lifting->scale;

  l = shift_up (length, level);
  h = shift_up (length, level + 1);

  /* scale and multiplex */
  for (i = 0; i < l - 1; i += 2) {
    v[i] = r[(i >> 1)] / scale;
    v[i + 1] = r[h + (i >> 1)] * scale;
  }
  if (l & 1)
    v[i] = r[(i >> 1)] / scale;

  /* invert lifting */
  for (j = count - 1; j >= 0; j--) {

    /* stage 1 : even samples lifting */
    s = step[2 * j + 1];
    if (l & 1) {
      v[0] -= s * 2 * v[1];
      for (i = 2; i < l - 1; i += 2)
	v[i] -= s * (v[i - 1] + v[i + 1]);
      if (l > 2)
	v[i] -= s * 2 * v[i - 1];
    }
    else {
      v[0] -= s * 2 * v[1];
      for (i = 2; i < l - 1; i += 2)
	v[i] -= s * (v[i - 1] + v[i + 1]);
    }

    /* stage 2 : odd samples lifting */
    s = step[2 * j];
    if (l & 1) {
      for (i = 1; i < l - 1; i += 2)
	v[i] -= s * (v[i - 1] + v[i + 1]);
    }
    else {
      for (i = 1; i < l - 1; i += 2)
	v[i] -= s * (v[i - 1] + v[i + 1]);
      v[i] -= s * 2 * v[i - 1];
    }
  }

  /* swap the vectors */
  wavelet->current = v;
  wavelet->next = r;

  return (0);
}

static void __wavelet_pack (vec flat, vec current, vec next, int levels)
{
  int  length;
  int  l;
  int  i;
  vec  v;

  /* merge the double buffer into a flat vector                       */
  /* with the lowest band stored in the 'current' buffer              */
  /* For example:                                                     */
  /* current: [ LLL | LLH |     ?     |           HHH           ]     */
  /* next   : [  ?  |  ?  |    LHH    |            ?            ]     */
  /* becomes                                                          */
  /* flat:    [ LLL | LLH |    LHH    |           HHH           ]     */

  length = vec_length (flat);

  /* low pass band */
  v = current;
  for (i = 0; i < shift_up (length, levels); i++)
    flat[i] = v[i];

  /* high pass bands */
  for (l = levels; l > 0; l--) {
    v = ((levels - l) & 1) ? next : current;
    for (; i < shift_up (length, l - 1); i++)
      flat[i] = v[i];
  }
}

static void __wavelet_unpack (vec flat, vec current, vec next, int levels)
{
  int  length;
  int  l;
  int  i;
  vec  v;

  /* split the flat vector into the double buffer,                    */
  /* with the lowest band stored in the 'current' buffer              */
  /* For example:                                                     */
  /* flat:    [ LLL | LLH |    LHH    |           HHH           ]     */
  /* becomes                                                          */
  /* current: [ LLL | LLH |     ?     |           HHH           ]     */
  /* next   : [  ?  |  ?  |    LHH    |            ?            ]     */
  /* why bother and not copy flat into both buffers? .. efficiency :) */

  length = vec_length (flat);

  /* low pass band */
  v = current;
  for (i = 0; i < shift_up (length, levels); i++)
    v[i] = flat[i];

  /* high pass bands */
  for (l = levels; l > 0; l--) {
    v = ((levels - l) & 1) ? next : current;
    for (; i < shift_up (length, l - 1); i++)
      v[i] = flat[i];
  }
}

/* compute the wavelet transform of the image */
static Vec __wavelet_transform (it_transform_t * transform, Vec __input)
{
  it_wavelet_t *wavelet = IT_WAVELET (transform);
  vec  flat;
  idx_t length;
  vec  input = (vec) __input;
  int  free_on_exit = 0;

  assert (input);
  /* check the input is actually a vec */
  assert (Vec_header (input).element_size == sizeof (double));

  length = vec_length (input);

  if (!wavelet->length) {
    it_transform_set_size (wavelet, length);
    free_on_exit = 1;
  }

  assert (wavelet->length == length);

  wavelet->level = 0;

  flat = vec_new (wavelet->length);

  vec_copy (wavelet->current, input);

  while (wavelet->level < wavelet->levels)
    __wavelet_split (wavelet);

  __wavelet_pack (flat, wavelet->current, wavelet->next, wavelet->levels);

  if (free_on_exit)
    it_transform_clear_size (wavelet);

  return (flat);
}

/* compute the inverse wavelet transform of the coefficients */
static Vec __wavelet_itransform (it_transform_t * transform, Vec __flat)
{
  it_wavelet_t *wavelet = IT_WAVELET (transform);
  vec  output;
  idx_t length;
  vec  flat = (vec) __flat;
  int  free_on_exit = 0;

  assert (flat);

  /* check the input is actually a vec */
  assert (Vec_header (flat).element_size == sizeof (double));

  length = vec_length (flat);

  if (!wavelet->length) {
    it_transform_set_size (wavelet, length);
    free_on_exit = 1;
  }

  wavelet->length = length;
  wavelet->level = wavelet->levels;

  __wavelet_unpack (flat, wavelet->current, wavelet->next, wavelet->levels);

  while (wavelet->level)
    __wavelet_merge (wavelet);

  output = vec_clone (wavelet->current);

  if (free_on_exit)
    it_transform_clear_size (wavelet);

  return (output);
}

static void __wavelet_get_output_size (it_transform_t * transform,
				       idx_t * input_size)
{
  /* this transform is critically sampled */
}

static void __wavelet_set_size (it_transform_t * transform, idx_t length)
{
  it_wavelet_t *wavelet = IT_WAVELET (transform);

  if (wavelet->length) {
    vec_delete (wavelet->current);
    vec_delete (wavelet->next);
  }

  if (length) {
    wavelet->current = vec_new (length);
    wavelet->next = vec_new (length);
  }

  wavelet->length = length;
}

static void __wavelet_get_size (it_transform_t * transform, idx_t * length)
{
  it_wavelet_t *wavelet = IT_WAVELET (transform);
  *length = wavelet->length;
}

static void __wavelet_destructor (it_object_t * it_this)
{
  it_wavelet_t *wavelet = IT_WAVELET (it_this);

  if (wavelet->length) {
    vec_delete (wavelet->current);
    vec_delete (wavelet->next);
  }

  /* call the parent destructor */
  wavelet->it_overloaded (destructor) (it_this);
}

it_instanciate (it_wavelet_t)
{
  it_new_args_start ();
  it_construct (it_transform_t);
  it_set_magic (it_this, it_wavelet_t);

  /* overload the virtual destructor */
  it_overload (it_this, it_object_t, destructor, __wavelet_destructor);

  IT_TRANSFORM (it_this)->transform = __wavelet_transform;
  IT_TRANSFORM (it_this)->itransform = __wavelet_itransform;
  IT_TRANSFORM (it_this)->get_output_size = __wavelet_get_output_size;
  IT_TRANSFORM (it_this)->set_size = __wavelet_set_size;
  IT_TRANSFORM (it_this)->get_size = __wavelet_get_size;

  it_this->level = 0;
  it_this->lifting = it_new_args_next (it_wavelet_lifting_t *);
  it_this->levels = it_new_args_next (int);
  it_this->length = 0;

  it_new_args_stop ();

  return (it_this);
}

/*--------------------------------------------------------------------*/
vec *it_wavelet_split (vec wav, int nb_levels)
{
  int  mid, nb, l;
  vec *subbands;

  assert (wav);

  nb = vec_length (wav);
  mid = (nb + 1) / 2;

  subbands = (vec *) malloc (sizeof (vec) * (nb_levels + 1));

  for (l = nb_levels; l > 0; l--) {
    subbands[l] = vec_get_subvector (wav, mid, nb - 1);
    nb = mid;
    mid = (nb + 1) / 2;
  }

  subbands[0] = vec_get_subvector (wav, 0, nb - 1);
  return subbands;
}

vec it_wavelet_merge (vec * subbands, int nb_levels)
{
  int  mid, nb, l;
  vec  wav;
  assert (subbands);

  nb = 0;
  for (l = 0; l <= nb_levels; l++)
    nb += vec_length (subbands[l]);

  wav = vec_new (nb);
  mid = (nb + 1) / 2;

  for (l = nb_levels; l > 0; l--) {
    vec_set_subvector (wav, subbands[l], mid);
    nb = mid;
    mid = (nb + 1) / 2;
  }

  vec_set_subvector (wav, subbands[0], 0);
  return wav;
}

vec it_dwt (vec v, it_wavelet_lifting_t const *lifting, int levels)
{
  it_wavelet_t *wavelet;
  vec  t;

  wavelet = it_wavelet_new (lifting, levels);

  t = (vec) it_transform (wavelet, v);

  it_delete (wavelet);

  return (t);
}

vec it_idwt (vec t, it_wavelet_lifting_t const *lifting, int levels)
{
  it_wavelet_t *wavelet;
  vec  v;

  wavelet = it_wavelet_new (lifting, levels);

  v = (vec) it_itransform (wavelet, t);

  it_delete (wavelet);

  return (v);
}
