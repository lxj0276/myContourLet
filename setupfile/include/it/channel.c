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
  Channels and modulation
  Copyright (C) 2005 Vivien Chappelier, Herve Jegou
*/

#include <it/vec.h>
#include <it/source.h>
#include <it/random.h>
#include <it/channel.h>

/* simple BSPK modulation */
vec modulate_bpsk (bvec b)
{
  idx_t i;
  vec  m = vec_new (bvec_length (b));

  for (i = 0; i < bvec_length (b); i++)
    m[i] = b[i] ? 1 : -1;

  return (m);
}


bvec channel_bsc (bvec v, double crossover_proba)
{
  idx_t i;
  bvec received = bvec_new (bvec_length (v));

  for (i = 0; i < bvec_length (v); i++)
    if (it_rand () < crossover_proba)
      received[i] = 1 - v[i];
    else
      received[i] = v[i];

  return received;
}


vec channel_awgn (vec v, double sigma)
{
  vec  received = source_gaussian (vec_length (v), 0.0, sigma);
  vec_add (received, v);
  return received;
}
