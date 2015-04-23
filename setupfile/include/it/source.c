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
  Source generators
  Copyright (C) 2005 Vivien Chappelier, Herve Jegou
*/


#include <it/random.h>
#include <it/source.h>
#include <it/source_func.h>
#include <it/io.h>
#include <it/math.h>

/* generate a vector of random binary values with the given */
/* probability of the 0 symbol */
bvec source_binary (idx_t size, double p0)
{
  bvec v = bvec_new (size);
  idx_t i;

  for (i = 0; i < size; i++)
    v[i] = (it_rand () >= p0);

  return (v);
}


/* generate a vector of values uniformly distributed in [a,b[ */
vec source_uniform (idx_t size, double a, double b)
{
  vec  v = vec_new (size);
  idx_t i;

  for (i = 0; i < size; i++)
    v[i] = a + (b - a) * it_rand ();

  return (v);
}

/* generate a vector of values uniformly distributed in [a,b[ */
ivec source_uniform_int (idx_t size, int a, int b)
{
  ivec  v = ivec_new (size);
  idx_t i;

  for (i = 0; i < size; i++)
    v[i] = floor ((double)a + (b - a) * it_rand ());

  return (v);
}

/* generate a vector of independent values drawn from a */
/* gaussian distribution of given mean and standard deviation */
vec source_gaussian (idx_t size, double mean, double std)
{
  vec  v = vec_new (size);
  idx_t i;

  for (i = 0; i < size; i++)
    v[i] = mean + std * it_randn ();

  return (v);
}


/* generate a stationnary random vector from a probability
   density function using the acceptance-rejection method.
   the pdf is assumed to be zero outside [a, b].
*/
vec source_pdf (idx_t size, double a, double b, it_function_t pdf,
		it_args_t args)
{
  idx_t i;
  vec  v = vec_new (size);

  for (i = 0; i < size; i++)
    v[i] = it_randpdf (a, b, pdf, args);

  return (v);
}


/*--------------------------------------------------------------------*/
ivec source_memoryless (idx_t size, vec pdf)
{
  idx_t t;
  double x;
  int  s, e, h;
  vec  pdfcs;
  ivec S = ivec_new (size);

  it_assert (is_valid_pdf (pdf, IT_EPSILON), "Invalid probability law");
  pdfcs = vec_cum_sum (pdf);
  vec_ins (pdfcs, 0, 0.0);

  for (t = 0; t < size; t++) {
    /* generate a uniformly distributed value  */
    /* and search the integral of the pdf for  */
    /* the abciss corresponding to that random */
    /* value.                                  */
    x = it_rand ();
    s = 0;
    e = vec_length (pdfcs) - 1;

    while (e > s + 1) {
      /* compare with the median of the range */
      /* and narrow down the range. */
      h = (s + e) / 2;

      if (x < pdfcs[h])
	e = h;
      else
	s = h;
    }

    if (x < pdfcs[e])
      S[t] = s;
    else
      S[t] = e;
  }
  vec_delete (pdfcs);
  return S;
}
