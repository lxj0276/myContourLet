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
  Various source functions
  Copyright (C) 2005 Herve Jegou
*/


#include <math.h>
#include <it/source_func.h>
#include <it/math.h>
#include <it/io.h>

/*--------------------------------------------------------------------*/
/* Return the stationnary entropy of a symbol distribution pdf */
double entropy (vec pdf)
{
  idx_t i;
  double H = 0.0;

  for (i = 0; i < vec_length (pdf); i++)
    if (pdf[i] > 0)
      H -= pdf[i] * log (pdf[i]);

  return H / log (2.0);
}


/*--------------------------------------------------------------------*/
double entropy_bin (double p)
{
  if (p <= 0 || p > 1)
    return 0;
  else
    return -p * log2 (p) + (p - 1) * log2 (1 - p);
}

/*--------------------------------------------------------------------*/
double source_expectation (vec pdf, vec symbols)
{
  assert (vec_length (symbols) == vec_length (pdf));
  return vec_inner_product (symbols, pdf);
}



/*--------------------------------------------------------------------*/
double source_variance (vec pdf, vec symbols)
{
  double v;
  idx_t i;
  double e = source_expectation (symbols, pdf);
  assert (is_valid_pdf (pdf, 1e-6));

  v = 0;
  for (i = 0; i < vec_length (symbols); i++)
    v += pdf[i] * (symbols[i] - e) * (symbols[i] - e);

  return v;
}

/*-----------------------------------------------------------------------*/
double entropy_markov (mat pt)
{
  idx_t omega = mat_width (pt), s, i;
  double entropy = 0;

  /* Compute the stationnary probability law */
  vec  pdf = markov_marg_pdf (pt);
  it_assert (is_valid_markov_matrix (pt, IT_EPSILON),
	     "Matrix does not define a valid probability transition matrix");

  for (s = 0; s < omega; s++)
    for (i = 0; i < omega; i++)
      if (pt[i][s] > 0)
	entropy -= pdf[s] * pt[i][s] * log (pt[i][s]);

  vec_delete (pdf);
  return (entropy / log (2.0));
}



/*-----------------------------------------------------------------------
  !!!TODO: use a SVD instead ! (Single value decomposition)             */
vec markov_marg_pdf (mat pt)
{
  int  i;
  mat  Mpow = mat_clone (pt), Mtmp;
  vec  pdf;

  if (mat_height (pt) != mat_width (pt)) {
    it_error ("Numbers of columns and rows don't match, return void vector");
    return vec_null;
  }

  /* Process the power of the matrix */
  for (i = 0; i < 7; i++) {
    Mtmp = mat_new_mul (Mpow, Mpow);
    mat_delete (Mpow);
    Mpow = Mtmp;
  }

  pdf = mat_rows_sum (Mpow);
  mat_delete (Mpow);
  vec_normalize (pdf, 1);
  return (pdf);
}


/*-----------------------------------------------------------------------
 * Return the histogram of the realization SN
 * The source is assumed to be an index source that takes its value
 * between 0 and omega-1
 */
ivec histogram (int omega, ivec S)
{
  idx_t t;
  ivec histo = ivec_new_zeros (omega);

  for (t = 0; t < ivec_length (S); t++)
    if (S[t] < omega)
      histo[S[t]]++;

  return histo;
}


/*-----------------------------------------------------------------------
   Return the histogram of the realization SN normalized to sum to 1
   The source is assumed to be an index source that takes its value
   between 0 and omega-1
 */
vec histogram_normalized (int omega, ivec S)
{
  idx_t t;
  double delta = 1.0 / ivec_length (S);
  vec  histo = vec_new_zeros (omega);

  for (t = 0; t < ivec_length (S); t++)
    if (S[t] < omega)
      histo[S[t]] += delta;

  return histo;
}


/*-----------------------------------------------------------------------*/
imat histogram_cond (int omega, ivec S)
{
  imat histo2D = imat_new_zeros (omega, omega);
  idx_t t;

  for (t = 1; t < ivec_length (S); t++)
    histo2D[S[t]][S[t - 1]]++;

  return histo2D;
}


/*----------------------------------------------------------------------*/
int is_valid_pdf (vec pdf, double tol)
{
  assert (pdf);
  return fabs (1.0 - vec_sum (pdf)) < tol;
}


/*-----------------------------------------------------------------------
 Return true if the input matrix is a conditionnal probability density function, 
 false otherwise */
int is_valid_markov_matrix (mat pt, double tol)
{
  idx_t j;
  assert (pt);
  for (j = 0; j < mat_width (pt); j++)
    if (fabs (1.0 - mat_col_sum (pt, j)) < tol)
      return 0;
  return 1;
}

/*-----------------------------------------------------------------------
 * Convert a source on an alphabet A to a source on the product alphabet A^d
 */
/* ivec to_product_alphabet( ivec S, int cardA, int d ) */
/* { */
/*   int NS = S.length(); */
/*   it_assert( NS % d == 0, "Conversion to the product alphabet A^d requires the length of the input vector to be a multiple of d." ); */

/*   int NP = S.length() / d; */
/*   ivec P( NP ); */
/*   P.zeros(); */

/*   // is and ip are respectively point in S and P */
/*   // c is the symbol with a product symbol */
/*   int is, ip = -1, c; */
/*   for( is = 0 ; is < NS ; is++ ) */
/*     { */
/*       c = is % d; */
/*       if( c == 0 ) */
/* 	ip++; */

/*       P( ip ) = P( ip ) * cardA + S[ is ]; */
/*     } */

/*   return P; */
/* } */


/* //----------------------------------------------------------------------- */
/* // */
/* ivec to_scalar_alphabet( ivec P, int cardA, int d ) */
/* { */
/*   int NP = P.length(); */
/*   int NS = NP * d; */

/*   ivec S( NS ); */

/*   int is = 0, ip, c, cwd; */
/*   for( ip = 0 ; ip < NP ; ip++ ) */
/*     { */
/*       cwd = P[ ip ]; */
/*       for( c = d - 1 ; c >= 0 ; c-- ) */
/* 	{ */
/* 	  S[ is + c ] = cwd % cardA; */
/* 	  cwd /= cardA; */
/* 	} */
/*       is += d; */
/*     } */

/*   return S; */
/* } */

/* //----------------------------------------------------------------------- */
/* // Return the stationnary transition probabilities of the product alphabet,  */
/* //  assuming a Markovian source */
/* vec alphabet_prod_pdf( mat pt, int d ) */
/* { */
/*   it_assert( pt.rows() == pt.cols(), "Matrix of probability is not square" ); */

/*   int n = pt.rows(); */
/*   int prod_n = (int) round( pow( n, d ) ); */

/*   vec prod_pdf( prod_n ); */
/*   prod_pdf.ones(); */
/*   vec pdf = marg_pdf( pt );   // The marginal law for the scalar alphabet */

/*   // cwd denotes a codeword of the product space  */
/*   int cwd_tmp, s, last_s, i, f; */
/*   double p; */
/*   for( int cwd = 0 ; cwd < prod_n ; cwd++ ) */
/*     { */
/*       last_s = -1; */
/*       cwd_tmp = cwd; */
/*       f = prod_n / n; */

/*       // i indexes the dimension */
/*       for( i = 0 ; i < d ; i++ ) */
/* 	{ */
/* 	  s = cwd_tmp / f; */
/* 	  cwd_tmp %= f; */
/* 	  f /= n; */

/* 	  // First symbol : no knowledge of past realizations */
/* 	  if( last_s == -1 ) */
/* 	    p = pdf( s ); */

/* 	  // Use previous symbol for the context */
/* 	  else */
/* 	    p = pt( s, last_s ); */

/* 	  prod_pdf( cwd ) *= p; */
/* 	  last_s = s; */
/* 	} */

/*     } */
/*   return prod_pdf; */
/* } */
