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
  Various Souce functions
  Copyright (C) 2005 Vivien Chappelier, Herve Jegou
*/


#ifndef __it_source_func_h
#define __it_source_func_h

#include <it/vec.h>
#include <it/mat.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Return the entropy of a memoryless discrete source */
double entropy( vec pdf );

/* Conditionnal entropy of a Markov chain defined by transition matrix pt */
double entropy_markov( mat pt );

/* Return the  entropy of a binary source source */
double entropy_bin( double p );

/* Return the histogram of the realization S. 
   The source is assumed to take its balues between 0 and omega-1 */
ivec histogram( int omega, ivec S );

/* Same as histogram return a normalized histogram, i.e. a pdf */
vec histogram_normalized( int omega, ivec S );

/* Return the conditionnal (i.e., bi-dimensional) histogram of a source realization */
imat histogram_cond( int omega, ivec S );

/* Return the expectation of a source defined by the pdf and the symbols values */
double source_expectation( vec pdf, vec symbols );

/* Return the variance of a source defined by its pdf and its symbols */
double source_variance( vec pdf, vec symbols );

/* Verify if the input vector is a valid probability density function */
int is_valid_pdf( vec pdf, double tol );

/* Check if pt is a valid matrix of transition probabilities (Markov chain) */
int is_valid_markov_matrix( mat pt, double tol );

/* Return the stationary probability of a Markov chain defined by transition matrix pt */
vec markov_marg_pdf( mat pt );


 

/*! \ingroup source
  \brief Convert a source on an alphabet A to a source on the product alphabet A^d.
         The size of the input vector should be a multiple of d.
  \param S This source must be a set of indexes between 0 and card(A)-1.
  \param cardA The cardinal of the base alphabet A
  \param d The dimension of the product alphabet considered here
*/
/*ivec to_product_alphabet( const ivec S, int cardA, int d );*/


/*! 
  \ingroup source
  \brief Convert a source on a product alphabet A^d to the source of indexes on A
  \param P This product source must be a set of indexes between 0 and card(A)^d-1.
  \param cardA The cardinal of the base alphabet A
  \param d The dimension of the product alphabet considered here
*/
/*ivec to_scalar_alphabet( const ivec P, int cardA, int d );*/


/*! \ingroup source
  \brief Return the stationnary transition probabilities of the product alphabet, 
  assuming a Markovian source. The order of the indexes is compatible with the one
  processed by functions scal_to_prod and prod_to_scal
  \param pt the transition probabilities of the source, assumed to be a first-order
  markov chain.
  \param d The power of the product space. 
  Ex: Let A be the alphabet of the source, of cardinal n. 
  The product alphabet is A^d, of cardinal n^d.
*/
/*vec alphabet_prod_pdf( const mat & pt, int d );*/

#ifdef __cplusplus
}
#endif /* extern "C" */
#endif

