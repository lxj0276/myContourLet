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
  Copyright (C) 2005 Vivien Chappelier, Herve Jegou
*/

#ifndef __it_distance_h
#define __it_distance_h

#include "it/types.h"
#include "it/vec.h"
#include "it/mat.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Return the Hamming distance between two vectors. 
   The Hamming distance is the number of elements that are not equal. If the vectors are 
   not of the same size, the distance is increased by the difference of lengths.          */
int vec_distance_hamming( vec v1, vec v2 );
int ivec_distance_hamming( ivec v1, ivec v2 );
int bvec_distance_hamming( bvec v1, bvec v2 );

/* Return the symbol error rate or bit error rate corresponding to the distance 
   between the input vector and the output vector. If v2 has a smaller size than 
   v1, then all the symbols which are not received are assumed to be in error. 
   If v2 has a greater size, the symbols are disgarded (but not added to increase the ser.*/
double vec_ser( vec v1, vec v2 );
double ivec_ser( ivec v1, ivec v2 );
double bvec_ber( bvec v1, bvec v2 );

/* Return the levenshtein distance between two vectors.
   The distance is computed with the algorithm of Wagner and Fischer.
   The cost of insertion, deletion and substitution operations are parameters. 
   Usually, they are all assigned 1.                                                      */
int ivec_distance_levenshtein( ivec v1, ivec v2, int cost_ins, int cost_del, int cost_sub );

/* Return the distance derived from a norm between two vectors. 
   This distance is given by the norm of the difference between the two vectors.
   If the vectors are not of the same size, the distance is increased assuming 
   that missing elements are equal to 0.
   Note that norm is often assigned 2 (square error)
*/
double vec_distance_norm( vec v1, vec v2, double norm );
double mat_distance_norm( mat v1, mat v2, double norm );


/* Return the Mean Square Error (MSE) between vector v2 and vector v1.
   Missing elements are assumed to be equal to the parameter 'rec_value'. */
double vec_distance_mse( vec v1, vec v2, double rec_value );
double mat_distance_mse( mat v1, mat v2, double rec_value );

double ivec_distance_mse( ivec v1, ivec v2, double rec_value );
double imat_distance_mse( imat v1, imat v2, double rec_value );

/* Return the Kullback-Leibler pseudo-distance between distribution pdf1 and pdf2.
   In other terms, the quantity \sum_{a in A} pdf1(a) log2(pdf1(a)/log(pdf1(a)/pdf2(a))) 
   is computed                                                                            */
double vec_distance_kullback_leibler( vec pdf1, vec pdf2 );

/* *** The following distance has to be implemented:
   - the symmetrized kullback distance
*/

#ifdef __cplusplus
}
#endif /* extern "C" */
#endif
