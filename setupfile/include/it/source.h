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
  Source definition functions
  Copyright (C) 2005 Vivien Chappelier, Herve Jegou
*/


#ifndef __it_source_h
#define __it_source_h

#include <it/vec.h>

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------------------*
 * Source Random number generator                                     *
 *--------------------------------------------------------------------*/

/* generate a vector of random binary values with the given */
/* probability for the 0 symbol */
bvec source_binary( size_t size, double p0 );

/* generate a vector of values uniformly distributed in [a,b[ */
vec source_uniform( size_t size, double a, double b );

/* generate a vector of independent values drawn from a */
/* gaussian distribution of given mean and standard deviation */
vec source_gaussian( size_t size, double mean, double std );

/* generate a stationnary random vector from a probability
   density function using the acceptance-rejection method.
   the pdf is assumed to be zero outside [a, b].
*/
vec source_pdf( size_t size, double a, double b, it_function_t pdf, it_args_t args );

/* Return a vector of size K of real values between 0 and 1           */
#define source_uniform_01( K ) source_uniform( K, 0, 1 )

/* Memoryless discrete source of length size defined by its stationary probabilities pdf */
ivec source_memoryless( size_t size, vec pdf );

#ifdef __cplusplus
}
#endif /* extern "C" */
#endif
