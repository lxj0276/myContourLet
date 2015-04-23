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

#ifndef __it_random_h
#define __it_random_h

#include <it/vec.h>

#ifdef __cplusplus
extern "C" {
#endif

/* initialize the random number generator (with a random seed)    */
/* Note: the seed is taken from the milliseconds of the current   */
/* time, which is not a serious option for security applications. */
/* In this case, always use it_seed with your favorite method     */
/* to obtain a good seed.                                         */
void it_randomize(void);

/* intializes the random generator from a seed */
void it_seed(unsigned int seed);

/* generate a value uniformly distributed in [0,1[ */
double it_rand(void);

/* generate a value uniformly distributed normally */
double it_randn(void);

/* generate a random variable from its probability
   density function using the acceptance-rejection method.
   the pdf is assumed to be zero outside [a, b].
*/
double it_randpdf(double a, double b, it_function_t pdf, it_args_t args);


#ifdef __cplusplus
}
#endif /* extern "C" */
#endif
