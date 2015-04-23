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
  Mathematical functions
  Copyright (C) 2005 Vivien Chappelier
*/

#ifndef __it_math_h
#define __it_math_h

#include <math.h>
#include <it/types.h>
#include <it/math.h>

#ifdef __cplusplus
extern "C" {
#endif


#ifndef M_PI
#define M_PI 3.141592653589793238
#endif

#ifndef M_E
#define M_E 2.71828182845905
#endif

/* logarithm in base 2 */
#define LOG2    0.69314718055994529
#define INVLOG2 1.44269504088896339
#define log2(x) (INVLOG2 * log(x))

/* some small constant */
#define IT_EPSILON (1e-10)

/*-----------------------------------------------------------------*/
/*                 Basic mathematical functions                    */
/*-----------------------------------------------------------------*/

/*---------------- some well-known functions ----------------------*/

/* identity, returns x */
it_function_t const itf_identity;

/*------------------------ Operators -----------------------------*/
/* These are just particular functions with one or more function(s)
   (and its parameters) passed as parameters. */

/* the differentiation operator */
it_function_args(itf_differentiate) {
  it_function_t function; /* function to differentiate */
  it_args_t args;         /* and its arguments */
};
it_function_t const itf_differentiate;

/* integration between [a, b] using the trapezoid method (N samples).
   the upper bound is the function argument whereas the lower bound
   and number of intervals are passed as parameters (along with the
   function to integrate and its arguments). This allows for a 
   consistent definition of the integral that can be differentiated
   again to obtain [an approximation of] the original function. */
it_function_args(itf_integrate_trapezoid) {
  double a;               /* lower bound */
  int N;                  /* number of intervals */
  it_function_t function; /* function to integrate */
  it_args_t args;         /* and its arguments */
};
it_function_t const itf_integrate_trapezoid;

/* integration between [a, b] using Romberg convergence acceleration method. */
/* the error of this method is O(((b-a)/N)^(2p+2))                           */
it_function_args(itf_integrate_romberg) {
  double a;               /* lower bound */
  int N;                  /* number of intervals */
  int p;                  /* precision order */
  it_function_t function; /* function to integrate */
  it_args_t args;         /* and its arguments */
};
it_function_t const itf_integrate_romberg;

/* default integration between [a, b]. This uses the Romberg method with
   a fixed precision of O(2^-12). The constant depends on the value
   of the tenth and further derivatives of the function, which are usually
   quite small for smooth enough functions.
 */
it_function_args(itf_integrate) {
  double a;               /* lower bound */
  it_function_t function; /* function to integrate */
  it_args_t args;         /* and its arguments */
};
it_function_t const itf_integrate;

/* compute the expectation (first moment of the function) */
it_function_args(itf_expectation)
{
  double a;               /* lower bound */
  it_function_t function; /* function to integrate */
  it_args_t args;         /* and its arguments */
};
it_function_t const itf_expectation;

/*------------------------- Binary operators ----------------------*/

/* function composition, returns f(g(x)) */
it_function_args(itf_compose)
{
  it_function_t f;
  it_args_t f_args;
  it_function_t g;
  it_args_t g_args;
};
it_function_t const itf_compose;

/* function sum, returns f(x) + g(x) */
it_function_args(itf_sum)
{
  it_function_t f;
  it_args_t f_args;
  it_function_t g;
  it_args_t g_args;
};
it_function_t const itf_sum;

/* function product, returns f(x) * g(x) */
it_function_args(itf_mul)
{
  it_function_t f;
  it_args_t f_args;
  it_function_t g;
  it_args_t g_args;
};
it_function_t const itf_mul;

#ifdef __cplusplus
}
#endif /* extern "C" */
#endif
