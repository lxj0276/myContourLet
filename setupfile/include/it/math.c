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
  Functions
  Copyright (C) 2005 Vivien Chappelier, Hervé Jégou
*/

#include <it/types.h>
#include <it/math.h>
#include <it/vec.h>
#include <it/io.h>

/*----------------------------------------------------------------------------*/
/* the differentiation operator */
it_function (itf_differentiate)
{
  volatile double round;
  double epsilon;
  double fmin, fmax;
  it_function_t f = it_this->function;
  it_args_t args = it_this->args;

  /* some small constant (~1/3 of double precision) */
  epsilon = 1e-6 * ((x > 0) ? (x) : (-x));
  if (epsilon < 1e-6)
    epsilon = 1e-6;

  /* round off to appropriate binary representation */
  round = x + epsilon;
  epsilon = round - x;

  /* add some precision */
  fmax = f (x + epsilon, args);
  fmin = f (x - epsilon, args);

  return ((fmax - fmin) / (2.0 * epsilon));
}


/*----------------------------------------------------------------------------*/
/* the 2nd-order differentiation operator */
it_function (itf_diff2)
{
  volatile double round;
  double epsilon;
  double fmin, fmid, fmax;
  it_function_t f = it_this->function;
  it_args_t args = it_this->args;

  /* some small constant (~1/6 of double precision) */
  epsilon = 1e-3 * ((x > 0) ? (x) : (-x));
  if (epsilon < 1e-3)
    epsilon = 1e-3;

  /* round off to appropriate binary representation */
  round = x + epsilon;
  epsilon = round - x;

  /* add some precision */
  fmax = f (x + epsilon, args);
  fmin = f (x - epsilon, args);
  fmid = f (x, args);

  return ((fmax - 2 * fmid + fmin) / (epsilon * epsilon));
}


/*----------------------------------------------------------------------------*/
/* integration using the trapezoid method */
it_function (itf_integrate_trapezoid)
{
  double b = x;
  double a = it_this->a;
  int  N = it_this->N;
  it_function_t f = it_this->function;
  it_args_t args = it_this->args;
  double step;
  double sum;
  int  i;

  step = (b - a) / N;

  sum = f (a, args) / 2;
  x = a + step;

  for (i = 1; i <= N - 1; i++) {
    sum += f (x, args);
    x += step;
  }

  sum += f (b, args) / 2;
  sum *= step;

  return (sum);
}


/*----------------------------------------------------------------------------*/
/* integration using the Romberg method (Richardson + trapezoid) */
it_function (itf_integrate_romberg)
{
  double b = x;
  double a = it_this->a;
  int  N = it_this->N;
  int  p = it_this->p;
  it_function_t f = it_this->function;
  it_args_t args = it_this->args;
  double step;
  double r;
  int  i, n;
  vec  S;

  step = (b - a) / N;

  S = vec_new (p);

  /* compute the first trapezoid integral */
  S[0] = f (a, args) / 2;
  x = a + step;

  for (i = 1; i <= N - 1; i++) {
    S[0] += f (x, args);
    x += step;
  }

  S[0] += f (b, args) / 2;
  S[0] *= step;

  /* compute p trapezoid integrals with step/2^p recursively */
  for (n = 1; n < p; n++) {

    step /= 2;
    N *= 2;
    x = a + step;
    S[n] = 0;

    for (i = 1; i <= N - 1; i += 2) {
      S[n] += f (x, args);
      x += 2 * step;
    }
    S[n] *= step;
    S[n] += S[n - 1] / 2;
  }

  /* use the Richardson method to improve the approximation order */
  for (n = p - 1; n >= 0; n--) {
    r = (1 << 2 * (p - n));	/* 4^j */
    for (i = 0; i < n; i++)
      S[i] = (r * S[i + 1] - S[i]) / (r - 1);
  }

  r = S[0];

  vec_delete (S);

  return (r);
}


/*----------------------------------------------------------------------------*/
/* default integral (romberg with fixed precision) */
it_function (itf_integrate)
{
  it_function_args (itf_integrate_romberg) integrate_romberg_args;
  double b = x;
  double a = it_this->a;
  int  N;
  int  p;

  /* compute N so that step < 0.5 */
  N = (int) (2 * (b - a) + 1);
  p = 10;			/* sane default */

  integrate_romberg_args.a = a;
  integrate_romberg_args.N = N;
  integrate_romberg_args.p = p;
  integrate_romberg_args.function = it_this->function;
  integrate_romberg_args.args = it_this->args;

  return (itf_integrate_romberg (x, &integrate_romberg_args));
}


/*----------------------------------------------------------------------------*/
/* compute the expectation (first moment of the function) */
it_function (itf_expectation)
{
  it_function_args (itf_mul) mul_args;
  it_function_args (itf_integrate) integrate_args;

  /* for integrating x * f(x) */
  mul_args.f = itf_identity;
  mul_args.f_args = NULL;
  mul_args.g = it_this->function;
  mul_args.g_args = it_this->args;
  integrate_args.function = itf_mul;
  integrate_args.args = &mul_args;
  integrate_args.a = it_this->a;

  return (itf_integrate (x, &integrate_args));
}


/*----------------------------------------------------------------------------*/
/* function composition, returns f(g(x)) */
it_function (itf_compose)
{
  it_function_t f = it_this->f;
  it_function_t g = it_this->g;
  it_args_t f_args = it_this->f_args;
  it_args_t g_args = it_this->g_args;

  return (f (g (x, g_args), f_args));
}


/*----------------------------------------------------------------------------*/
/* function sum, returns f(x) + g(x) */
it_function (itf_sum)
{
  it_function_t f = it_this->f;
  it_function_t g = it_this->g;
  it_args_t f_args = it_this->f_args;
  it_args_t g_args = it_this->g_args;

  return (f (x, f_args) + g (x, g_args));
}


/*----------------------------------------------------------------------------*/
/* function product, returns f(x) * g(x) */
it_function (itf_mul)
{
  it_function_t f = it_this->f;
  it_function_t g = it_this->g;
  it_args_t f_args = it_this->f_args;
  it_args_t g_args = it_this->g_args;

  return (f (x, f_args) * g (x, g_args));
}


/*----------------------------------------------------------------------------*/
/* identity */
it_function (itf_identity)
{
  return (x);
}


/*----------------------------------------------------------------------------*/
/* Gaussian distribution */
it_function (itf_gaussian)
{
  double sigma = it_this->sigma;

  return (1. / (sqrt (2. * M_PI) * sigma) *
	  exp (-x * x / (2. * sigma * sigma)));
}


/*----------------------------------------------------------------------------*/
/* Laplacian distribution */
it_function (itf_laplacian)
{
  double lambda = it_this->lambda;	/* variance is 2*lambda^2 */

  return (1. / (2. * lambda) * exp (-fabs (x) / lambda));
}


/*----------------------------------------------------------------------------*/
/* Generalized Gaussian distribution */
#ifdef WIN32
it_function (itf_generalized_gaussian)
{
  it_fprintf (stderr, "undefined function lgamma()\n");
  return (NAN);
}
#else
it_function (itf_generalized_gaussian)
{
  double alpha = it_this->alpha;
  double beta = it_this->beta;

  assert (alpha > 0);
  assert (beta > 0);

  return (beta / (2.0 * alpha) *
	  exp (-lgamma (1.0 / beta) - pow (fabs (x) / alpha, beta)));
}
#endif

/*----------------------------------------------------------------------------*/
#define erfinv_a3 -0.140543331
#define erfinv_a2 0.914624893
#define erfinv_a1 -1.645349621
#define erfinv_a0 0.886226899

#define erfinv_b4 0.012229801
#define erfinv_b3 -0.329097515
#define erfinv_b2 1.442710462
#define erfinv_b1 -2.118377725
#define erfinv_b0 1

#define erfinv_c3 1.641345311
#define erfinv_c2 3.429567803
#define erfinv_c1 -1.62490649
#define erfinv_c0 -1.970840454

#define erfinv_d2 1.637067800
#define erfinv_d1 3.543889200
#define erfinv_d0 1

#ifdef WIN32
double erfinv (double x)
{
  it_fprintf (stderr, "undefined function erf()\n");
  return (NAN);
}
#else
double erfinv (double x)
{
  double x2, r, y;
  int  sign_x;

  if (x < -1 || x > 1)
    return NAN;

  if (x == 0)
    return 0;

  if (x > 0)
    sign_x = 1;
  else {
    sign_x = -1;
    x = -x;
  }

  if (x <= 0.7) {

    x2 = x * x;
    r =
      x * (((erfinv_a3 * x2 + erfinv_a2) * x2 + erfinv_a1) * x2 + erfinv_a0);
    r /= (((erfinv_b4 * x2 + erfinv_b3) * x2 + erfinv_b2) * x2 +
	  erfinv_b1) * x2 + erfinv_b0;
  }
  else {
    y = sqrt (-log ((1 - x) / 2));
    r = (((erfinv_c3 * y + erfinv_c2) * y + erfinv_c1) * y + erfinv_c0);
    r /= ((erfinv_d2 * y + erfinv_d1) * y + erfinv_d0);
  }

  r = r * sign_x;
  x = x * sign_x;

  r -= (erf (r) - x) / (2 / sqrt (M_PI) * exp (-r * r));
  r -= (erf (r) - x) / (2 / sqrt (M_PI) * exp (-r * r));

  return r;
}
#endif

#undef erfinv_a3
#undef erfinv_a2
#undef erfinv_a1
#undef erfinv_a0

#undef erfinv_b4
#undef erfinv_b3
#undef erfinv_b2
#undef erfinv_b1
#undef erfinv_b0

#undef erfinv_c3
#undef erfinv_c2
#undef erfinv_c1
#undef erfinv_c0

#undef erfinv_d2
#undef erfinv_d1
#undef erfinv_d0


/*----------------------------------------------------------------------------*/
static int nchoosek_tmp (int n, int k)
{
  if (k == 0)
    return 1;

  return (n * nchoosek_tmp (n - 1, k - 1)) / k;
}


int nchoosek (int n, int k)
{
  if (k > n || k < 0)
    return 0;

  if (k > n / 2)
    k = n - k;

  return nchoosek_tmp (n, k);
}


/*----------------------------------------------------------------------------*/
double lognchoosek (int n, int k)
{
  int  i;
  double r = 0;
  if (k > n || k < 0)
    return 0;

  if (k > n / 2)
    k = n - k;

  for (i = 1; i <= k; i++)
    r -= log (i);

  for (i = n - k + 1; i <= n; i++)
    r += log (i);

  return r;
}

/*----------------------------------------------------------------------------*/
double it_integrate (it_function_t function, it_args_t args, double a,
		     double b)
{
  it_function_args (itf_integrate) integrate_args;

  integrate_args.function = function;
  integrate_args.args = args;
  integrate_args.a = a;
  return (itf_integrate (b, &integrate_args));
}

double it_differentiate (it_function_t function, it_args_t args, double a)
{
  it_function_args (itf_differentiate) differentiate_args;

  differentiate_args.function = function;
  differentiate_args.args = args;
  return (itf_differentiate (a, &differentiate_args));
}



/*----------------------------------------------------------------------------*/
double log_sum(double log_a, double log_b)
{
  if (log_a < log_b)
    return log_b + log (1 + exp(log_a - log_b));
  else
    return log_a + log (1 + exp(log_b - log_a));
}


/* Compute the logarithm of the Gamma function */
double log_gamma (double x)
{
  double z=1/(x*x);

    x=x+6;
    z=(((-0.000595238095238*z+0.000793650793651)
	*z-0.002777777777778)*z+0.083333333333333)/x;
    z=(x-0.5)*log(x)-x+0.918938533204673+z-log(x-1)-
	log(x-2)-log(x-3)-log(x-4)-log(x-5)-log(x-6);
    return z;
}


/* Sigmoid and functions: 1 / (1 + exp(-lambda * x)) */
double sigmoid (double x, double lambda)
{
  return 1 / (1 + exp (-lambda * x));
}


double invsigmoid (double x, double lambda)
{
  return - log(1 / x - 1) / lambda;
}
