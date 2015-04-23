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
  Source fitting functions
  Copyright (C) 2005 Vivien Chappelier, Herve Jegou
*/


#include <math.h>
#include <it/source_fit.h>
#include <it/distance.h>
#include <it/math.h>
#include <it/io.h>

/*--------------------------------------------------------------------*/
static double GenGaussian (double alpha, double beta, double x)
{
  return exp (-pow (fabs (x) / alpha, beta));
}


/*--------------------------------------------------------------------*/
/* Return the discrete pdf law associated with the Generalized 
   Gaussian of parameters gamma, beta and C                           */
vec source_pdf_GG (vec symbols, double alpha, double beta)
{
  idx_t i;
  vec  pdf = vec_new_zeros (vec_length (symbols));
  for (i = 0; i < vec_length (symbols); i++)
    pdf[i] = GenGaussian (alpha, beta, symbols[i]);

  vec_normalize (pdf, 1);
  return pdf;
}



/*--------------------------------------------------------------------*/
void source_estim_GG_from_histo (vec pdf, vec symbols, double *palpha,
				 double *pbeta)
{
  double step_alpha, step_beta;
  double dK, last_dK;
  int  i;

  /* derivatives of the Kullback distance on alpha and beta */
  double derdK_alpha, derdK_beta;
  vec  pdf_GG, pdf_GGbisalpha, pdf_GGbisbeta;

  /* Kind of a gradient descent with the derivative of kullback distance.   */
  step_alpha = *palpha / 10;
  step_beta = *pbeta / 10;
  last_dK = 666;

  /* Initial Kullback distance */
  pdf_GG = source_pdf_GG (symbols, *palpha, *pbeta);
  dK = vec_distance_kullback_leibler (pdf, pdf_GG);
  vec_delete (pdf_GG);

  /* "Geometric" decreasing of the convergence step */
  for (i = 0; i < 5; i++) {
    do {
      pdf_GGbisalpha = source_pdf_GG (symbols, *palpha + IT_EPSILON, *pbeta);
      derdK_alpha =
	(vec_distance_kullback_leibler (pdf, pdf_GGbisalpha) -
	 dK) / IT_EPSILON;
      vec_delete (pdf_GGbisalpha);

      pdf_GGbisbeta = source_pdf_GG (symbols, *palpha, *pbeta + IT_EPSILON);
      derdK_beta =
	(vec_distance_kullback_leibler (pdf, pdf_GGbisbeta) -
	 dK) / IT_EPSILON;
      vec_delete (pdf_GGbisbeta);

      *palpha = *palpha - derdK_alpha * step_alpha;
      *pbeta = *pbeta - derdK_beta * step_beta;

      last_dK = dK;

      pdf_GG = source_pdf_GG (symbols, *palpha, *pbeta);
      dK = vec_distance_kullback_leibler (pdf, pdf_GG);
      vec_delete (pdf_GG);
      //              printf( "alpha = %f, beta = %f, dK = %f\n", *palpha, *pbeta, dK );
    }
    while (dK + IT_EPSILON * 1e3 < last_dK);

    step_alpha *= 0.5;
    step_beta *= 0.5;
  }
}
