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


#ifndef __it_source_fit_h
#define __it_source_fit_h

#include <it/vec.h>

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------------------*
 * Fit a given histogram realization with a model                     *
 *--------------------------------------------------------------------*/


/* Return the discrete pdf law associated with the Generalized 
   Gaussian of parameters gamma, beta and C                           */
vec source_pdf_GG( vec symbols, double alpha, double beta );

/* Estimate the parameters of a generalized Gaussian in the Kullback distance sense */
void source_estim_GG_from_histo( vec pdf, vec symbols, double * alpha, double * beta );

#ifdef __cplusplus
}
#endif /* extern "C" */
#endif
