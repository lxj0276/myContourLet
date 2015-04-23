/*
   contourlet - Implementation of the contourlet transform for image coding
   Copyright (C) 2005 Vivien Chappelier - IRISA/University of Rennes 1

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public
   License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU General Public
   License along with this program; if not, write to the Free
   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#ifndef __PYRAMID_H
#define __PYRAMID_H

#include "it/mat.h"

/* split an image into two laplacian pyramid bands */
void laplacian_pyramid_split(mat image, /* input band                 */
			     mat low,   /* low-frequency band         */
			     mat high,  /* high-frequency band        */
			     vec H0,    /* symmetric analysis  filter */
			     vec G0);   /* symmetric synthesis filter */


/* merge two laplacian pyramid bands into an image */
void laplacian_pyramid_merge(mat image, /* output band                */
			     mat low,   /* low-frequency band         */
			     mat high,  /* high-frequency band        */
			     vec H0,    /* symmetric analysis  filter */
			     vec G0);   /* symmetric synthesis filter */

#endif
