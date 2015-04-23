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
  Separable 2D transform from a 1D transform
  Copyright (C) 2005 Vivien Chappelier
*/

#include <it/types.h>
#include <it/transform.h>
#include <it/transform2D.h>
#include <it/separable2D.h>
#include <it/io.h>

/* compute the transform of the image */
static Mat __it_separable2D_transform (it_transform2D_t * transform2D,
				       Mat image)
{
  it_separable2D_t *separable = IT_SEPARABLE2D (transform2D);
  it_transform_t *transform;
  int  x, y;
  int  width, height;
  Mat  tmp, output;
  Vec  v;
  idx_t elem_size;

  assert (image);
  assert (image[0]);

  width = Mat_width (image);
  height = Mat_height (image);
  transform = separable->transform;
  elem_size = Vec_element_size (image[0]);

  /* create a temporary matrix to store the transposed coefficients */
  tmp = __Mat_new (elem_size, width, height);

  /* apply on all rows */
  it_transform_set_size (transform, width);
  for (y = 0; y < height; y++) {
    v = it_transform (transform, image[y]);
    /* store in columns */
    Mat_set_col (tmp, y, v);
    Vec_delete (v);
  }

  /* create the output matrix */
  output = __Mat_new (elem_size, height, width);

  /* apply on all columns (stored in rows of 'tmp') */
  it_transform_set_size (transform, height);
  for (x = 0; x < width; x++) {
    v = it_transform (transform, tmp[x]);
    /* store in columns */
    Mat_set_col (output, x, v);
    Vec_delete (v);
  }
  Mat_delete (tmp);

  return (output);
}

static Mat __it_separable2D_itransform (it_transform2D_t * transform2D,
					Mat input)
{
  it_separable2D_t *separable = IT_SEPARABLE2D (transform2D);
  it_transform_t *transform;
  int  x, y;
  int  width, height;
  Mat  tmp, output;
  Vec  v;
  idx_t elem_size;

  assert (input);
  assert (input[0]);

  width = Mat_width (input);
  height = Mat_height (input);
  transform = separable->transform;
  elem_size = Vec_element_size (input[0]);

  /* create a temporary matrix to store the transposed coefficients */
  tmp = __Mat_new (elem_size, width, height);

  /* apply on all rows */
  it_transform_set_size (transform, width);
  for (y = 0; y < height; y++) {
    v = it_itransform (transform, input[y]);
    /* store in columns */
    Mat_set_col (tmp, y, v);
    Vec_delete (v);
  }

  /* create the output matrix */
  output = __Mat_new (elem_size, height, width);

  /* apply on all columns (stored in rows of 'tmp') */
  it_transform_set_size (transform, height);
  for (x = 0; x < width; x++) {
    v = it_itransform (transform, tmp[x]);
    /* store in columns */
    Mat_set_col (output, x, v);
    Vec_delete (v);
  }
  Mat_delete (tmp);

  return (output);
}


static void __separable2D_get_output_size (it_transform2D_t * transform,
					   idx_t * width, idx_t * height)
{
  it_separable2D_t *separable = IT_SEPARABLE2D (transform);

  it_transform_get_output_size (separable->transform, width);
  it_transform_get_output_size (separable->transform, height);
}

static void __separable2D_set_size (it_transform2D_t * transform,
				    idx_t width, idx_t height)
{
  it_separable2D_t *separable = IT_SEPARABLE2D (transform);

  it_transform_set_size (separable->transform, width);
  it_transform_set_size (separable->transform, height);
}

static void __separable2D_get_size (it_transform2D_t * transform,
				    idx_t * width, idx_t * height)
{
  it_separable2D_t *separable = IT_SEPARABLE2D (transform);

  it_transform_get_size (separable->transform, width);
  it_transform_get_size (separable->transform, height);
}


it_instanciate (it_separable2D_t)
{
  it_new_args_start ();
  it_construct (it_transform2D_t);
  it_set_magic (it_this, it_separable2D_t);

  IT_TRANSFORM2D (it_this)->transform = __it_separable2D_transform;
  IT_TRANSFORM2D (it_this)->itransform = __it_separable2D_itransform;
  IT_TRANSFORM2D (it_this)->get_output_size = __separable2D_get_output_size;
  IT_TRANSFORM2D (it_this)->set_size = __separable2D_set_size;
  IT_TRANSFORM2D (it_this)->get_size = __separable2D_get_size;

  it_this->transform = it_new_args_next (it_transform_t *);

  it_new_args_stop ();

  return (it_this);
}
