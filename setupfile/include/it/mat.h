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
  Matrices functions
  Copyright (C) 2005 Vivien Chappelier, Herve Jegou
*/

#ifndef __it_mat_h
#define __it_mat_h

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "it/types.h"
#include "it/vec.h"

#ifdef __cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------*/
/*                          Matrices                               */
/*-----------------------------------------------------------------*/

/* Like vectors, the pointer corresponding to a given matrix 
   represents the first component ([0][0]) of a given matrix.
   M being a matrix, the indexation convention is M[row][col]. 
   As the vector length, the number of rows is stored before the components.    */

typedef Vec * Mat;
typedef vec * mat;
typedef ivec * imat;
typedef bvec * bmat;
typedef cvec * cmat;

/*-----------------------------------------------------------------*/
/*                Allocation, initialisation and free                           */

/* Return the length of the vector and the number of allocated elements         */
#define Mat_height( m )     Vec_length(m)
#define Mat_height_max( m ) Vec_length_max(m)

#define Mat_width( m )      Vec_length(m[0])
#define Mat_width_max( m )  Vec_length_max(m[0])

/* Allocate a matrix with explicit allocation of the maximum height and width   */
#define Mat_new_alloc( type_t, h, w, hmax, wmax ) \
          ((type_t **) __Mat_new_alloc( sizeof(type_t), h, w, hmax, wmax ))

static inline Mat __Mat_new_alloc(size_t elem_size, size_t h, size_t w, size_t hmax, size_t wmax)
{
  size_t i;
  Mat m;
  
  m = Vec_new_alloc( Vec, h, hmax );
  for( i = 0 ; i < hmax ; i++ )
    m[i] = __Vec_new_alloc( elem_size, w, wmax );
  return(m);
}

/* Allocate a matrix. This function do not return the matrix pointer.
   For this purpose, use specialized matrix instead                             */
#define Mat_new( type_t, h, w ) ((type_t **) __Mat_new( sizeof(type_t), h, w ))
static inline Mat __Mat_new( size_t elem_size, size_t h, size_t w)
{
  return(__Mat_new_alloc(elem_size, h, w, h, w));
}

/* Free the vector                                                              */
#define Mat_delete( m ) __Mat_delete((Mat) m)
static inline void __Mat_delete(Mat m)
{
  size_t i;
  for( i = 0 ; i < Mat_height_max(m) ; i++ )
    Vec_delete(m[i]);
  Vec_delete(m);
}

/* Modify the height of the matrix (but never reduce the allocated memory)      */
#define Mat_set_height( m, h )     do {                                          \
                                      if( h > Mat_height( m ) ) {                \
                                        Mat_set_height_max( m, h );              \
                                        Vec_length(m) = (h);                     \
                                      }                                          \
                                      else                                       \
                                        Vec_length(m) = (h);                     \
                                   } while(0)
 
#define Mat_set_height_max( m, hmax ) do {                                       \
                                     size_t oldhmax = Mat_height_max(m), i;      \
                                     assert( hmax >= Vec_length(m) );            \
                                     if( hmax > Mat_height(m) ) {                \
                                       Vec_set_length_max(m,hmax);               \
	                                 for( i = oldhmax ; i < hmax ; i++ )         \
                                       Vec_set_length(m[i],Mat_width(m));        \
                                     } } while(0)


/* Initialize the matrix from a buffer of the same type            */
#define Mat_init( m, buf, w, h ) do {                                            \
                                    size_t i;                                    \
                                    Mat_set_height(m, h);                        \
                                    for( i=0 ; i < h ; i++ )                     \
                                      Vec_init( m[i], buf+(w)*(i), w );          \
                                 } while(0)

/* Set all the elements of the vector to value val                 */
#define Mat_set( m, val ) do {                                                  \
                              size_t i,j;                                       \
                              assert( m );                                      \
                              for( i = 0 ; i < Mat_height(m) ; i++ )            \
                                for( j = 0 ; j < Mat_width(m) ; j++ )           \
                                  m[i][j] = val;                                \
                             } while(0)

/* to use the entity end in a mat function, this macro is called */
#define MAT_END_COL_PARAM( v, param ) do { if( param == (size_t) end ) param = Mat_width( v ) - 1; } while(0)
#define MAT_END_ROW_PARAM( v, param ) do { if( param == (size_t) end ) param = Mat_height( v ) - 1; } while(0)
     
/*-----------------------------------------------------------------*/
/* double matrices functions                                       */

static inline mat mat_new( size_t h, size_t w ) { return(Mat_new(double, h, w)); }
static inline imat imat_new( size_t h, size_t w ) { return(Mat_new(int, h, w)); }
static inline bmat bmat_new( size_t h, size_t w ) { return(Mat_new(byte, h, w)); }
static inline cmat cmat_new( size_t h, size_t w ) { return(Mat_new(cplx, h, w)); }

static inline mat mat_new_alloc( size_t h, size_t w, size_t hmax, size_t wmax ) { return(Mat_new_alloc( double, h, w, hmax, wmax)); }
static inline imat imat_new_alloc( size_t h, size_t w, size_t hmax, size_t wmax ) { return(Mat_new_alloc( int, h, w, hmax, wmax)); }
static inline bmat bmat_new_alloc( size_t h, size_t w, size_t hmax, size_t wmax ) { return(Mat_new_alloc( byte, h, w, hmax, wmax)); }
static inline cmat cmat_new_alloc( size_t h, size_t w, size_t hmax, size_t wmax ) { return(Mat_new_alloc( cplx, h, w, hmax, wmax)); }

static inline void mat_delete( mat m ) { Mat_delete(m); }
static inline void imat_delete( imat m ) { Mat_delete(m); }
static inline void bmat_delete( bmat m ) { Mat_delete(m); }
static inline void cmat_delete( cmat m ) { Mat_delete(m); }

static inline size_t mat_height( mat m ) { return Mat_height(m); }
static inline size_t mat_height_max( mat m ) { return Mat_height_max(m); }
static inline size_t mat_width( mat m ) { return Mat_width(m); } 
static inline size_t mat_width_max( mat m ) { return Mat_width_max(m); }

static inline size_t imat_height( imat m ) { return Mat_height(m); }
static inline size_t imat_height_max( imat m ) { return Mat_height_max(m); }
static inline size_t imat_width( imat m ) { return Mat_width(m); } 
static inline size_t imat_width_max( imat m ) { return Mat_width_max(m); }

static inline size_t bmat_height( bmat m ) { return Mat_height(m); }
static inline size_t bmat_height_max( bmat m ) { return Mat_height_max(m); }
static inline size_t bmat_width( bmat m ) { return Mat_width(m); } 
static inline size_t bmat_width_max( bmat m ) { return Mat_width_max(m); }

static inline size_t cmat_height( cmat m ) { return Mat_height(m); }
static inline size_t cmat_height_max( cmat m ) { return Mat_height_max(m); }
static inline size_t cmat_width( cmat m ) { return Mat_width(m); } 
static inline size_t cmat_width_max( cmat m ) { return Mat_width_max(m); }


static inline mat _mat_set_height( mat m, size_t h ) { Mat_set_height(m, h); return m; } 
static inline mat _mat_set_height_max( mat m, size_t hmax ){ Mat_set_height_max(m,hmax); return m; }
#define mat_set_height( m, h ) do { m = _mat_set_height(m, h); } while(0)
#define mat_set_height_max( m, h ) do { m = _mat_set_height_max(m, h); } while(0)

static inline imat _imat_set_height( imat m, size_t h ) { Mat_set_height(m, h); return m; } 
static inline imat _imat_set_height_max( imat m, size_t hmax ){ Mat_set_height_max(m,hmax); return m; }
#define imat_set_height( m, h ) do { m = _imat_set_height(m, h); } while(0)
#define imat_set_height_max( m, h ) do { m = _imat_set_height_max(m, h); } while(0)

static inline bmat _bmat_set_height( bmat m, size_t h ) { Mat_set_height(m, h); return m; } 
static inline bmat _bmat_set_height_max( bmat m, size_t hmax ){ Mat_set_height_max(m,hmax); return(m); }
#define bmat_set_height( m, h ) do { m = _bmat_set_height(m, h); } while(0)
#define bmat_set_height_max( m, h ) do { m = _bmat_set_height_max(m, h); } while(0)

static inline cmat _cmat_set_height( cmat m, size_t h ) { Mat_set_height(m, h); return m; } 
static inline cmat _cmat_set_height_max( cmat m, size_t hmax ){ Mat_set_height_max(m,hmax); return(m); }
#define cmat_set_height( m, h ) do { m = _cmat_set_height(m, h); } while(0)
#define cmat_set_height_max( m, h ) do { m = _cmat_set_height_max(m, h); } while(0)

static inline mat  _mat_init( mat m, double * buf, size_t w, size_t h ) { Mat_init(m, buf, w, h); return(m); }
static inline imat _imat_init( imat m, int * buf, size_t w, size_t h ) { Mat_init(m, buf, w, h); return(m); }
static inline bmat _bmat_init( bmat m, byte * buf, size_t w, size_t h ) { Mat_init(m, buf, w, h); return(m); }
static inline cmat _cmat_init( cmat m, cplx * buf, size_t w, size_t h ) { Mat_init(m, buf, w, h); return(m); }

#define mat_init(m, buf, w, h) _mat_init(m, buf, w, h) 
#define imat_init(m, buf, w, h) _imat_init(m, buf, w, h) 
#define bmat_init(m, buf, w, h) _bmat_init(m, buf, w, h) 
#define cmat_init(m, buf, w, h) _cmat_init(m, buf, w, h) 

static inline mat  mat_set( mat m, double val ) { Mat_set(m,val); return m; }
static inline imat imat_set( imat m, int val )  { Mat_set(m,val); return m; }
static inline bmat bmat_set( bmat m, byte val ) { Mat_set(m,val); return m;}
static inline cmat cmat_set( cmat m, cplx val ) { Mat_set(m,val); return m;}

void mat_copy( mat dest, mat orig );
void imat_copy( imat dest, imat orig );
void bmat_copy( bmat dest, bmat orig );
void cmat_copy( cmat dest, cmat orig );

mat  mat_clone( mat m );
imat imat_clone( imat m );
bmat bmat_clone( bmat m );
cmat cmat_clone( cmat m );

/*------------------------------------------------------------------*/
/* Measures for matrices                                            */
double mat_sum( mat m );
int    imat_sum( imat m );
cplx   cmat_sum( cmat m );

double mat_row_sum( mat m, size_t c );
double mat_col_sum( mat m, size_t c );
int    imat_row_sum( imat m, size_t c );
int    imat_col_sum( imat m, size_t c );
cplx   cmat_row_sum( cmat m, size_t c );
cplx   cmat_col_sum( cmat m, size_t c );

vec  mat_rows_sum( mat m );
vec  mat_cols_sum( mat m );
ivec imat_rows_sum( imat m );
ivec imat_cols_sum( imat m );
cvec cmat_rows_sum( cmat m );
cvec cmat_cols_sum( cmat m );

double mat_mean( mat m );
double imat_mean( imat m );
cplx   cmat_mean( cmat m );

double mat_max( mat m );
int imat_max( imat m );

double mat_min( mat m );
int imat_min( imat m );

double mat_variance( mat m );
void mat_normalize( mat m );

/* Matrice norms */
double mat_norm_1( mat m );
double mat_norm_inf( mat m );

/*------------------------------------------------------------------*/
/* Arithmetic matrix functions                                      */

/* Operations with a scalar value                                               */
mat  mat_incr( mat m, double a );       /* add a to the components              */
mat  mat_decr( mat m, double a );       /* substract a to the components        */
mat  mat_mul_by( mat m, double a );     /* multiply the components per a        */
mat  mat_div_by( mat m, double a );     /* divide the components per a          */

imat imat_incr( imat m, int a );        /* add a to the components              */
imat imat_decr( imat m, int a );        /* substract a to the components        */
imat imat_mul_by( imat m, int a );      /* multiply the components per a        */
imat imat_div_by( imat m, int a );      /* divide the components per a          */

mat  mat_col_incr( mat m, size_t col, double a ); 
mat  mat_col_decr( mat m, size_t col, double a ); 
mat  mat_col_mul_by( mat m, size_t col, double a );
mat  mat_col_div_by( mat m, size_t col, double a ); 

imat imat_col_incr( imat m, size_t col, int a );  
imat imat_col_decr( imat m, size_t col, int a );  
imat imat_col_mul_by( imat m, size_t col, int a ); 
imat imat_col_div_by( imat m, size_t col, int a ); 

mat  mat_row_incr( mat m, size_t row, double a ); 
mat  mat_row_decr( mat m, size_t row, double a ); 
mat  mat_row_mul_by( mat m, size_t row, double a );
mat  mat_row_div_by( mat m, size_t row, double a ); 

imat imat_row_incr( imat m, size_t row, int a );  
imat imat_row_decr( imat m, size_t row, int a );  
imat imat_row_mul_by( imat m, size_t row, int a ); 
imat imat_row_div_by( imat m, size_t row, int a ); 


/* Components per components operations (vectors must be of same size).
   The only difference between mat_add and mat_elem_add is that the former 
   allocate a new matrix, where in the latter return the address of m1          */
mat  mat_elem_add( mat m1, mat m2 );      
mat  mat_elem_sub( mat m1, mat m2 );
mat  mat_elem_mul( mat m1, mat m2 );
mat  mat_elem_div( mat m1, mat m2 );
     
imat imat_elem_add( imat m1, imat m2 );      
imat imat_elem_sub( imat m1, imat m2 );
imat imat_elem_mul( imat m1, imat m2 );
imat imat_elem_div( imat m1, imat m2 );

/* Matrix to matrix, matrix to vector and vector to matrix 
   addition, substraction and multiplication                        */
mat  mat_add( mat m1, mat m2 );      
mat  mat_sub( mat m1, mat m2 );
mat  mat_mul( mat m1, mat m2 );
     
imat imat_add( imat m1, imat m2 );      
imat imat_sub( imat m1, imat m2 );
imat imat_mul( imat m1, imat m2 );

vec  mat_vec_mul( mat m, vec v );
vec  mat_ivec_mul( mat m, ivec v );
vec  imat_vec_mul( imat m, vec v );
ivec imat_ivec_mul( imat m, ivec v );

vec  vec_mat_mul( vec v, mat m );
ivec ivec_imat_mul( ivec v, imat m );

/* Transposition of the matrix. Note: the functions _mat_transpose 
   and _imat_transpose can be called directly only for symmetric matrices. 
   We advise to use only the macros instead (work in the general case). */
mat _mat_transpose( mat m );
imat _imat_transpose( imat m );
bmat _bmat_transpose( bmat m );

#define mat_transpose( m ) do { 		\
  if( mat_height( m ) == mat_width( m ) )	\
    m = _mat_transpose( m );			\
  else {					\
    mat t = mat_new_transpose( m );		\
    mat_delete( m );				\
    m = t;					\
  } } while(0)

#define imat_transpose( m ) do { 		\
  if( imat_height( m ) == imat_width( m ) )	\
    m = _imat_transpose( m );			\
  else {					\
    imat t = imat_new_transpose( m );		\
    imat_delete( m );				\
    m = t;					\
  } } while(0)

#define bmat_transpose( m ) do { 		\
  if( bmat_height( m ) == bmat_width( m ) )	\
    m = _bmat_transpose( m );			\
  else {					\
    bmat t = bmat_new_transpose( m );		\
    bmat_delete( m );				\
    m = t;					\
  } } while(0)


mat mat_new_transpose( mat m );
imat imat_new_transpose( imat m );
bmat bmat_new_transpose( bmat m );



/* Inversion of a matrix using the Gaussian method */
mat mat_inv( mat m );
mat mat_new_inv( mat m );

/*------------------------------------------------------------------*/
/* Retrieve part or assemble matrices                               */
mat  mat_get_submatrix( mat m, size_t r1, size_t c1, size_t r2, size_t c2 );
imat imat_get_submatrix( imat m, size_t r1, size_t c1, size_t r2, size_t c2 );
bmat bmat_get_submatrix( bmat m, size_t r1, size_t c1, size_t r2, size_t c2 );

mat  mat_set_submatrix( mat m, mat s, size_t r, size_t c );
imat imat_set_submatrix( imat m, imat s, size_t r, size_t c );
bmat bmat_set_submatrix( bmat m, bmat s, size_t r, size_t c );

static inline vec mat_copy_row( vec v, mat m, size_t r) {
  assert(m);
  MAT_END_ROW_PARAM( m, r );
  assert(r < mat_height(m));
  assert(vec_length(v) == mat_width(m));
  return(vec_copy(v, m[r]));
}
static inline ivec imat_copy_row( ivec v, imat m, size_t r) {
  assert(m);
  MAT_END_ROW_PARAM( m, r );
  assert(r < imat_height(m));
  assert(ivec_length(v) == imat_width(m));
  return(ivec_copy(v, m[r]));
}
static inline bvec bmat_copy_row( bvec v, bmat m, size_t r) {
  assert(m);
  MAT_END_ROW_PARAM( m, r );
  assert(r < bmat_height(m));
  assert(bvec_length(v) == bmat_width(m));
  return(bvec_copy(v, m[r]));
}
static inline cvec cmat_copy_row( cvec v, cmat m, size_t r) {
  assert(m);
  MAT_END_ROW_PARAM( m, r );
  assert(r < cmat_height(m));
  assert(cvec_length(v) == cmat_width(m));
  return(cvec_copy(v, m[r]));
}

vec mat_copy_col( vec v, mat m, size_t c);
ivec imat_copy_col( ivec v, imat m, size_t c);
bvec bmat_copy_col( bvec v, bmat m, size_t c);
cvec cmat_copy_col( cvec v, cmat m, size_t c);

static inline vec mat_get_row( mat m, size_t r) { vec v; assert(m); v = vec_new(mat_width(m)); return(mat_copy_row(v, m, r)); }
static inline ivec imat_get_row( imat m, size_t r) { ivec v; assert(m); v = ivec_new(imat_width(m)); return(imat_copy_row(v, m, r)); }
static inline bvec bmat_get_row( bmat m, size_t r) { bvec v; assert(m); v = bvec_new(bmat_width(m)); return(bmat_copy_row(v, m, r)); }
static inline cvec cmat_get_row( cmat m, size_t r) { cvec v; assert(m); v = cvec_new(cmat_width(m)); return(cmat_copy_row(v, m, r)); }

static inline vec mat_get_col( mat m, size_t c) { vec v; assert(m); v = vec_new(mat_height(m)); return(mat_copy_col(v, m, c)); }
static inline ivec imat_get_col( imat m, size_t c) { ivec v; assert(m); v = ivec_new(imat_height(m)); return(imat_copy_col(v, m, c)); }
static inline bvec bmat_get_col( bmat m, size_t c) { bvec v; assert(m); v = bvec_new(bmat_height(m)); return(bmat_copy_col(v, m, c)); }
static inline cvec cmat_get_col( cmat m, size_t c) { cvec v; assert(m); v = cvec_new(cmat_height(m)); return(cmat_copy_col(v, m, c)); }

void Mat_set_col(Mat m, size_t c, Vec v);
void mat_set_col(mat m, size_t c, vec v);
void imat_set_col(imat m, size_t c, ivec v);
void bmat_set_col(bmat m, size_t c, bvec v);
void cmat_set_col(cmat m, size_t c, cvec v);

/*------------------------------------------------------------------*/
/* Application of general functions                                 */
void mat_apply_function( mat v, it_function_t function, it_args_t args ); 
mat mat_new_apply_function( mat v, it_function_t function, it_args_t args );
#define mat_eval(v, f, a) mat_apply_function(v, f, a)
#define mat_new_eval(v, f, a) mat_new_apply_function(v, f, a)

void imat_apply_function( imat v, it_ifunction_t function, it_args_t args ); 
imat imat_new_apply_function( imat v, it_ifunction_t function, it_args_t args );
#define imat_eval(v, f, a) imat_apply_function(v, f, a)
#define imat_new_eval(v, f, a) imat_new_apply_function(v, f, a)


/*------------------------------------------------------------------
 Special matrices
 Two versions of the following functions exist. The first set 
 modify already existing matrices. The other allocate new matrices  */


mat mat_void();
imat imat_void();
bmat bmat_void();
cmat cmat_void();

mat mat_zeros( mat m );
imat imat_zeros( imat m );
bmat bmat_zeros( bmat m );
cmat cmat_zeros( cmat m );

mat mat_ones( mat m );
imat imat_ones( imat m );
bmat bmat_ones( bmat m );
cmat cmat_ones( cmat m );

/* Set a matrix to identity. This function also works if the matrix is not square */
mat mat_eye( mat m );
imat imat_eye( imat m );
bmat bmat_eye( bmat m );
cmat cmat_eye( cmat m );

/* Construct a diagonal matrix */
mat mat_diag( mat m, vec v );
imat imat_diag( imat m, ivec v );
bmat bmat_diag( bmat m, bvec v );

/*------------------------------------------------------------------*/

mat mat_new_void();
imat imat_new_void();
bmat bmat_new_void();
cmat cmat_new_void();

mat mat_new_zeros( size_t h, size_t w );
imat imat_new_zeros( size_t h, size_t w );
bmat bmat_new_zeros( size_t h, size_t w );
cmat cmat_new_zeros( size_t h, size_t w );

mat mat_new_ones( size_t h, size_t w );
imat imat_new_ones( size_t h, size_t w );
bmat bmat_new_ones( size_t h, size_t w );
cmat cmat_new_ones( size_t h, size_t w );

mat mat_new_eye( size_t n );
imat imat_new_eye( size_t n );
bmat bmat_new_eye( size_t n );
cmat cmat_new_eye( size_t n );

mat mat_new_diag( vec v );
imat imat_new_diag( ivec v );
bmat bmat_new_diag( bvec v );

/*------------------------------------------------------------------*/
/* Conversion functions                                             */

/* Vectorize the matrix into a vector                               */
vec mat_to_vec( mat m );
ivec imat_to_ivec( imat m );
mat imat_to_mat( imat m );

/* Transform a vector into a matrix, assuming a given matrix width  */
mat vec_to_mat( vec v, size_t width );
imat ivec_to_imat( ivec v, size_t width );


#ifdef __cplusplus
}
#endif /* extern "C" */
#endif




