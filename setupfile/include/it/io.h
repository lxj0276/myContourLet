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
  I/O functions for the libit library
  Copyright (C) 2005 Herve Jegou
*/

#ifndef __it_io_h
#define __it_io_h

#include <stdarg.h>
#include <assert.h>
#include "it/mat.h"

#ifdef __cplusplus
extern "C" {
#endif
/*----------------------------------------------------------------------*/
/* Extending printing of formatted string (includes vector and matrix)  */
void it_vprintf( char *fmt, va_list args );
void it_fprintf( FILE * output, char *fmt, ...);
void it_printf( char *fmt, ... );
void it_vfprintf( FILE * output, char *fmt, va_list ap );

/* Error functions                                                      */
void it_warning( char *fmt, ... );
void it_error( char *fmt, ... );

void __it_assert( int a, const char * filename, int line, char * message );
#ifdef NDEBUG
#define it_assert(a,expr) ((void) (0))
#else
#if defined(__GNUC__)
#define it_assert(a,expr) __it_assert((int)(a), __FILE__, __LINE__, expr)
#else
#define it_assert(a,expr) __it_assert((int)(a), "", 0, expr)
#endif
#endif

/* Modify the default format for matrix printing (called by format $v)  */
void it_set_vec_default_fmt( char * fmt );

/* Modify the default format for vector printing (called by format $v)  */
void it_set_mat_default_fmt( char * fmt );


/* Read integers and double from an NOT null-terminated string */
char * it_read_int( char * s, int * p_val );
char * it_read_double( char * s, double * p_val );
char * it_read_cplx( char * s, cplx * p_val );

/* Read a vector from a string and return the position after the vector */
char * it_read_vec( char * s, vec * v );
char * it_read_ivec( char * s, ivec * v );
char * it_read_bvec( char * s, bvec * v );
char * it_read_cvec( char * s, cvec * v );

/* Create a vector from a string */
vec vec_new_string( char * s );
ivec ivec_new_string( char * s );
bvec bvec_new_string( char * s );
cvec cvec_new_string( char * s );

/*----------------------------------------------------------------------*/
/* Determines the type of a PNM file, based on magic numbers.
   The returned value is a character between '1' and '6'.
   If an error occured, the returned value the value 0                  */
char pnm_type( const char * filename );

/* Retrieve some information about an pnm file.                         
   String comments must be pre-allocated                                */
int pnm_info( const char * filename, char * p_pnm_type, int * p_width, 
	       int * p_height, int * p_max_val, char * comments );

/* Read a pgm file and return the corresponding matrix of double values */
mat mat_pgm_read( const char* filename );

/* Read a pgm file and return the corresponding matrix of int values */
imat imat_pgm_read( const char* filename );

/* Write a matrix of double as a pgm image file                         */
int mat_pgm_write( const char* filename, mat m );

/* Write a matrix of integers as a pgm file                             */
int imat_pgm_write( const char* filename, imat m );

/*----------------------------------------------------------------------*/
/* WAV file handling functions */
int wav_info( const char * filename, int * p_channels, int *p_srate, int *p_depth, int *p_length);
mat mat_wav_read( const char* filename );
imat imat_wav_read( const char* filename );
void mat_wav_write( const char* filename, mat m, int srate, int depth );
void imat_wav_write( const char* filename, imat m, int srate, int depth );

#ifdef __cplusplus
}
#endif /* extern "C" */
#endif 
