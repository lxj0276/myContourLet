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
  Parser for Parameter files and command line arguments
  Copyright (C) 2005 Herve Jegou
*/

#ifndef __it_parser_h
#define __it_parser_h

#include <it/vec.h>
#include <it/mat.h>
#include <it/cplx.h>

#ifdef __cplusplus
extern "C" {
#endif


typedef char * parser_t;

/*----------------------------------------------------------------------*/

/* Init a parser from a parameter file                                  */
parser_t * parser_init_file( const char * filename );

/* Init a parser with the command line arguments. Since the first element 
 of argv is the program name, this element is discarded                 */
parser_t * parser_init_params( int argc, char ** argv );

/* Init a parser from a string                                          */
parser_t * parser_init_string( char * s );

/* The following functions add another source to an existing parser     */
parser_t * _parser_add_file( parser_t * p, const char * filename );
parser_t * _parser_add_params( parser_t * p, int argc, char ** argv );
parser_t * _parser_add_string( parser_t * p, char * s );

#define parser_add_file( p, filename )     do { p = _parser_add_file( p, filename );}   while(0)
#define parser_add_params( p, argc, argv ) do { p = _parser_add_params( p, argc, argv );} while(0)
#define parser_add_string( p, s )          do { p = _parser_add_string( p, s );} while(0)

/* Free the memory allocated for a parser                               */
void parser_delete( parser_t * p );

/* Print the content of the parser                                      */
void parser_print( parser_t * p );

/*----------------------------------------------------------------------*/

/* Retrieve the variable whose name is given by varname                 */
int parser_get_int( const parser_t * p, const char * varname );
double parser_get_double( const parser_t * p, const char * varname );
byte parser_get_byte( const parser_t * p, const char * varname ); 
char * parser_get_string( const parser_t * p, const char * varname ); 
cplx parser_get_cplx( const parser_t * p, const char * varname ); 

/* Retrieve the vectors whose name is given by varname                  */
vec  parser_get_vec( const parser_t * p, const char * varname );
ivec parser_get_ivec( const parser_t * p, const char * varname );
bvec parser_get_bvec( const parser_t * p, const char * varname );
cvec parser_get_cvec( const parser_t * p, const char * varname );

/* Retrieve the matrices whose name is given by varname                 */
mat parser_get_mat( const parser_t * p, const char * varname );
imat parser_get_imat( const parser_t * p, const char * varname );
bmat parser_get_bmat( const parser_t * p, const char * varname );
cmat parser_get_cmat( const parser_t * p, const char * varname );

/* Return 1 if the string has been found as a valid identifier, 0 otherwise */
int parser_exists( const parser_t * p, const char * varname );

/*----------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* extern "C" */
#endif 
