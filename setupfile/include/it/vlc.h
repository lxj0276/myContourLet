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
  Variable Length Codes
  Copyright (C) 2005 Herve Jegou
*/

#ifndef __it_vlc_h
#define __it_vlc_h

#include <it/vec.h>

#ifdef __cplusplus
extern "C" {
#endif

#define DUMM_NODE (-1)

/*----------------------------------------------------------------------*/
/* The vlc_t structure 
   Note: the height of the VLC codetree is limited to sizeof(int)-1     */
typedef struct {
  int nb_symb;         /* Number of symbols defined by the VLC codetree */
  int nb_nodes;        /* The total number of nodes                     */

  ivec cwd;            /* The codewords, stored in a compact form       */
  ivec cwd_length;     /* The length of the corresponding symbols       */
  ivec map[ 2 ];       /* The mapping between nodes                     */
  
  int node_root;       /* The index of the node root                    */
} vlc_t;


/*----------------------------------------------------------------------*/
/*! Initialization of a new VLC code                                    */
vlc_t * vlc_new( int n );

/*! Create a clone of a given VLC                                       */ 
vlc_t * vlc_clone( const vlc_t * vlc );

/*! Delete a VLC code                                                   */
void vlc_delete( vlc_t * vlc );

/*! Return the integer code associated with a node.                     */
static inline int vlc_get_cwd( const vlc_t * vlc, int node ) { return vlc->cwd[ node ]; }

/*! Return the length of integer code associated with a node.           */
static inline int vlc_get_cwd_length( const vlc_t * vlc, int node ) { return vlc->cwd_length[ node ]; }

/*! The child node of s if the next bit transition is a 0               */
static inline int vlc_get_child0( const vlc_t * vlc, int s ) { return vlc->map[ 0 ][ s ]; }

/*! The child node of s if the next bit transition is a 1               */
static inline int vlc_get_child1( const vlc_t * vlc, int s ) { return vlc->map[ 1 ][ s ]; }

/* Return 1 if s is a leaf, 0 otherwise                                 */
static inline int vlc_is_leaf( const vlc_t * vlc, int s ) { return s < vlc->nb_symb; }

/* Return 1 if s is a internal node, 0 otherwise                        */
static inline int vlc_is_node( const vlc_t * vlc, int s ) { return s >= vlc->nb_symb; }


/*----------------------------------------------------------------------*/
/* VLC design                                                           */

/* Sort the VLC so that it is the "more" lexicographic as possible      */
void vlc_quasi_lexicographic( const vlc_t * vlc, 
			      const vec pdf, 
			      const vec symb );

/* Return the Huffman code corresponding to the probability law pdf     */
vlc_t * vlc_huffman( const vec pdf );

/* Return the Hu-Tucker code corresponding to the probability law pdf,
 assuming that these symbols are sorted in the natural order            */
vlc_t * vlc_hu_tucker( const vec pdf );

/* Return a Fixed length Code for the number of bits                    */
vlc_t * vlc_flc( int nb_bits );

/* Read a variable length code from a string                            */
vlc_t * vlc_read( const char * s );

/* Assign the codewords according to the definition of variable map     */
void vlc_affect_codewords( vlc_t * vlc );


/*----------------------------------------------------------------------*/
/* VLC performance and measures                                         */

/* The length of the shortest codeword                                  */
int vlc_minh( const vlc_t * vlc );

/* The length of the longest codeword                                   */
int vlc_maxh( const vlc_t * vlc );

/* Kraft sum. (sum_s 2^{-L(s)})                                         */
double vlc_kraft_sum( const vlc_t * vlc );

/* Return the average length of the code (i.e. the expectation 
   of the length) if the distribution of the source is pdf              */
double vlc_mdl( const vlc_t * vlc, const vec pdf );                        

/* Return the respective probabilities of the nodes                     */
vec vlc_nodes_pdf( const vlc_t * vlc, vec pdf );

/* Return the respective expectations of the nodes                      */
vec vlc_nodes_expectation( const vlc_t * vlc, vec pdf, vec symb );

/* Return the value of variance when being in a node                    */
double vlc_node_variance( const vlc_t * vlc, int node, vec pdf, vec symbols );
vec vlc_nodes_variance( const vlc_t * vlc, vec pdf, vec symb );

/* Return the binary entropy of the internal nodes of the VLC           */
vec vlc_nodes_entropy( const vlc_t * vlc, vec pdf );

/* Expectation of the decrease of MSE corresponding to internal nodes   */
vec vlc_nodes_delta_energy( const vlc_t * vlc, vec pdf, vec symb );

/* Return the order between node so that each transmitted bit is the one
   that minimizes the rate-distortion curve (assuming that the VLC 
   is encoded with an optimal entropy coder afterward                   */
ivec vlc_energy_order( const vlc_t * vlc, vec pdf, vec symb );

/*----------------------------------------------------------------------*/
/* VLC I/O                                                              */

void vlc_print( const vlc_t * vlc );
void vlc_print_all( const vlc_t * vlc, vec pdf, vec symb );

#ifdef __cplusplus
}
#endif /* extern "C" */
#endif
