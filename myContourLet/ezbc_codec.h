/* libezbc - EZBC subband coding/decoding
      Copyright (C) 2004 Vivien Chappelier

      This library is inspired from code originaly developped by Yongjun Wu Jan in MC-EZBC. This rewritten version may still contain code similar to his original implementation.

      Reference paper: "Embedded Image Coding Using Zeroblocks Of Subband/Wavelet Coefficient And Context Modeling", Shih-Ta Hsiang and John W. Woods, Data Compression Conference (DCC '01), March 27 - 29, 2001, Snowbird, Utah 
*/

#ifndef __EZBC_CODEC_H
#define __EZBC_CODEC_H

#include "ezbc_types.h"
#include "ezbc_tables.h"
#include "arithmetic_codec.h"

#define GET_PARENT_MODELS

#define LAST_NSIG_IN_GROUP (!sig && (x & 1) && (y & 1))

/* after BINARY_MODEL_MAX events, rescale the counters */
#define BINARY_MODEL_MAX 512

/* context model */
typedef struct {
  int c0; /* count of zeros */
  int c1; /* count of ones */
} binary_model_t;

#define binary_model_init(m) { (m).c0 = (m).c1 = 1; }
#define binary_model_update(m, b) { if(b) (m).c1 ++; else (m).c0 ++; };
#define binary_model_p0(m) (((m).c0 + (m).c1)?((double) (m).c0 / ((m).c0 + (m).c1)):0.5)
#define binary_model_p1(m) (((m).c0 + (m).c1)?((double) (m).c1 / ((m).c0 + (m).c1)):0.5)
#define binary_model_scale(m) while((m).c0 + (m).c1 > BINARY_MODEL_MAX) { (m).c0 = ((m).c0 + 1) >> 1; (m).c1 = ((m).c1 + 1) >> 1; }
#define binary_model_copy(dst, src) do {	\
  (dst).c0 = (src).c0;				\
  (dst).c1 = (src).c1;				\
} while(0);

/* contexts */
typedef struct {
  int offset;
  int count;
  ezbc_byte_t *table;
} context_t;

/* node stack element */
typedef struct {
  int x;
  int y;
  int l;
} ezbc_node_t;

typedef struct _ezbc_subband_codec_ {
  int width;             /* width of the subband */
  int height;            /* height of the subband */
  int depth;             /* depth of the quad tree */
  int bpp;
  int target_rate;       /* target rate in bits */

  ezbc_coeff_t *coeffs;   /* subband coefficients */

  ezbc_point_t *LIP;     /* list of insignificant pixels */
  ezbc_point_t *LIP_end; /* pointer to the end of the list of insignificant pixels */
  ezbc_point_t *LSP;     /* list of significant pixels */
  ezbc_point_t *LSP_end; /* pointer to the end of the list of significant pixels */
  ezbc_point_t **LIS;     /* list of insignificant sets (qtree) */
  ezbc_point_t **LIS_end; /* pointer to the end of the list of insignificant sets (one per level) */

  ezbc_node_t *LIS_stack;      /* stack of insignificant sets to code */
  ezbc_node_t *LIS_stack_top;  /* top of the stack insignificant sets to code */
  ezbc_point_t *LSP_split;     /* pointer to the first refined pixel */
  ezbc_point_t *LSP_mark;      /* number of pixels in each bitplane */
  int LSP_plane;               /* current refinement plane */

  ezbc_point_t **LSP_bit_index_marks;
  ezbc_point_t ***LSP_ids;

  context_t context_LIP;    /* list of insignificant pixels contexts */
  context_t context_LSP;    /* list of significant pixels contexts */
  context_t context_sign;   /* sign context quad tree */
  context_t *context_sig;    /* significant context quad tree */
  context_t *context_jsig;   /* just significant context quad tree */
  context_t *context_node;   /* node context quad tree */
  int context_count;
  binary_model_t *context_models;

  int **base_state;  /* current state of significant pixels */
  int **sign_state;  /* current state of sign */
  int ***node_state; /* current state of qtree nodes */
  ezbc_coeff_t ***nodes;      /* state of qtree nodes */

  ezbc_context_tables_t tables;

  arithmetic_codec_t *arithmetic_codec; /* arithmetic coder/decoder */

  struct _ezbc_subband_codec_ *parent_codec; /* codec of the parent band */
  struct _ezbc_subband_codec_ *child_codec; /* codec of the child band */

  /* the following variables are necessary to know when to
     stop decoding in case a VLC is used. */
  int *binbook; /* VLC: binarization book */
  int *hash;    /* VLC: index to the binarization book sorted by length */

} ezbc_subband_codec_t;

/* integer part of the log2 of v (number of bits - 1 needed to code v) */
static int inline intlog2(int v)
{
  int i;

  for(i = 0; v >>= 1; i++);

  return(i);
}

#define EZBC_ORIENT_DIAGONAL           0
#define EZBC_ORIENT_HORIZONTAL         1
#define EZBC_ORIENT_VERTICAL           2
#define EZBC_ORIENT_HORIZONTALVERTICAL 3
#define EZBC_ORIENT_OTHER 4
#define EZBC_ORIENT_NONE 5

void update_node_contexts(ezbc_subband_codec_t *codec);

/* update the zero coding state of neighbors */
void update_zc_state(int **base_state, int x, int y);

/* update the sign coding state of neighbors */
void update_sc_state(int **sign_state, int x, int y, int sign);

/* update the zero coding parent state of child (in other subband) */
void update_child_state(int **child_state, int x, int y);

void ezbc_subband_init(ezbc_subband_codec_t *codec,
		       ezbc_subband_codec_t *parent_codec,
		       ezbc_subband_codec_t *child_codec,
		       arithmetic_codec_t *arithmetic_codec,
		       ezbc_coeff_t *coeffs,
		       int width,
		       int height,
		       int orientation,
		       int bpp,
		       int target_rate);

void ezbc_subband_cleanup(ezbc_subband_codec_t *codec);

void ezbc_subband_start(ezbc_subband_codec_t *codec);

void ezbc_sort_list(ezbc_point_t *list, ezbc_coeff_t **coeffs, int length, int index);

void ezbc_init();

#endif
