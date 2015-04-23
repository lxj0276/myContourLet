/* libezbc - EZBC subband coding/decoding
      Copyright (C) 2004 Vivien Chappelier

      This library is inspired from code originaly developped by Yongjun Wu Jan in MC-EZBC. This rewritten version may still contain code similar to his original implementation.

      Reference paper: "Embedded Image Coding Using Zeroblocks Of Subband/Wavelet Coefficient And Context Modeling", Shih-Ta Hsiang and John W. Woods, Data Compression Conference (DCC '01), March 27 - 29, 2001, Snowbird, Utah 
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ezbc_types.h"
#include "ezbc_tables.h"
#include "ezbc_codec.h"

int vlc_is_leaf(ezbc_subband_codec_t *codec, int code);

#define qtree_C(l, x, y) ((1 << (2*l)) | (y << l) | (x))

/* decode a bit using a binary arithmetic codec */
static int decode_symbol(ezbc_subband_codec_t *codec, binary_model_t *model)
{
  int bit;
  arithmetic_codec_t *arithmetic_codec = codec->arithmetic_codec;

  bit = arithmetic_codec_decode(arithmetic_codec, binary_model_p0(*model));

  binary_model_update(*model, bit);

  return(bit);
}

/* decode the List of Insignificant Pixels */
void decode_LIP(ezbc_subband_codec_t *codec, int index)
{
  int width = codec->width;
  int height = codec->height;
  int x, y;
  ezbc_point_t cur_coord;
  ezbc_coeff_t coeff;
  ezbc_coeff_t ***nodes = codec->nodes;
  ezbc_point_t *LIP_cur, *LIP_end, *LIP_old_end, *LIP;
  ezbc_point_t *LSP_end, *LSP;

  int **base_state = codec->base_state;
  binary_model_t *LIP_models = &codec->context_models[codec->context_LIP.offset];
  ezbc_byte_t *LIP_context_table = codec->context_LIP.table;

  binary_model_t *sign_models = &codec->context_models[codec->context_sign.offset];
  ezbc_byte_t *sign_context_table = codec->context_sign.table;
  int **sign_state = codec->sign_state;
  int sign_predict, sign_bit;
  ezbc_coeff_t threshold;
  int **child_node_state = NULL; /* for good measure */
  int **child_base_state = NULL;

  if(codec->child_codec) {
    child_node_state = codec->child_codec->node_state[1];
    child_base_state = codec->child_codec->base_state;
  }

  /* compare to this threshold to determine signifiance */
  threshold = 1 << index; 

  LIP = codec->LIP;
  LIP_cur = LIP_end = LIP + width * height;
  LIP_old_end = codec->LIP_end;

  LSP = codec->LSP;
  LSP_end = codec->LSP_end;

  while(LIP_cur != LIP_old_end){
    cur_coord = *(--LIP_cur);
    x = ezbc_coord_x(cur_coord);
    y = ezbc_coord_y(cur_coord);

    coeff = nodes[0][y][x];
    if(decode_symbol(codec, &LIP_models[LIP_context_table[base_state[y][x] & ZC_MASK]])) {
      /* coefficient is significant */
      nodes[0][y][x] |= threshold;
      nodes[0][y][x] |= EZBC_MSB - index + 1; /* set the length */

      if(++LSP_end == LIP_old_end) /* LSP_pos runs beyond the end of the LSP */
        *(LIP_cur++) = *(LIP_old_end++);
      *LSP_end = cur_coord;

      /* update contexts */
      update_zc_state(base_state, x, y);

      if(codec->child_codec) {
	update_child_state(child_node_state, x, y);
	update_child_state(child_base_state, 2*x+0, 2*y+0); /* predict     */
	update_child_state(child_base_state, 2*x+1, 2*y+0); /* children    */
	update_child_state(child_base_state, 2*x+0, 2*y+1); /* are         */
	update_child_state(child_base_state, 2*x+1, 2*y+1); /* significant */
      }

      /* decode the sign */
      sign_predict = sign_context_table[sign_state[y][x] & SC_MASK];
      sign_bit = decode_symbol(codec, &sign_models[EZBC_SIGN_CONTEXT(sign_predict)]);
      sign_bit ^= EZBC_SIGN_PREDICTOR(sign_predict);
      update_sc_state(sign_state, x, y, sign_bit);

      if(sign_bit)
	nodes[0][y][x] |= EZBC_SIGN;
    } else {
      /* coefficient is still insignificant */
      *(--LIP_end) = cur_coord;
    }
  }

  codec->LSP_end = LSP_end;
  codec->LIP_end = LIP_end;
  codec->LSP_ids[EZBC_MSB - index][0] = codec->LSP_end;
}


/* decode a significant leaf */
int decode_sig_leaf_coeff(ezbc_subband_codec_t *codec, int x, int y, int index, int sig)
{
  int width = codec->width;
  int height = codec->height;
  int **base_state =  codec->base_state;
  ezbc_coeff_t **coeffs = codec->nodes[0];
  binary_model_t *sig_models = &codec->context_models[codec->context_sig[0].offset];
  ezbc_byte_t *sig_table;

  binary_model_t *sign_models = &codec->context_models[codec->context_sign.offset];
  ezbc_byte_t *sign_table = codec->context_sign.table;
  int **sign_state = codec->sign_state;
  int sign_predict, sign_bit;
  ezbc_coeff_t threshold;
  int **child_node_state = NULL; /* for good measure */
  int **child_base_state = NULL;

  if(codec->child_codec) {
    child_node_state = codec->child_codec->node_state[1];
    child_base_state = codec->child_codec->base_state;
  }

  if(sig)
    sig_table = codec->context_jsig[(x & 1) + ((y & 1) << 1)].table;
  else
    sig_table = codec->context_jsig[3].table;

  /* compare to this threshold to determine signifiance */
  threshold = 1 << index; 

  if(x < width && y < height) {

    /* last child need not be coded if all other children are insignificant
       then it is the one causing significance of its parent */
    if(LAST_NSIG_IN_GROUP ||
       decode_symbol(codec, &sig_models[sig_table[base_state[y][x] & ZC_MASK]])) {
      /* coefficient is significant */
      coeffs[y][x] |= threshold;
      coeffs[y][x] |= EZBC_MSB - index + 1; /* set the length */

      /* update the zero coding state */
      update_zc_state(base_state, x, y);

      if(codec->child_codec) {
	update_child_state(child_node_state, x, y);
	update_child_state(child_base_state, 2*x+0, 2*y+0); /* predict     */
	update_child_state(child_base_state, 2*x+1, 2*y+0); /* children    */
	update_child_state(child_base_state, 2*x+0, 2*y+1); /* are         */
	update_child_state(child_base_state, 2*x+1, 2*y+1); /* significant */
      }

      /* sign coding */
      sign_predict = sign_table[sign_state[y][x] & SC_MASK];
      sign_bit = decode_symbol(codec, &sign_models[EZBC_SIGN_CONTEXT(sign_predict)]);
      sign_bit ^= EZBC_SIGN_PREDICTOR(sign_predict);

      /* update the sign coding context */
      update_sc_state(sign_state, x, y, sign_bit);

      if(sign_bit)
	coeffs[y][x] |= EZBC_SIGN;
      
      /* add to the List of Significant Pixels */
      *++codec->LSP_end = ezbc_point(x, y);

      return(sig | 1);
    } else {
      /* add to the list of Insignificant Pixels */
      *--codec->LIP_end = ezbc_point(x, y);

      return(sig);
    }
  }
  return(0);
}

/* decode leaves of the List of Insignificant Sets */
void decode_LIS_leaves(ezbc_subband_codec_t *codec, int index)
{
  int i;
  int x, y;
  int W, H;
  ezbc_point_t cur_coord, *LIS_cur, *LIS_end, *LIS_end_old;
  ezbc_coeff_t **node;
  int **node_state;
  binary_model_t *sign_models;
  binary_model_t *leaf_models;
  ezbc_byte_t *leaf_table;
  ezbc_coeff_t threshold;
  int **child_node_state = NULL; /* for good measure */
  int sig;

  if(codec->child_codec)
    child_node_state = codec->child_codec->node_state[2];

  /* compare to this threshold to determine signifiance */
  threshold = 1 << index; 

  W = codec->width;
  H = codec->height;

  node_state = codec->node_state[1];
  node = codec->nodes[1];
  leaf_models = &codec->context_models[codec->context_node[1].offset];
  leaf_table = codec->context_node[1].table;

  /* scale sign models */
  sign_models = &codec->context_models[codec->context_sign.offset];
  for(i = 0; i < codec->context_sign.count; i++)
    binary_model_scale(sign_models[i]);

  LIS_cur = LIS_end = codec->LIS[1] - 1;
  LIS_end_old = codec->LIS_end[1];

  while(LIS_cur <  LIS_end_old) {
    cur_coord = *(++LIS_cur);
    y = ezbc_coord_y(cur_coord);
    x = ezbc_coord_x(cur_coord);

    if(decode_symbol(codec, &leaf_models[leaf_table[node_state[y][x] & ZC_MASK]])) {
      /* update the zero coding state for that leaf */
      update_zc_state(node_state, x, y);

      if(codec->child_codec)
	update_child_state(child_node_state, x, y);

      /* decode each coefficient of a significant leaf */
      sig = decode_sig_leaf_coeff(codec, 2*x + 0, 2*y + 0, index, 0);
      sig = decode_sig_leaf_coeff(codec, 2*x + 1, 2*y + 0, index, sig);
      sig = decode_sig_leaf_coeff(codec, 2*x + 0, 2*y + 1, index, sig);
      sig = decode_sig_leaf_coeff(codec, 2*x + 1, 2*y + 1, index, sig);
    } else {
      /* add back to the end of List of Insignificant Sets */
      *++LIS_end = cur_coord;
    }
  }
  codec->LIS_end[1] = LIS_end;
  codec->LSP_ids[EZBC_MSB - index][1] = codec->LSP_end;
}

/* decode a significant node of the quad tree at coord x,y and level l */
int decode_sig_node_child(ezbc_subband_codec_t *codec, int x, int y, int l, int index, int sig)
{
  int w, h;
  int ***node_state = codec->node_state;
  binary_model_t *sig_models = &codec->context_models[codec->context_sig[l].offset];
  ezbc_byte_t *sig_table;
  ezbc_coeff_t threshold;
  int **child_node_state = NULL; /* for good measure */

  if(codec->child_codec)
    child_node_state = codec->child_codec->node_state[l+1];

  if(sig)
    sig_table = codec->context_jsig[4 * l + 3].table;
  else
    sig_table = codec->context_jsig[4 * l + (x & 1) + ((y & 1) << 1)].table;

  /* compare to this threshold to determine signifiance */
  threshold = 1 << index; 

  w = shift_ceil(codec->width, l);
  h = shift_ceil(codec->height, l);

  if(x < w && y < h) {
    /* last child need not be coded if all other children are insignificant
       then it is the one causing significance of its parent */
    if(LAST_NSIG_IN_GROUP ||
       decode_symbol(codec, &sig_models[sig_table[node_state[l][y][x] & ZC_MASK]])) {
      /* update zero coding state */
      update_zc_state(node_state[l], x, y);

     if(codec->child_codec)
	update_child_state(child_node_state, x, y);

      /* add to the stack of insignificant sets */
      codec->LIS_stack_top++;
      codec->LIS_stack_top->x = x;
      codec->LIS_stack_top->y = y;
      codec->LIS_stack_top->l = l;

      return(sig | 1);
    } else {
      /* add back to the list of insignificant nodes */
      *(++codec->LIS_end[l]) = ezbc_point(x, y);

      return(sig);
    }
  }
  return(0);
}

void decode_LIS_stack(ezbc_subband_codec_t *codec, int index)
{
  int x, y, l;
  int W, H;
  int sig;

  while(codec->LIS_stack <= codec->LIS_stack_top)
  {
    x = codec->LIS_stack_top->x;
    y = codec->LIS_stack_top->y;
    l = codec->LIS_stack_top->l; 
    W = shift_ceil(codec->width, l-1);
    H = shift_ceil(codec->height, l-1);
    codec->LIS_stack_top--;

    if(l > 1) {
      sig = decode_sig_node_child(codec, 2*x+0, 2*y+0, l - 1, index, 0);
      sig = decode_sig_node_child(codec, 2*x+1, 2*y+0, l - 1, index, sig);
      sig = decode_sig_node_child(codec, 2*x+0, 2*y+1, l - 1, index, sig);
      sig = decode_sig_node_child(codec, 2*x+1, 2*y+1, l - 1, index, sig);
    } else {
      sig = decode_sig_leaf_coeff(codec, 2*x+0, 2*y+0, index, 0);
      sig = decode_sig_leaf_coeff(codec, 2*x+1, 2*y+0, index, sig);
      sig = decode_sig_leaf_coeff(codec, 2*x+0, 2*y+1, index, sig);
      sig = decode_sig_leaf_coeff(codec, 2*x+1, 2*y+1, index, sig);
    }
  }
}

void decode_qtree_level(ezbc_subband_codec_t *codec, int l, int index)
{
  int i;
  int x, y;
  int W, H;
  ezbc_point_t coord, *LIS_cur, *LIS_end, *LIS_end_old;
  ezbc_coeff_t **nodes;
  int **node_state = codec->node_state[l];
  binary_model_t *node_models = &codec->context_models[codec->context_node[l].offset];
  ezbc_byte_t *node_table = codec->context_node[l].table;
  binary_model_t *sign_models = &codec->context_models[codec->context_sign.offset];
  ezbc_coeff_t threshold;
  int **child_node_state = NULL; /* for good measure */
  int sig;

  if(codec->child_codec)
    child_node_state = codec->child_codec->node_state[l+1];

  /* compare to this threshold to determine signifiance */
  threshold = 1 << index; 

  W = shift_ceil(codec->width, l-1);
  H = shift_ceil(codec->height, l-1);

  /* rescale sign models */
  for(i = 0; i < codec->context_sign.count - 1; i++)
      binary_model_scale(sign_models[i]);

  LIS_cur = LIS_end = codec->LIS[l] - 1;
  LIS_end_old = codec->LIS_end[l];
  nodes = codec->nodes[l];

  while(LIS_cur <  LIS_end_old){
    coord = *(++LIS_cur);
    y = ezbc_coord_y(coord);
    x = ezbc_coord_x(coord);

    if(decode_symbol(codec, &node_models[node_table[node_state[y][x] & ZC_MASK]])) {
      update_zc_state(node_state, x, y);

      if(codec->child_codec)
	update_child_state(child_node_state, x, y);

      /* decode the significant node */
      sig = decode_sig_node_child(codec, 2*x+0, 2*y+0, l - 1, index, 0);
      sig = decode_sig_node_child(codec, 2*x+1, 2*y+0, l - 1, index, sig);
      sig = decode_sig_node_child(codec, 2*x+0, 2*y+1, l - 1, index, sig);
      sig = decode_sig_node_child(codec, 2*x+1, 2*y+1, l - 1, index, sig);

      decode_LIS_stack(codec, index);

    } else {
      /* leave the node in the List of Insignificant Sets */
      *++LIS_end = coord;
    }
  }

  codec->LIS_end[l] = LIS_end;
  codec->LSP_ids[EZBC_MSB - index][l] = codec->LSP_end;
}
void decode_LSP_and_bit_index(ezbc_subband_codec_t *codec, int index)
{
  ezbc_coord_t x, y;
  ezbc_point_t coord;
  ezbc_point_t *last, *sp;
  int LSP_plane = codec->LSP_plane;
  ezbc_point_t *LSP_mark = codec->LSP_mark;

  ezbc_point_t **LSP_bit_index_marks = codec->LSP_bit_index_marks;
  ezbc_point_t ***LSP_ids = codec->LSP_ids;
  int i;
  int depth = codec->depth;

  binary_model_t *LSP_models = &codec->context_models[codec->context_LSP.offset];
  ezbc_byte_t *LSP_table = codec->context_LSP.table;
  int **base_state = codec->base_state;
  ezbc_coeff_t **coeffs = codec->nodes[0];
  ezbc_coeff_t mag;
  int bit;
  int LSP_level;
  int LSP_offset;

  LSP_mark[LSP_plane] = codec->LSP_end - codec->LSP;
  LSP_bit_index_marks[EZBC_MSB - index] = codec->LSP_end;
  
  if(index < EZBC_MSB) {
    
    sp = LSP_bit_index_marks[EZBC_MSB - index - 1];

    for(LSP_level = depth - 1; LSP_level > 1; LSP_level--) {
      last = LSP_ids[EZBC_MSB - index - 1][LSP_level - 1];
      for(; sp > last; sp--) {
	coord = *sp;
	y = ezbc_coord_y(coord);
	x = ezbc_coord_x(coord);
	mag = (coeffs[y][x] & EZBC_ABS) >> index;
	
	LSP_offset = LSP_table[base_state[y][x]];
	
	if(!vlc_is_leaf(codec, coeffs[y][x])) {
	  /* decode refinement bit */
	  bit = decode_symbol(codec, &LSP_models[LSP_offset]);
	  coeffs[y][x] |= bit << index;
	  coeffs[y][x]++; /* increase the length */
	}
      }

      /* scale LSP models */
      for(i = 0; i < codec->context_LSP.count; i++)
	binary_model_scale(LSP_models[i]);
    }

    if(index < EZBC_MSB - 1) {
      /* from leaves ?? */
      last = LSP_ids[EZBC_MSB - index - 1][0];

      for(; sp > last; sp--) {

	coord = *sp;
	y = ezbc_coord_y(coord);
	x = ezbc_coord_x(coord);
	mag = (coeffs[y][x] & EZBC_ABS) >> index;
	LSP_offset = LSP_table[base_state[y][x]];

	if(!vlc_is_leaf(codec, coeffs[y][x])) {
	  /* decode refinement bit */
	  bit = decode_symbol(codec, &LSP_models[LSP_offset]);
	  coeffs[y][x] |= bit << index;
	  coeffs[y][x]++; /* increase the length */
	}
      }

      /* scale LSP models */
      for(i = 0; i < codec->context_LSP.count; i++)
	binary_model_scale(LSP_models[i]);

      last = LSP_bit_index_marks[EZBC_MSB - index - 2];
      for(; sp > last; sp--) {
	coord = *sp;
	y = ezbc_coord_y(coord);
	x = ezbc_coord_x(coord);
	mag = (coeffs[y][x] & EZBC_ABS) >> index;
	
	LSP_offset = LSP_table[base_state[y][x]];

	if(!vlc_is_leaf(codec, coeffs[y][x])) {
	  /* decode refinement bit */
	  bit = decode_symbol(codec, &LSP_models[LSP_offset]);
	  coeffs[y][x] |= bit << index;
	  coeffs[y][x]++; /* increase the length */
	}
      }
    }
  }

  if(index < EZBC_MSB - 1) {

    /* scale LSP models */
    for(i = 0; i < codec->context_LSP.count; i++)
      binary_model_scale(LSP_models[i]);

    for(LSP_level = depth - 1; LSP_level >= 0; LSP_level--) { 
      if(LSP_level)
	last = LSP_ids[EZBC_MSB - index - 2][LSP_level - 1];
      else {
	if(index == EZBC_MSB - 2)
	  last = codec->LSP - 1;
	else
	  last = LSP_bit_index_marks[EZBC_MSB - index - 3];
      }
	    
      for(; sp > last; sp--) {
	coord = *sp;
	y = ezbc_coord_y(coord);
	x = ezbc_coord_x(coord);
	mag = (coeffs[y][x] & EZBC_ABS) >> index;
	
	LSP_offset = LSP_table[base_state[y][x]];

	if(!vlc_is_leaf(codec, coeffs[y][x])) {
	  /* decode refinement bit */
	  bit = decode_symbol(codec, &LSP_models[LSP_offset]);
	  coeffs[y][x] |= bit << index;
	  coeffs[y][x]++; /* increase the length */
	}
      }

      /* scale LSP models */
      for(i = 0; i < codec->context_LSP.count; i++)
	binary_model_scale(LSP_models[i]);
    }
  }

  if(index < EZBC_MSB - 2) {
    
    /* scale LSP models */
    for(i = 0; i < codec->context_LSP.count; i++)
      binary_model_scale(LSP_models[i]);

    last = codec->LSP - 1;
    for(; sp > last; sp--){
	coord = *sp;
	y = ezbc_coord_y(coord);
	x = ezbc_coord_x(coord);
	mag = (coeffs[y][x] & EZBC_ABS) >> index;
	
	LSP_offset = LSP_table[base_state[y][x]];

	if(!vlc_is_leaf(codec, coeffs[y][x])) {
	  /* decode refinement bit */
	  bit = decode_symbol(codec, &LSP_models[LSP_offset]);
	  coeffs[y][x] |= bit << index;
	  coeffs[y][x]++; /* increase the length */
	}
    }
  }

  codec->LSP_plane++;
  codec->LSP_split = codec->LSP;
}

void ezbc_stop_decoding(ezbc_subband_codec_t *codec)
{
  int x, y;
  int length;
  int LSP_plane = codec->LSP_plane;
  ezbc_point_t coord;
  ezbc_point_t *sp;
  ezbc_point_t *LSP = codec->LSP;
  ezbc_point_t *LSP_end = codec->LSP_end;
  ezbc_point_t *LSP_split = codec->LSP_split;
  ezbc_coeff_t **coeffs = codec->nodes[0];
  ezbc_point_t *LSP_mark = codec->LSP_mark;

  LSP_plane --;

  if(LSP_plane >= 0)
    LSP_split = LSP + LSP_mark[LSP_plane] + 1;
  else
    return;

  for(sp = LSP; sp < LSP_split; sp++) {
      coord = *sp;
      y = ezbc_coord_y(coord);
      x = ezbc_coord_x(coord);
      length = coeffs[y][x] & EZBC_LENGTH;
      if(length && length < LSP_plane) continue;
      coeffs[y][x] &= ~EZBC_LENGTH;
      coeffs[y][x] |= LSP_plane;
  }

  LSP_plane ++;

  for(; sp <= LSP_end; sp++) {
      coord = *sp;
      y = ezbc_coord_y(coord);
      x = ezbc_coord_x(coord);
      length = coeffs[y][x] & EZBC_LENGTH;
      if(length && length < LSP_plane) continue;
      coeffs[y][x] &= ~EZBC_LENGTH;
      coeffs[y][x] |= LSP_plane;
  }
}

static int decoded_length(ezbc_subband_codec_t *codec)
{
  return(8 * bitbuffer_length(codec->arithmetic_codec->bitbuffer));
}

static int marker(ezbc_subband_codec_t *codec)
{
  return(arithmetic_codec_decode(codec->arithmetic_codec, 0.5));
}

void ezbc_subband_decode(ezbc_subband_codec_t *codec, int bit)
{
  int depth = codec->depth;
  int W, H, i, l;
  int sig;
  int target_rate = codec->target_rate;

  if(!marker(codec) || decoded_length(codec) >= target_rate) return;

  update_node_contexts(codec);

  /* decode List Of Insignificant Pixel */
  codec->LSP_split = codec->LSP_end + 1;
  decode_LIP(codec, bit);

  if(!marker(codec) || decoded_length(codec) >= target_rate) return;

  /* decode leaves of the signifiance quad tree */
  codec->LSP_split = codec->LSP_end + 1;
  decode_LIS_leaves(codec, bit);

  if(!marker(codec) || decoded_length(codec) >= target_rate) return;

  /* decode the signifiance quad tree */
  codec->LSP_split = codec->LSP_end + 1;
  for(l = 1; l < depth; l++) {
    if(bit == EZBC_MSB && l == depth - 1) {

      /* empty lists in the top plane */
      for(i = 0; i < depth; i++)
	codec->LSP_ids[EZBC_MSB - bit][i] = codec->LSP_end;

      /* decode root */
      W = shift_ceil(codec->width, l-1);
      H = shift_ceil(codec->height, l-1);

      sig = decode_sig_node_child(codec, 0, 0, l - 1, bit, 0);
      sig = decode_sig_node_child(codec, 0, 1, l - 1, bit, sig);
      sig = decode_sig_node_child(codec, 1, 0, l - 1, bit, sig);
      /* make sure the last coefficient of the group is coded */
      sig = 1;
      sig = decode_sig_node_child(codec, 1, 1, l - 1, bit, sig);
      decode_LIS_stack(codec, bit);
    } else {
      decode_qtree_level(codec, l, bit);
    }
  }

  if(!marker(codec) || decoded_length(codec) >= target_rate) return;

  /* decode List of Significant Pixels */
  decode_LSP_and_bit_index(codec, bit);

  /* reset context models */
  for(l = 0; l < codec->context_count; l++)
    binary_model_scale(codec->context_models[l]);
}

/* initialize vlc decoding */
void ezbc_init_vlc_decoding(ezbc_subband_codec_t *codec,
			    int *binbook,
			    int *hash)
{
  codec->binbook = binbook;
  codec->hash = hash;
 }

int vlc_is_leaf(ezbc_subband_codec_t *codec, int code)
{
  int s, e, c, v, h, l;
  int *binbook = codec->binbook;
  int *hash = codec->hash;

  /* flc case */
  if(hash == NULL || binbook == NULL) return(0);

  /* get the current length of the code */
  l = code & EZBC_LENGTH;
  v = code & EZBC_ABS;

  /* search the sorted binarization book for 'code' */
  s = hash[l];   /* start */
  e = hash[l+1]; /* end   */
  if(s == e) return(0); /* no leaf of that size exist */

  while(e > s + 1) {
    /* compare with the median of the range */
    h = (s + e) / 2; 
    c = binbook[h] & EZBC_ABS;

    /* narrow down the range */
    if(v < c)
      e = h;
    else
      s = h;
  }

  /* find the closest vector */
  if(v == (binbook[s] & EZBC_ABS)) return(1);
  if(v == (binbook[e] & EZBC_ABS)) return(1);

  return(0);
}
