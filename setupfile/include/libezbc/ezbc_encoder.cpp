/* libezbc - EZBC subband coding/decoding
      Copyright (C) 2004 Vivien Chappelier

      This library is inspired from code originaly developped by Yongjun Wu Jan in MC-EZBC. This rewritten version may still contain code similar to his original implementation.

      Reference paper: "Embedded Image Coding Using Zeroblocks Of Subband/Wavelet Coefficient And Context Modeling", Shih-Ta Hsiang and John W. Woods, Data Compression Conference (DCC '01), March 27 - 29, 2001, Snowbird, Utah 
*/

#include <math.h>
#include <stdio.h>
#include "ezbc_types.h"
#include "ezbc_tables.h"
#include "ezbc_codec.h"

#define qtree_C(l, x, y) ((1 << (2*l)) | (y << l) | (x))

#define log2(x) (log(x)/log(2.))

int bitcount_sign = 0;
int bitcount_LSP = 0;
int bitcount_LIP = 0;
int bitcount_leaf = 0;
int bitcount_node = 0;
int bitcount_jsig = 0;

float cost_sign = 0;
float cost_LSP = 0;
float cost_LIP = 0;
float cost_leaf = 0;
float cost_node = 0;
float cost_jsig = 0;

void dump_ezbc_cost()
{
  /* dump costs */
  fprintf(stderr, "sign : count = %d cost = %f (%2.2f %%)\n", bitcount_sign, cost_sign, 100 * cost_sign / bitcount_sign);
  fprintf(stderr, "LSP : count = %d cost = %f (%2.2f %%)\n", bitcount_LSP, cost_LSP, 100 * cost_LSP / bitcount_LSP);
  fprintf(stderr, "LIP : count = %d cost = %f (%2.2f %%)\n", bitcount_LIP, cost_LIP, 100 * cost_LIP / bitcount_LIP);
  fprintf(stderr, "leaf : count = %d cost = %f (%2.2f %%)\n", bitcount_leaf, cost_leaf, 100 * cost_leaf / bitcount_leaf);
  fprintf(stderr, "node : count = %d cost = %f (%2.2f %%)\n", bitcount_node, cost_node, 100 * cost_node / bitcount_node);
  fprintf(stderr, "jsig : count = %d cost = %f (%2.2f %%)\n", bitcount_jsig, cost_jsig, 100 * cost_jsig / bitcount_jsig);
  fprintf(stderr, "total : count = %d cost = %f (%d bytes)\n",
	  bitcount_sign + bitcount_LSP + bitcount_LIP + bitcount_leaf + bitcount_node + bitcount_jsig,
	  cost_sign + cost_LSP + cost_LIP + cost_leaf + cost_node + cost_jsig,
	  (int) ((cost_sign + cost_LSP + cost_LIP + cost_leaf + cost_node + cost_jsig) / 8.0));
}

static float cost(int bit, binary_model_t *model)
{
  double entropy;

  if(bit)
    entropy = -log2(binary_model_p1(*model));
  else
    entropy = -log2(binary_model_p0(*model));

  return(entropy);
}

/* encode a bit using a binary arithmetic coder */
void encode_symbol(ezbc_subband_codec_t *codec, int bit, binary_model_t *model)
{
  double entropy;
  arithmetic_codec_t *arithmetic_codec = codec->arithmetic_codec;

  if(bit)
    entropy = -log2(binary_model_p1(*model));
  else
    entropy = -log2(binary_model_p0(*model));

  arithmetic_codec_encode(arithmetic_codec, binary_model_p0(*model), bit);

  binary_model_update(*model, bit);
}

/* encode the List of Insignificant Pixels */
void encode_LIP(ezbc_subband_codec_t *codec, int index)
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
    if(coeff & threshold) {
      /* coefficient becomes significant */
      if(++LSP_end == LIP_old_end) /* LSP_pos runs beyond the end of the LSP */
        *(LIP_cur++) = *(LIP_old_end++);
      *LSP_end = cur_coord;

      // TEMP
      bitcount_LIP++;
      cost_LIP += cost(1, &LIP_models[LIP_context_table[base_state[y][x] & ZC_MASK]]);


      /* encode a 1 */
      encode_symbol(codec, 1, &LIP_models[LIP_context_table[base_state[y][x] & ZC_MASK]]);

      /* update contexts */
      update_zc_state(base_state, x, y);

      if(codec->child_codec) {
	update_child_state(child_node_state, x, y);
	update_child_state(child_base_state, 2*x+0, 2*y+0); /* predict     */
	update_child_state(child_base_state, 2*x+1, 2*y+0); /* children    */
	update_child_state(child_base_state, 2*x+0, 2*y+1); /* are         */
	update_child_state(child_base_state, 2*x+1, 2*y+1); /* significant */
      }

      /* encode the sign */
      sign_predict = sign_context_table[sign_state[y][x] & SC_MASK];
      sign_bit = (coeff & EZBC_SIGN)?1:0;
      sign_bit ^= EZBC_SIGN_PREDICTOR(sign_predict);

      // TEMP
      bitcount_sign++;
      cost_sign += cost(sign_bit, &sign_models[EZBC_SIGN_CONTEXT(sign_predict)]);

      encode_symbol(codec, sign_bit, &sign_models[EZBC_SIGN_CONTEXT(sign_predict)]);

      update_sc_state(sign_state, x, y, (coeff & EZBC_SIGN)?1:0);
    } else {
      /* coefficient is still insignificant */
      *(--LIP_end) = cur_coord;

      // TEMP
      bitcount_LIP++;
      cost_LIP += cost(0, &LIP_models[LIP_context_table[base_state[y][x] & ZC_MASK]]);

      /* encode a 0 */
      encode_symbol(codec, 0, &LIP_models[LIP_context_table[base_state[y][x] & ZC_MASK]]);


    }
  }

  codec->LSP_end = LSP_end;
  codec->LIP_end = LIP_end;
  codec->LSP_ids[EZBC_MSB - index][0] = codec->LSP_end;
}

/* encode a significant leaf */
/* return 1 if leaf is significant */
/* sig=1 if at least one previous coefficient in the group is significant */
int encode_sig_leaf_coeff(ezbc_subband_codec_t *codec, int x, int y, int index, int sig)
{
  int width = codec->width;
  int height = codec->height;
  ezbc_coeff_t coeff;
  int **base_state =  codec->base_state;
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
    coeff = codec->nodes[0][y][x];
    if(coeff & threshold) {

      /* last child need not be coded if all other children are insignificant
	 then it is the one causing significance of its parent */
      if(!LAST_NSIG_IN_GROUP) {
	// TEMP
	bitcount_jsig++;
	cost_jsig += cost(1, &sig_models[sig_table[base_state[y][x] & ZC_MASK]]);

	/* coefficient is significant */
	encode_symbol(codec, 1, &sig_models[sig_table[base_state[y][x] & ZC_MASK]]);
      }

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
      sign_bit = (coeff & EZBC_SIGN)?1:0;
      sign_bit ^= EZBC_SIGN_PREDICTOR(sign_predict);

      // TEMP
      bitcount_sign++;
      cost_sign += cost(sign_bit, &sign_models[EZBC_SIGN_CONTEXT(sign_predict)]);
      
      encode_symbol(codec, sign_bit, &sign_models[EZBC_SIGN_CONTEXT(sign_predict)]);

      /* update the sign coding context */
      update_sc_state(sign_state, x, y, (coeff & EZBC_SIGN)?1:0);
      

      /* add to the List of Significant Pixels */
      *++codec->LSP_end = ezbc_point(x, y);

      return(sig | 1);
    } else {
      // TEMP
      bitcount_jsig++;
      cost_jsig += cost(0, &sig_models[sig_table[base_state[y][x] & ZC_MASK]]);

      /* coefficient is insignificant */
      encode_symbol(codec, 0, &sig_models[sig_table[base_state[y][x] & ZC_MASK]]);

      /* add to the list of Insignificant Pixels */
      *--codec->LIP_end = ezbc_point(x, y);

      return(sig);
    }
  }

  return(0);
}

/* encode leaves of the List of Insignificant Sets */
void encode_LIS_leaves(ezbc_subband_codec_t *codec, int index)
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

  W = codec->width;
  H = codec->height;

  while(LIS_cur <  LIS_end_old) {
    cur_coord = *(++LIS_cur);
    y = ezbc_coord_y(cur_coord);
    x = ezbc_coord_x(cur_coord);

    if(node[y][x] & threshold) {

      // TEMP
      bitcount_leaf++;
      cost_leaf += cost(1, &leaf_models[leaf_table[node_state[y][x] & ZC_MASK]]);
      /* leaf is significant */
      encode_symbol(codec, 1, &leaf_models[leaf_table[node_state[y][x] & ZC_MASK]]);

      /* update the zero coding state for that leaf */
      update_zc_state(node_state, x, y);

      if(codec->child_codec)
	update_child_state(child_node_state, x, y);

      /* encode each coefficient of a significant leaf */
      sig = encode_sig_leaf_coeff(codec, 2*x + 0, 2*y + 0, index, 0);
      sig = encode_sig_leaf_coeff(codec, 2*x + 1, 2*y + 0, index, sig);
      sig = encode_sig_leaf_coeff(codec, 2*x + 0, 2*y + 1, index, sig);
      sig = encode_sig_leaf_coeff(codec, 2*x + 1, 2*y + 1, index, sig);
    } else {
      // TEMP
      bitcount_leaf++;
      cost_leaf += cost(0, &leaf_models[leaf_table[node_state[y][x] & ZC_MASK]]);

      /* encode a zero */
      encode_symbol(codec, 0, &leaf_models[leaf_table[node_state[y][x] & ZC_MASK]]);

      /* add back to the end of List of Insignificant Sets */
      *++LIS_end = cur_coord;
    }
  }
  codec->LIS_end[1] = LIS_end;
  codec->LSP_ids[EZBC_MSB - index][1] = codec->LSP_end;
}

/* encode a significant node of the quad tree at coord x,y and level l */
int encode_sig_node_child(ezbc_subband_codec_t *codec, int x, int y, int l, int index, int sig)
{
  int w, h;
  ezbc_coeff_t ***nodes = codec->nodes;
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
    if(nodes[l][y][x] & threshold) {

      /* last child need not be coded if all other children are insignificant
	 then it is the one causing significance of its parent */
      if(!LAST_NSIG_IN_GROUP) {
	// TEMP
	bitcount_jsig++;
	cost_jsig += cost(1, &sig_models[sig_table[node_state[l][y][x] & ZC_MASK]]);

	/* node is significant */
	encode_symbol(codec, 1, &sig_models[sig_table[node_state[l][y][x] & ZC_MASK]]);
      }

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
      // TEMP
      bitcount_jsig++;
      cost_jsig += cost(0, &sig_models[sig_table[node_state[l][y][x] & ZC_MASK]]);

      /* node is insignificant */
      encode_symbol(codec, 0, &sig_models[sig_table[node_state[l][y][x] & ZC_MASK]]);

      /* add back to the list of insignificant nodes */
      *(++codec->LIS_end[l]) = ezbc_point(x, y);

      return(sig);
    }
  }
  return(0);
}

void encode_LIS_stack(ezbc_subband_codec_t *codec, int index)
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
      sig = encode_sig_node_child(codec, 2*x+0, 2*y+0, l - 1, index, 0);
      sig = encode_sig_node_child(codec, 2*x+1, 2*y+0, l - 1, index, sig);
      sig = encode_sig_node_child(codec, 2*x+0, 2*y+1, l - 1, index, sig);
      sig = encode_sig_node_child(codec, 2*x+1, 2*y+1, l - 1, index, sig);
    } else {
      sig = encode_sig_leaf_coeff(codec, 2*x+0, 2*y+0, index, 0);
      sig = encode_sig_leaf_coeff(codec, 2*x+1, 2*y+0, index, sig);
      sig = encode_sig_leaf_coeff(codec, 2*x+0, 2*y+1, index, sig);
      sig = encode_sig_leaf_coeff(codec, 2*x+1, 2*y+1, index, sig);
    }
  }

}

void encode_qtree_level(ezbc_subband_codec_t *codec, int l, int index)
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

  /* rescale sign models */
  for(i = 0; i < codec->context_sign.count - 1; i++)
      binary_model_scale(sign_models[i]);

  LIS_cur = LIS_end = codec->LIS[l] - 1;
  LIS_end_old = codec->LIS_end[l];
  nodes = codec->nodes[l];
  W = shift_ceil(codec->width, l-1);
  H = shift_ceil(codec->height, l-1);

  while(LIS_cur <  LIS_end_old){
    coord = *(++LIS_cur);
    y = ezbc_coord_y(coord);
    x = ezbc_coord_x(coord);

    if(nodes[y][x] & threshold){

      // TEMP
      bitcount_node++;
      cost_node += cost(1, &node_models[node_table[node_state[y][x] & ZC_MASK]]);

      /* node becomes significant */
      encode_symbol(codec, 1, &node_models[node_table[node_state[y][x] & ZC_MASK]]);

      update_zc_state(node_state, x, y);

      if(codec->child_codec)
	update_child_state(child_node_state, x, y);

      /* encode the significant node */
      sig = encode_sig_node_child(codec, 2*x+0, 2*y+0, l - 1, index, 0);
      sig = encode_sig_node_child(codec, 2*x+1, 2*y+0, l - 1, index, sig);
      sig = encode_sig_node_child(codec, 2*x+0, 2*y+1, l - 1, index, sig);
      sig = encode_sig_node_child(codec, 2*x+1, 2*y+1, l - 1, index, sig);

      encode_LIS_stack(codec, index);

    } else {
      // TEMP
      bitcount_node++;
      cost_node += cost(0, &node_models[node_table[node_state[y][x] & ZC_MASK]]);

      /* node is insignificant */
      encode_symbol(codec, 0, &node_models[node_table[node_state[y][x] & ZC_MASK]]);

      /* leave the node in the List of Insignificant Sets */
      *++LIS_end = coord;
    }
  }

  codec->LIS_end[l] = LIS_end;
  codec->LSP_ids[EZBC_MSB - index][l] = codec->LSP_end;
}

void encode_LSP_and_bit_index(ezbc_subband_codec_t *codec, int index)
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
  int length;
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

	length = coeffs[y][x] & EZBC_LENGTH;

	if(!length || length >= EZBC_MSB + 1 - index) {
	  bitcount_LSP++;
	  cost_LSP += cost(mag & 1, &LSP_models[LSP_offset]);
	  encode_symbol(codec, mag & 1, &LSP_models[LSP_offset]);
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

	length = coeffs[y][x] & EZBC_LENGTH;

	if(!length || length >= EZBC_MSB + 1 - index) {
	  bitcount_LSP++;
	  cost_LSP += cost(mag & 1, &LSP_models[LSP_offset]);
	  encode_symbol(codec, mag & 1, &LSP_models[LSP_offset]);
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

	length = coeffs[y][x] & EZBC_LENGTH;

	if(!length || length >= EZBC_MSB + 1 - index) {
	  bitcount_LSP++;
	  cost_LSP += cost(mag & 1, &LSP_models[LSP_offset]);
	  encode_symbol(codec, mag & 1, &LSP_models[LSP_offset]);
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
	
	length = coeffs[y][x] & EZBC_LENGTH;

	if(!length || length >= EZBC_MSB + 1 - index) {
	  bitcount_LSP++;
	  cost_LSP += cost(mag & 1, &LSP_models[LSP_offset]);
	  encode_symbol(codec, mag & 1, &LSP_models[LSP_offset]);
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
    for(; sp > last; sp--) {
      coord = *sp;
      y = ezbc_coord_y(coord);
      x = ezbc_coord_x(coord);
      mag = (coeffs[y][x] & EZBC_ABS) >> index;
      
      LSP_offset = LSP_table[base_state[y][x]];

      length = coeffs[y][x] & EZBC_LENGTH;      

      if(!length || length >= EZBC_MSB + 1 - index) {
	bitcount_LSP++;
	cost_LSP += cost(mag & 1, &LSP_models[LSP_offset]);
	encode_symbol(codec, mag & 1, &LSP_models[LSP_offset]);
      }
    }
  }

  codec->LSP_plane++;
}

static int encoded_length(ezbc_subband_codec_t *codec)
{
  return(8 * bitbuffer_length(codec->arithmetic_codec->bitbuffer));
}

static void marker(ezbc_subband_codec_t *codec, int bit)
{
  arithmetic_codec_encode(codec->arithmetic_codec, 0.5, bit);
}

int ezbc_subband_encode(ezbc_subband_codec_t *codec, int bit, int *truncation_list)
{
  int depth = codec->depth;
 int W, H, i, l;
  int target_rate = codec->target_rate;
  int sig;
  int n = 0;

  if(encoded_length(codec) >= target_rate) {
    marker(codec, 0);
    return(n);
  } else
    marker(codec, 1);

  update_node_contexts(codec);

  /* encode List Of Insignificant Pixel */
  encode_LIP(codec, bit);
  if(truncation_list) truncation_list[n++] = encoded_length(codec);
  if(encoded_length(codec) >= target_rate) {
    marker(codec, 0);
    return(n);
  } else
    marker(codec, 1);


  /* encode leaves of the signifiance quad tree */
  encode_LIS_leaves(codec, bit);
  if(truncation_list) truncation_list[n++] = encoded_length(codec);
  if(encoded_length(codec) >= target_rate) {
    marker(codec, 0);
    return(n);
  } else
    marker(codec, 1);

  /* encode signifiance quad tree */
  for(l = 1; l < depth; l++) {

    if(bit == EZBC_MSB && l == depth - 1) {

      /* empty lists in the top plane */
      for(i = 0; i < depth; i++)
	codec->LSP_ids[EZBC_MSB - bit][i] = codec->LSP_end;

      /* encode root */
      W = shift_ceil(codec->width, l-1);
      H = shift_ceil(codec->height, l-1);

      // TODO
      sig = encode_sig_node_child(codec, 0, 0, l - 1, bit, 0);
      sig = encode_sig_node_child(codec, 0, 1, l - 1, bit, sig);
      sig = encode_sig_node_child(codec, 1, 0, l - 1, bit, sig);
      /* make sure the last coefficient of the group is coded */
      sig = 1;
      sig = encode_sig_node_child(codec, 1, 1, l - 1, bit, sig);

      encode_LIS_stack(codec, bit);
    } else {
      encode_qtree_level(codec, l, bit);
    }
  }

  if(truncation_list) truncation_list[n++] = encoded_length(codec);
  if(encoded_length(codec) >= target_rate) {
    marker(codec, 0);
    return(n);
  } else
    marker(codec, 1);

  /* encode List of Significant Pixels */
  encode_LSP_and_bit_index(codec, bit);

  /* update context models */
  for(l = 0; l < codec->context_count; l++)
    binary_model_scale(codec->context_models[l]);

  if(truncation_list) truncation_list[n++] = encoded_length(codec);
  return(n);
}
