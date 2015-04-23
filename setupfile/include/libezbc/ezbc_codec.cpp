/* libezbc - EZBC subband coding/decoding
      Copyright (C) 2004 Vivien Chappelier

      This library is inspired from code originaly developped by Yongjun Wu Jan in MC-EZBC. This rewritten version may still contain code similar to his original implementation.

      Reference paper: "Embedded Image Coding Using Zeroblocks Of Subband/Wavelet Coefficient And Context Modeling", Shih-Ta Hsiang and John W. Woods, Data Compression Conference (DCC '01), March 27 - 29, 2001, Snowbird, Utah 
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "ezbc_types.h"
#include "ezbc_tables.h"
#include "ezbc_codec.h"

void update_node_contexts(ezbc_subband_codec_t *codec)
{
  ezbc_point_t *list_cur, *list_end, coord_cur;
  int *base_cxt_sp, **base_cxt;
  int *context_node_sp, **context_node;
  int l;

  /* update list of insignificant pixels */
  list_cur = codec->LIP;
  list_end = codec->LIP_end;
  base_cxt = codec->base_state;

  while(list_cur > list_end){
    coord_cur = *(--list_cur);
    base_cxt_sp = &base_cxt[ezbc_coord_y(coord_cur)][ezbc_coord_x(coord_cur)];

    if(*base_cxt_sp & TC_SIG)
      *base_cxt_sp |= TC2_SIG; // TODO: ??

    if(*base_cxt_sp & CL_SIG)
      *base_cxt_sp |= CL2_SIG;

    if(*base_cxt_sp & CR_SIG)
      *base_cxt_sp |= CR2_SIG;

    if(*base_cxt_sp & BC_SIG)
      *base_cxt_sp |= BC2_SIG;
  }

  /* update list of insignificant sets */
  int depth = codec->depth;

  for(l = 1; l < depth; l++) {
    list_cur = codec->LIS[l];
    list_end = codec->LIS_end[l];
    context_node = codec->node_state[l];
    while(list_cur < list_end) {
     coord_cur = *(++list_cur);
     context_node_sp = context_node[ezbc_coord_y(coord_cur)] + ezbc_coord_x(coord_cur);

     if(*context_node_sp & TC_SIG)
       *context_node_sp |= TC2_SIG;

     if(*context_node_sp & CL_SIG)
       *context_node_sp |= CL2_SIG;

     if(*context_node_sp & CR_SIG)
       *context_node_sp |= CR2_SIG;

     if(*context_node_sp & BC_SIG)
       *context_node_sp |= BC2_SIG;
    }
  }

  /* update list of significant pixels */
  list_cur = codec->LSP;
  list_end = codec->LSP_end;

  while(list_cur < list_end){
    coord_cur = *(++list_cur);
    base_cxt_sp = &base_cxt[ezbc_coord_y(coord_cur)][ezbc_coord_x(coord_cur)];

    if(*base_cxt_sp & TC_SIG)
      *base_cxt_sp |= TC2_SIG;

    if(*base_cxt_sp & CL_SIG)
      *base_cxt_sp |= CL2_SIG;

    if(*base_cxt_sp & CR_SIG)
      *base_cxt_sp |= CR2_SIG;

    if(*base_cxt_sp & BC_SIG)
      *base_cxt_sp |= BC2_SIG;
  }
}

/* initialize zero coding, sign coding and quad tree state */
static void init_state(ezbc_subband_codec_t *codec)
{
  int ***node_state = codec->node_state;
  int **base_state = codec->base_state;
  int **sign_state = codec->sign_state;
  int depth = codec->depth;
  int x, y, l;
  int width = codec->width;
  int height = codec->height;

  /* pixel state initialization */
  for(y = -EZBC_CONTEXT_SIZE; y < height + EZBC_CONTEXT_SIZE; y++)
    for(x = -EZBC_CONTEXT_SIZE; x < width + EZBC_CONTEXT_SIZE; x++)
      if(x >= 0 && x < width && y >= 0 && y < height)
	base_state[y][x] = 0; /* initialize context to 0 inside the frame */
      else
	base_state[y][x] = OUT_OF_BOUNDS; /* initialize context to 'out of bounds' outside */

  /* node state initialization */
  for(l = 0; l < depth; l++) {
    int w, h;

    w = (( width + (1 << l) - 1) >> l);
    h = ((height + (1 << l) - 1) >> l);

    for(y = -EZBC_CONTEXT_SIZE; y < h + EZBC_CONTEXT_SIZE; y++)
      for(x = -EZBC_CONTEXT_SIZE; x < w + EZBC_CONTEXT_SIZE; x++)
	if(x >= 0 && x < w && y >= 0 && y < h)	
	  node_state[l][y][x] = 0;
	else
	  node_state[l][y][x] = OUT_OF_BOUNDS;
  }

  /* initialize sign state */
  for(y = -EZBC_CONTEXT_SIZE; y < height + EZBC_CONTEXT_SIZE; y++)
    for(x = -EZBC_CONTEXT_SIZE; x < width + EZBC_CONTEXT_SIZE; x++)
      sign_state[y][x] = 0;
}

/* initialize the subband quad tree with the maximum absolute values */
static void init_node_qtree(ezbc_subband_codec_t *codec)
{
  ezbc_coeff_t ***nodes = codec->nodes;
  int depth = codec->depth;
  int x, y, p, l;
  int width = codec->width;
  int height = codec->height;
  ezbc_coeff_t *coeffs;

  coeffs = codec->coeffs;

  p = width;

  /* initialize node values */
  /* WARNING: for decoding, the coefficients must be zeroed beforehand */
  for(l = 1; l < depth; l++) {
    int w, h;
    int W, H;

    w = shift_ceil(width, l);
    h = shift_ceil(height, l);
    W = shift_ceil(width, l-1);
    H = shift_ceil(height, l-1);

    for(y = 0; y < h; y++) {
      for(x = 0; x < w; x++) {
	
	/* propagate up the maximum absolute value */
        nodes[l][y][x]  = nodes[l-1][2*y+0][2*x+0];
	if(2*x+1 < W) nodes[l][y][x] |= nodes[l-1][2*y+0][2*x+1];
	if(2*y+1 < H) nodes[l][y][x] |= nodes[l-1][2*y+1][2*x+0];
	if(2*x+1 < W && 2*y+1 < H) nodes[l][y][x] |= nodes[l-1][2*y+1][2*x+1];
      }
    }
  }
}

/* update the zero coding state of neighbors */
void update_zc_state(int **base_state, int x, int y)
{
  base_state[y - 1][x - 1] |= BR_SIG;
  base_state[y - 1][x    ] |= BC_SIG;
  base_state[y - 1][x + 1] |= BL_SIG;
  base_state[y    ][x - 1] |= CR_SIG;
  base_state[y    ][x + 1] |= CL_SIG;
  base_state[y + 1][x - 1] |= TR_SIG;
  base_state[y + 1][x    ] |= TC_SIG;
  base_state[y + 1][x + 1] |= TL_SIG;
}

/* update the sign coding state of neighbors */
void update_sc_state(int **sign_state, int x, int y, int sign)
{
  if(sign) {
    sign_state[y - 1][x - 1] |= NW_NVE_MASK;
    sign_state[y - 1][x    ] |= V_NVE_MASK;
    sign_state[y - 1][x + 1] |= NE_NVE_MASK;
    sign_state[y    ][x - 1] |= H_NVE_MASK;
    sign_state[y    ][x + 1] |= H_NVE_MASK;
    sign_state[y + 1][x - 1] |= NE_NVE_MASK;
    sign_state[y + 1][x    ] |= V_NVE_MASK;
    sign_state[y + 1][x + 1] |= NW_NVE_MASK;
  } else {
    sign_state[y - 1][x - 1] |= NW_PVE_MASK;
    sign_state[y - 1][x    ] |= V_PVE_MASK;
    sign_state[y - 1][x + 1] |= NE_PVE_MASK;
    sign_state[y    ][x - 1] |= H_PVE_MASK;
    sign_state[y    ][x + 1] |= H_PVE_MASK;
    sign_state[y + 1][x - 1] |= NE_PVE_MASK;
    sign_state[y + 1][x    ] |= V_PVE_MASK;
    sign_state[y + 1][x + 1] |= NW_PVE_MASK;
  }
}

void update_child_state(int **child_state, int x, int y)
{
  child_state[y][x] |= PA_SIG;
}

void ezbc_subband_init(ezbc_subband_codec_t *codec,
		       ezbc_subband_codec_t *parent_codec,
		       ezbc_subband_codec_t *child_codec,
		       arithmetic_codec_t *arithmetic_codec,
		       ezbc_coeff_t *coeffs,
		       int width, int height,
		       int orientation,
		       int bpp,
		       int target_rate)
{
  int i, l;
  int depth;
  int total;

  codec->width = width;
  codec->height = height;
  codec->depth = depth = max(intlog2(width - 1), intlog2(height - 1)) + 2;
  codec->coeffs = coeffs;
  codec->bpp = bpp;
  codec->target_rate = target_rate;
  codec->arithmetic_codec = arithmetic_codec;
  codec->parent_codec = parent_codec; /* codec of the parent subband (if any)*/
  codec->child_codec = child_codec; /* codec of the child subband (if any)*/

  /* lists */
  codec->LIP = (ezbc_point_t *) malloc(width * height * sizeof(ezbc_point_t));
  codec->LIP_end = codec->LIP + width * height;
  codec->LSP = codec->LIP;
  codec->LSP_end = codec->LSP - 1;
  codec->LIS = (ezbc_point_t **) malloc(depth  * sizeof(ezbc_point_t *));
  codec->LIS_end = (ezbc_point_t **) malloc(depth * sizeof(ezbc_point_t *));

  for(l = 0; l < depth; l++) {
    int w, h;
  
    w = (( width + (1 << l) - 1) >> l); /* round up */
    h = ((height + (1 << l) - 1) >> l); /* round up */

    codec->LIS[l] = (ezbc_point_t *) malloc(w * h * sizeof(ezbc_point_t));
    codec->LIS_end[l] = codec->LIS[l] - 1;
  }
  codec->LSP_plane = 0;
  codec->LSP_mark = (ezbc_point_t*) malloc(bpp * sizeof(int));
  codec->LSP_split = codec->LSP;

  codec->LSP_bit_index_marks = (ezbc_point_t **) malloc(bpp * sizeof(ezbc_point_t *));
  memset(codec->LSP_bit_index_marks, 0, bpp * sizeof(ezbc_point_t *));

  codec->LSP_ids = (ezbc_point_t ***) malloc(bpp * sizeof(ezbc_point_t *));
  codec->LSP_ids[0] = (ezbc_point_t **) malloc(depth * bpp * sizeof(ezbc_point_t *));
  memset(codec->LSP_ids[0], 0, depth * bpp * sizeof(ezbc_point_t *));

  for(i = 1; i < bpp; i++)
    codec->LSP_ids[i] = codec->LSP_ids[i - 1] + depth;

  /* stack of insignificant sets */
  codec->LIS_stack = (ezbc_node_t *) malloc(4 * depth * sizeof(ezbc_node_t));
  codec->LIS_stack_top = codec->LIS_stack - 1;

  /* contexts */
  codec->context_sig = (context_t *) malloc(depth * sizeof(context_t));
  codec->context_jsig = (context_t *) malloc(4 * depth * sizeof(context_t));
  codec->context_node = (context_t *) malloc(depth * sizeof(context_t));

  /* create context models */
  switch(orientation) {
    case EZBC_ORIENT_DIAGONAL:
      ezbc_attach_cxts_diag(&codec->tables);
      break;
    case EZBC_ORIENT_HORIZONTAL:
    case EZBC_ORIENT_VERTICAL:
      ezbc_attach_cxts_main(&codec->tables);
      break;
    case EZBC_ORIENT_HORIZONTALVERTICAL:
      ezbc_attach_cxts_hv(&codec->tables);
      break;
    case EZBC_ORIENT_OTHER:
      ezbc_attach_cxts_other(&codec->tables);
      break;
    case EZBC_ORIENT_NONE:
      ezbc_attach_cxts_none(&codec->tables);
      break;
    default:
      fprintf(stderr, "unknown subband orientation %d\n", orientation);
      exit(1);
  }

  /* signifiance coding context models */
  total = 0;
  for(l = 0; l < depth; l++) {
    codec->context_sig[l].offset = total;
    for(i = 0; i < 4; i++) {
#if 0 // TEMP: reduce drastically the number of contexts
      switch(l) {
        case 0: /* level 0: uses jsig0 tables */
	  codec->context_jsig[4 * l + i].offset = total;
	  codec->context_jsig[4 * l + i].table = codec->tables.jsig0_zc_cxts[i];
	  codec->context_jsig[4 * l + i].count = codec->tables.jsig0_zc_cxts_count[i];
	  break;
        case 1: /* levels 1-2: uses jsig tables */
        case 2:
	  codec->context_jsig[4 * l + i].offset = total;
	  codec->context_jsig[4 * l + i].table = codec->tables.jsig_zc_cxts[i];
	  codec->context_jsig[4 * l + i].count = codec->tables.jsig_zc_cxts_count[i];
	  break;
        default: /* levels >2: uses same tables as level 2 */
	  codec->context_jsig[4 * l + i].offset = codec->context_jsig[4 * 2 + i].offset;
	  codec->context_jsig[4 * l + i].table = codec->context_jsig[4 * 2 + i].table;
	  codec->context_jsig[4 * l + i].count = 0;
	  break;
      }
#else
      /* keep sane on the number of contexts */
      if(l == 0) {
	  codec->context_jsig[i].offset = total;
	  codec->context_jsig[i].table = codec->tables.jsig0_zc_cxts[i];
	  codec->context_jsig[i].count = codec->tables.jsig0_zc_cxts_count[i];
      } else {
	  codec->context_jsig[4*l+i].offset = codec->context_jsig[i].offset;
	  codec->context_jsig[4*l+i].table = codec->context_jsig[i].table;
	  codec->context_jsig[4*l+i].count = 0;
      }
#endif
      total += codec->context_jsig[4 * l + i].count;
    }
    codec->context_sig[l].count = total - codec->context_sig[l].offset;
  }

  /* List of Insignificant Pixels context model */
  codec->context_LIP.offset = total;
  codec->context_LIP.table = codec->tables.LIP_zc_cxts;
  codec->context_LIP.count = codec->tables.LIP_zc_cxts_count;
  total += codec->context_LIP.count;

  /* qtree nodes context models */
  for(l = 1; l < depth; l++) {
    codec->context_node[l].offset = total;
    codec->context_node[l].table = codec->tables.node_zc_cxts;
    if(l < 3)
      codec->context_node[l].count = codec->tables.node_zc_cxts_count;
    else
      codec->context_node[l].count = 0;
    total += codec->context_node[l].count;
  }

  /* List of Significant Pixels */
  codec->context_LSP.offset = total;
  codec->context_LSP.count = 79;
  codec->context_LSP.table = codec->tables.LSP_zc_cxts;
  total += codec->context_LSP.count;

  /* Sign coding */
  codec->context_sign.offset = total;
  codec->context_sign.table = codec->tables.sc_cxts;
  codec->context_sign.count = codec->tables.sc_cxts_count;
  total += codec->context_sign.count;
  codec->context_count = total;

  codec->context_models = (binary_model_t *) malloc(total * sizeof(binary_model_t));
  for(i = 0; i < total; i++)
    binary_model_init(codec->context_models[i]);

  /* Allocate storage for current state */
  codec->base_state = (int **) malloc((height + 2 * EZBC_CONTEXT_SIZE) * sizeof(int *));
  codec->base_state[0] = (int *) malloc((width + 2 * EZBC_CONTEXT_SIZE) * 
					  (height + 2 * EZBC_CONTEXT_SIZE) * sizeof(int));
  codec->base_state[0] += EZBC_CONTEXT_SIZE; // align properly
  for(i = 1; i < (height + 2 * EZBC_CONTEXT_SIZE); i++)
    codec->base_state[i] = codec->base_state[i-1] + (width + 2 * EZBC_CONTEXT_SIZE);
  codec->base_state += EZBC_CONTEXT_SIZE; // align properly

  codec->sign_state = (int **) malloc((height + 2 * EZBC_CONTEXT_SIZE) * sizeof(int *));
  codec->sign_state[0] = (int *) malloc((width + 2 * EZBC_CONTEXT_SIZE) * 
					  (height + 2 * EZBC_CONTEXT_SIZE) * sizeof(int));
  codec->sign_state[0] += EZBC_CONTEXT_SIZE; // align properly
  for(i = 1; i < (height + 2 * EZBC_CONTEXT_SIZE); i++)
    codec->sign_state[i] = codec->sign_state[i-1] + (width + 2 * EZBC_CONTEXT_SIZE);
  codec->sign_state += EZBC_CONTEXT_SIZE; // align properly

  codec->node_state = (int ***) malloc((depth + 1) * sizeof(int **));
  for(l = 0; l < depth + 1; l++) {
    int w, h;
  
    w = (( width + (1 << l) - 1) >> l) + 2 * EZBC_CONTEXT_SIZE; /* round up + room for borders */
    h = ((height + (1 << l) - 1) >> l) + 2 * EZBC_CONTEXT_SIZE; /* round up + room for borders */

    codec->node_state[l] = (int **) malloc(h * sizeof(int *));
    codec->node_state[l][0] = (int *) malloc(w * h * sizeof(int));
    codec->node_state[l][0] += EZBC_CONTEXT_SIZE; // align properly
    for(i = 1; i < h; i++)
      codec->node_state[l][i] = codec->node_state[l][i-1] + w;
    codec->node_state[l] += EZBC_CONTEXT_SIZE; // align properly
  }

  /* allocate coefficient quad tree. */
  codec->nodes = (ezbc_coeff_t ***) malloc((depth + 1) * sizeof(ezbc_coeff_t **));
  for(l = 1; l < depth + 1; l++) {
    int w, h;
  
    w = (( width + (1 << l) - 1) >> l); /* round up */
    h = ((height + (1 << l) - 1) >> l); /* round up */

    codec->nodes[l] = (ezbc_coeff_t **) malloc(h * sizeof(ezbc_coeff_t *));
    codec->nodes[l][0] = (ezbc_coeff_t *) malloc(w * h * sizeof(ezbc_coeff_t));
    for(i = 1; i < h; i++)
      codec->nodes[l][i] = codec->nodes[l][i-1] + w;
  }

  codec->nodes[0] = (ezbc_coeff_t **) malloc(height * sizeof(ezbc_coeff_t *));
  codec->nodes[0][0] = (ezbc_coeff_t *) coeffs;
  for(i = 1; i < height; i++)
    codec->nodes[0][i] = codec->nodes[0][i-1] + width;

  /* initialize the codec state */
  init_state(codec);

  /* init the node quad tree */
  init_node_qtree(codec);

  codec->binbook = NULL;
  codec->hash = NULL;
}

void ezbc_subband_cleanup(ezbc_subband_codec_t *codec)
{
  int l;
  int depth = codec->depth;

  free(codec->LIP);

  for(l = 0; l < depth; l++)
    free(codec->LIS[l]);

  free(codec->LIS);
  free(codec->LIS_end);

  free(codec->LSP_mark);

  free(codec->LSP_bit_index_marks);
  free(codec->LSP_ids[0]);
  free(codec->LSP_ids);

  free(codec->LIS_stack);

  free(codec->context_sig);
  free(codec->context_jsig);
  free(codec->context_node);

  free(codec->context_models);

  codec->base_state -= EZBC_CONTEXT_SIZE; // align properly
  codec->base_state[0] -= EZBC_CONTEXT_SIZE; // align properly
  free(codec->base_state[0]);
  free(codec->base_state);

  codec->sign_state -= EZBC_CONTEXT_SIZE; // align properly
  codec->sign_state[0] -= EZBC_CONTEXT_SIZE; // align properly
  free(codec->sign_state[0]);
  free(codec->sign_state);

  for(l = 0; l < depth + 1; l++) {
    codec->node_state[l] -= EZBC_CONTEXT_SIZE; // align properly
    codec->node_state[l][0] -= EZBC_CONTEXT_SIZE; // align properly
    free(codec->node_state[l][0]);
    free(codec->node_state[l]);
  }
  free(codec->node_state);

  for(l = 1; l < depth + 1; l++) {
    free(codec->nodes[l][0]);
    free(codec->nodes[l]);
  }
  free(codec->nodes[0]);

  free(codec->nodes);
}

/* sorting functions */
typedef struct {
  ezbc_coeff_t coeff;
  ezbc_point_t point;
} ezbc_point_coeff_t;

static int sort_list_comp(const void *e1, const void *e2)
{
  ezbc_coeff_t c1 = ((const ezbc_point_coeff_t *) e1)->coeff;
  ezbc_coeff_t c2 = ((const ezbc_point_coeff_t *) e2)->coeff;

  /* sort in descending order */
  return(!(c1 < c2));
}

void ezbc_sort_list(ezbc_point_t *list, ezbc_coeff_t **coeffs, int length, int index)
{
  int i;
  ezbc_point_coeff_t *L;

  if(!length) return;

  /* TODO: far from being the fastest way of sorting */
  /* ideally lists should be pointer lists when x,y are not necessary */
  /* (as is the case for LSP at least) */
  L = (ezbc_point_coeff_t *) malloc(length * sizeof(ezbc_point_coeff_t));

  /* associate points to coefficients */
  for(i = 0; i < length; i++) {
    int x, y;

    L[i].point = list[i];
    x = ezbc_coord_x(list[i]);
    y = ezbc_coord_y(list[i]);
    L[i].coeff = (coeffs[y][x] & EZBC_ABS) >> (index + 1);
  }

  qsort(L, length, sizeof(ezbc_point_coeff_t), sort_list_comp);

  for(i = 0; i < length; i++)
    list[i] = L[i].point;

  free(L);
}

void ezbc_init()
{
  ezbc_initialize_zc_cxts();
  ezbc_initialize_sc_cxts();
}

void ezbc_subband_start(ezbc_subband_codec_t *codec)
{
  int i, l, count, depth;
  binary_model_t *models;
  binary_model_t *parent_models;

  /* retrieve the context models from the parents */
  if(codec->parent_codec) {

    /* sign models */
    models = &codec->context_models[codec->context_sign.offset];
    parent_models = &codec->parent_codec->context_models[codec->parent_codec->context_sign.offset];
    count = codec->context_sign.count;
    
    for(i = 0; i < count; i++) {
      binary_model_copy(models[i], parent_models[i]);
      binary_model_scale(models[i]);
    }

    depth = codec->parent_codec->depth;
    for(l = 0; l < depth; l++) {
      models = &codec->context_models[codec->context_sig[l].offset];
      parent_models = &codec->parent_codec->context_models[codec->parent_codec->context_sig[l].offset];
      count = codec->context_sig[l].count;

      for(i = 0; i < count; i++) {
	binary_model_copy(models[i], parent_models[i]);
	binary_model_scale(models[i]);
      }
    }
  }
}
