/*
   libit - Library for basic source and channel coding functions
   Copyright (C) 2005-2006 Vivien Chappelier, Herve Jegou

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
  Variable length codes
  Copyright (C) 2005-2006 Herve Jegou
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <it/vlc.h>
#include <it/io.h>
#include <it/math.h>

/*------------------------------------------------------------------------------*/
vlc_t *vlc_new (int n)
{
  vlc_t *vlc = (vlc_t *) malloc (sizeof (vlc_t));
  assert (vlc);

  vlc->nb_symb = n;
  vlc->nb_nodes = 2 * n - 1;	/* By default, we assume that the kraft equality is verified */
  vlc->node_root = 2 * n - 2;

  vlc->cwd_length = ivec_new_zeros (vlc->nb_nodes);
  vlc->cwd = ivec_new_zeros (vlc->nb_nodes);
  vlc->map[0] = ivec_new (vlc->nb_nodes);
  vlc->map[1] = ivec_new (vlc->nb_nodes);
  ivec_set (vlc->map[0], DUMM_NODE);
  ivec_set (vlc->map[1], DUMM_NODE);

  return vlc;
}


/*------------------------------------------------------------------------------*/
vlc_t *vlc_clone (const vlc_t * vlc)
{
  vlc_t *v = (vlc_t *) malloc (sizeof (vlc_t *));
  assert (v);

  v->nb_symb = vlc->nb_symb;
  v->nb_nodes = vlc->nb_nodes;
  v->node_root = vlc->node_root;

  v->cwd_length = ivec_clone (v->cwd_length);
  v->cwd = ivec_clone (v->cwd);
  v->map[0] = ivec_clone (v->map[0]);
  v->map[1] = ivec_clone (v->map[0]);

  return v;
}


/*------------------------------------------------------------------------------*/
void vlc_delete (vlc_t * vlc)
{
  assert (vlc);
  ivec_delete (vlc->cwd_length);
  ivec_delete (vlc->cwd);
  ivec_delete (vlc->map[0]);
  ivec_delete (vlc->map[1]);
  free (vlc);
}


/*------------------------------------------------------------------------------*/
vlc_t *vlc_flc (int nb_bits)
{
  vlc_t *flc = vlc_new (1 << nb_bits);
  int  l, l_start, l_start_next;
  int  nl, child0, child1;

  for (l = 0; l < nb_bits; l++) {

    l_start = (1 << l) - 1;
    l_start_next = (1 << (l + 1)) - 1;

    for (nl = l_start; nl < l_start_next; nl++) {
      child1 = (nl - l_start) * 2 + l_start_next;
      child0 = (nl - l_start) * 2 + 1 + l_start_next;

      flc->map[0][flc->node_root - nl] = flc->node_root - child0;
      flc->map[1][flc->node_root - nl] = flc->node_root - child1;
    }

  }

  vlc_affect_codewords (flc);
  return flc;
}

/*------------------------------------------------------------------------------*/
vlc_t *vlc_huffman (const vec pdf)
{
  int  i, n, argmin1, argmin2;
  double min1, min2;
  vlc_t *vlc = vlc_new (Vec_length (pdf));
  vec  probas = vec_new_zeros (vlc->nb_nodes);	/* The probability of the nodes are included in probas */
  /* Active_symbols allows to know which symbol/node must be added to the tree */
  bvec active_symbols = bvec_new_zeros (vlc->nb_nodes);
  bvec_set_between (active_symbols, 0, vlc->nb_symb - 1, 1);

  for (i = 0; i < vlc->nb_symb; i++)
    probas[i] = pdf[i];

  /* Construction of the tree */
  for (n = vlc->nb_symb; n < vlc->nb_nodes; n++) {
    argmin1 = -1;
    argmin2 = -1;
    min1 = 1;
    min2 = 2;
    for (i = 0; i < vlc->nb_nodes; i++) {
      if (active_symbols[i]) {
	if (probas[i] < min1) {
	  min2 = min1;
	  argmin2 = argmin1;
	  min1 = probas[i];
	  argmin1 = i;
	}
	else if (probas[i] < min2) {
	  min2 = probas[i];
	  argmin2 = i;
	}
      }
    }

    vlc->map[0][n] = argmin1;
    vlc->map[1][n] = argmin2;
    active_symbols[argmin1] = 0;
    active_symbols[argmin2] = 0;

    /* Process the proba of the new symbol */
    probas[n] = probas[argmin1] + probas[argmin2];

    /* Let the new symbol be active and the symbols min1 and min2 be inactive */
    active_symbols[n] = 1;
  }
  vec_delete (probas);
  bvec_delete (active_symbols);
  vlc_affect_codewords (vlc);
  return vlc;
}


/*--------------------------------------------------------------*/
/* Create the VLC code according to the Hu-Tucker algorithm
   The algorithm assumes that the probability distribution function is pdf.
   Note: this implementation could (and should) be optimized in terms of computionnal cost */
vlc_t *vlc_hu_tucker (const vec pdf)
{
  int  i, j, k, i1, i2, n, node1, node2;
  vlc_t *vlc = vlc_new (vec_length (pdf));
  ivec active_symbols = ivec_new_arithm (0, 1, vlc->nb_symb);
  vec  probas = vec_new_zeros (vlc->nb_nodes);	/* The probability of the nodes are included in probas */
  int  nrof_act = vlc->nb_symb;	/* number of active symbols */
  vec  dumm_vec_symb = vec_new_1N (vlc->nb_symb);

  for (i = 0; i < vlc->nb_symb; i++)
    probas[i] = pdf[i];

  for (n = vlc->nb_symb; n < vlc->nb_nodes; n++) {
    /* At this point, all couples are assumed compatible    */
    ivec merg_comp = ivec_new_ones (nrof_act * nrof_act);

    /* Store the probability that is the smallest among compatible nodes and the corresponding indexes */
    ivec best_for_i = ivec_new (nrof_act);
    vec  best_for_i_prob = vec_new (nrof_act);
    ivec_set (best_for_i, DUMM_NODE);
    vec_set (best_for_i_prob, 1.1);

    /* Process the merging compatibilities */
    for (i = 0; i < nrof_act; i++) {
      for (j = i + 1; j < nrof_act; j++) {
	/* Test if there is an original block between i and j */
	for (k = i + 1; k < j; k++)
	  if (vlc_is_leaf (vlc, active_symbols[k])) {
	    merg_comp[i * nrof_act + j] = 0;
	    merg_comp[j * nrof_act + i] = 0;
	  }
      }
      /* A node can not be merged with itself */
      merg_comp[i * nrof_act + i] = 0;
    }

    /* Retrieve the couple corresponding to smallest probabilities */
    for (i = 0; i < nrof_act; i++)
      for (j = 0; j < nrof_act; j++)

	if (merg_comp[i * nrof_act + j] == 1)
	  if (probas[active_symbols[j]] < best_for_i_prob[i]) {
	    best_for_i[i] = j;
	    best_for_i_prob[i] = probas[active_symbols[j]];
	  }

    /* Retrieve at least one of the candidate couples */
    /* The new index is assigned the index n          */
    i1 = -1, i2 = -1;
    for (i = 0; i < nrof_act && i1 < 0; i++)
      for (j = 0; j < nrof_act && i1 < 0; j++)
	if (merg_comp[i * nrof_act + j]
	    && (best_for_i[i] == j) && (best_for_i[j] == i)) {
	  i1 = i;
	  i2 = j;
	}

    node1 = active_symbols[i1];
    node2 = active_symbols[i2];
    vlc->map[0][n] = node1;
    vlc->map[1][n] = node2;

    /* Process the proba of the new symbol */
    probas[n] = probas[node1] + probas[node2];

    /* Update the list of active symbols */
    active_symbols[i1] = n;
    ivec_del (active_symbols, i2);
    ivec_delete (merg_comp);
    ivec_delete (best_for_i);
    vec_delete (best_for_i_prob);

    nrof_act--;
  }

  vlc_affect_codewords (vlc);
  vlc_quasi_lexicographic (vlc, pdf, dumm_vec_symb);

  vec_delete (dumm_vec_symb);
  vec_delete (probas);
  ivec_delete (active_symbols);
  return vlc;
}


/*------------------------------------------------------------------------------*/
static void vlc_affect_codewords_subtree (vlc_t * vlc, int s)
{
  /* Symbol is a dumm symbol */
  if (s == DUMM_NODE)
    return;

  /* Symbol is a leaf */
  if (!vlc_is_leaf (vlc, s)) {
    /* Otherwise, internal node */
    int  child0 = vlc_get_child0 (vlc, s);
    int  child1 = vlc_get_child1 (vlc, s);

    /* The child0 is a valid symbols/leaf */
    if (child0 != DUMM_NODE) {
      vlc->cwd_length[child0] = vlc->cwd_length[s] + 1;
      vlc->cwd[child0] = vlc->cwd[s];
      vlc_affect_codewords_subtree (vlc, child0);
    }

    /* Same story */
    if (child1 != DUMM_NODE) {
      vlc->cwd_length[child1] = vlc->cwd_length[s] + 1;
      vlc->cwd[child1] = vlc->cwd[s] + (1 << (vlc->cwd_length[s]));
      vlc_affect_codewords_subtree (vlc, child1);
    }
  }
}


/*------------------------------------------------------------------------------*/
/* To call this function, the map table must be correctly set, and the 
   pointer of variables codes and cwd_length must be already allocated 
   (It is the case if Function vlc_new has been called)                         */
void vlc_affect_codewords (vlc_t * vlc)
{
  vlc->cwd[vlc->node_root] = 0;
  vlc->cwd_length[vlc->node_root] = 0;
  vlc_affect_codewords_subtree (vlc, vlc->node_root);
}


/*------------------------------------------------------------------------------*/
int vlc_minh (const vlc_t * vlc)
{
  int  the_min = 1000;
  int  i;
  for (i = 0; i < vlc->nb_symb; i++)
    if (the_min > vlc->cwd_length[i])
      the_min = vlc->cwd_length[i];
  return the_min;
}


/*------------------------------------------------------------------------------*/
int vlc_maxh (const vlc_t * vlc)
{
  int  the_max = 0;
  int  i;
  for (i = 0; i < vlc->nb_symb; i++)
    if (the_max < vlc->cwd_length[i])
      the_max = vlc->cwd_length[i];
  return the_max;
}


/*------------------------------------------------------------------------------*/
double vlc_mdl (const vlc_t * vlc, const vec pdf)
{
  double mdl = 0;
  int  i;
  for (i = 0; i < vlc->nb_symb; i++)
    mdl += vlc->cwd_length[i] * pdf[i];
  return mdl;
}


/*------------------------------------------------------------------------------*/
double vlc_kraft_sum (const vlc_t * vlc)
{
  double the_sum = 0;
  int  i;
  for (i = 0; i < vlc->nb_symb; i++) {
    the_sum += pow (2, -vlc->cwd_length[i]);
  }

  return the_sum;
}


/*------------------------------------------------------------------------------*/
/* Tranform the tree in order to generate a quasi-lexicographic order           */
void vlc_quasi_lexicographic (const vlc_t * vlc, const vec pdf,
			      const vec symb)
{
  int  i, l, curnode = vlc->nb_symb, child0, child1;
  ivec tmp_idx, tmp_idx_sub;	/* Node/Leaf index  */
  vec  tmp_rec, tmp_rec_sub;	/* Reconstruction value for nodes/symbols */
  int  nbtmpsymb;		/* Number of symbols/node in this layer   */
  int  maxh;

  vec  probas = vec_new_zeros (vlc->nb_nodes);	/* The probability of the nodes are included in probas */
  vec  symbols = vec_new_zeros (vlc->nb_nodes);

  assert (vlc);
  assert (vlc_kraft_sum (vlc) == 1);

  for (i = 0; i < vlc->nb_symb; i++) {
    probas[i] = pdf[i];
    symbols[i] = symb[i];
  }
  ivec_set_between (vlc->cwd_length, vlc->nb_symb, vlc->nb_nodes - 1, 0);

  maxh = vlc_maxh (vlc);

  /* The highest size is allocated by default                                   */
  tmp_idx = ivec_new (vlc->nb_nodes);
  tmp_rec = vec_new (vlc->nb_nodes);

  for (l = maxh; l >= 0; l--) {
    ivec sort_order;

    /* Retrieve the nodes of the current layer */
    nbtmpsymb = 0;

    for (i = 0; i < vlc->nb_nodes; i++)
      if (vlc->cwd_length[i] == l) {
	/* Is it a node/leaf of the focused layer? */
	tmp_idx[nbtmpsymb] = i;
	tmp_rec[nbtmpsymb] = symbols[i];
	nbtmpsymb++;
      }

    /* We have to sort the vectors tmp_idx and tmp_rec according to
       the order generated by values of tmp_rec
     */
    Vec_length (tmp_rec) = nbtmpsymb;
    sort_order = vec_sort_index (tmp_rec);
    Vec_length (tmp_rec) = vlc->nb_nodes;

    tmp_idx_sub = ivec_index_by (tmp_idx, sort_order);
    tmp_rec_sub = vec_index_by (tmp_rec, sort_order);

    memcpy (tmp_idx, tmp_idx_sub, nbtmpsymb * sizeof (*tmp_idx));
    memcpy (tmp_rec, tmp_rec_sub, nbtmpsymb * sizeof (*tmp_idx));

    /* Aggregation of symbols of the current layer by pairs */
    for (i = 0; i < nbtmpsymb / 2; i++) {
      child0 = tmp_idx[i * 2];
      child1 = tmp_idx[i * 2 + 1];

      /* Update the map table */
      vlc->map[0][curnode] = child0;
      vlc->map[1][curnode] = child1;

      /* probas, expectation, length,  of this node */
      probas[curnode] = probas[child0] + probas[child1];
      symbols[curnode] = (probas[child0] * symbols[child0]
			  +
			  probas[child1] * symbols[child1]) /
	(probas[curnode]);
      vlc->cwd_length[curnode] = vlc->cwd_length[child0] - 1;

      /* update the structure used for the next sort */
      tmp_idx[curnode] = curnode;
      tmp_rec[curnode] = symbols[curnode];

      curnode++;
    }

    ivec_delete (tmp_idx_sub);
    vec_delete (tmp_rec_sub);
    ivec_delete (sort_order);
  }

  vlc->cwd[vlc->node_root] = 0;

  for (i = vlc->nb_nodes - 1; i >= vlc->nb_symb; i--) {
    child0 = vlc->map[0][i];
    child1 = vlc->map[1][i];
    vlc->cwd[child0] = vlc->cwd[i];
    vlc->cwd[child1] = vlc->cwd[i] + (1 << (vlc->cwd_length[i]));
  }

  ivec_delete (tmp_idx);
  vec_delete (tmp_rec);
  vec_delete (probas);
  vec_delete (symbols);
}


/*------------------------------------------------------------------------------*/
/* Return the respective probabilities of the nodes                              */

static void vlc_nodes_pdf_tmp (const vlc_t * vlc, int node, vec pdf)
{
  int  child0, child1;

  if (node < vlc->nb_symb)
    return;

  child0 = vlc_get_child0 (vlc, node);
  child1 = vlc_get_child1 (vlc, node);

  vlc_nodes_pdf_tmp (vlc, child0, pdf);
  vlc_nodes_pdf_tmp (vlc, child1, pdf);

  if (child0 == DUMM_NODE) {
    pdf[node] = pdf[child1];
  }
  else if (child1 == DUMM_NODE) {
    pdf[node] = pdf[child0];
  }
  else {
    pdf[node] = pdf[child0] + pdf[child1];
  }
}


vec vlc_nodes_pdf (const vlc_t * vlc, vec pdf)
{
  int  i;
  vec  probas = vec_new_zeros (vlc->nb_nodes);	/* The probability of the nodes are included in probas */

  for (i = 0; i < vlc->nb_symb; i++)
    probas[i] = pdf[i];

  vlc_nodes_pdf_tmp (vlc, vlc->node_root, probas);
  return probas;
}


vec vlc_nodes_proba0 (const vlc_t * vlc, vec pdf)
{
  int  i, child0, child1;
  vec  nodes_pdf = vlc_nodes_pdf (vlc, pdf);	/* The probability of the nodes are included in probas */
  vec  proba0 = vec_new (vlc->nb_nodes - vlc->nb_symb);
  double pchild0, pchild1;

  for (i = vlc->nb_symb; i < vlc->nb_nodes; i++) {
    child0 = vlc->map[0][i];
    child1 = vlc->map[1][i];

    if (child0 != DUMM_NODE)
      pchild0 = nodes_pdf[child0];
    else
      pchild0 = 0;

    if (child1 != DUMM_NODE)
      pchild1 = nodes_pdf[child1];
    else
      pchild1 = 0;

    proba0[i - vlc->nb_symb] = pchild0 / (pchild0 + pchild1);
  }

  vec_delete (nodes_pdf);
  return proba0;
}


/*------------------------------------------------------------------------------*/
/* Return the respective expectations of the nodes                              */

static void vlc_nodes_expectation_tmp (const vlc_t * vlc,
				       int node, vec pdf, vec symbols)
{
  int  child0, child1;

  if (node < vlc->nb_symb)
    return;

  child0 = vlc_get_child0 (vlc, node);
  child1 = vlc_get_child1 (vlc, node);

  vlc_nodes_expectation_tmp (vlc, child0, pdf, symbols);
  vlc_nodes_expectation_tmp (vlc, child1, pdf, symbols);

  if (child0 == DUMM_NODE) {
    pdf[node] = pdf[vlc_get_child1 (vlc, node)];
    symbols[node] = symbols[vlc_get_child1 (vlc, node)];
  }
  else if (child1 == DUMM_NODE) {
    pdf[node] = pdf[child0];
    symbols[node] = symbols[child0];
  }
  else {
    pdf[node] = pdf[child0] + pdf[child1];
    symbols[node] = (symbols[child0] * pdf[child0]
		     + symbols[child1] * pdf[child1]) / pdf[node];
  }
}


vec vlc_nodes_expectation (const vlc_t * vlc, vec pdf, vec symb)
{
  int  i;
  vec  probas = vec_new_zeros (vlc->nb_nodes);	/* The probability of the nodes are included in probas */
  vec  symbols = vec_new_zeros (vlc->nb_nodes);

  for (i = 0; i < vlc->nb_symb; i++) {
    probas[i] = pdf[i];
    symbols[i] = symb[i];
  }

  vlc_nodes_expectation_tmp (vlc, vlc->node_root, probas, symbols);
  vec_delete (probas);
  return symbols;
}


/*------------------------------------------------------------------*/
vec vlc_nodes_entropy (const vlc_t * vlc, vec pdf)
{
  int  i, child0, child1;
  vec  probas = vlc_nodes_pdf (vlc, pdf);
  vec  H = vec_new_zeros (vec_length (probas));
  double p0, p1;

  for (i = vlc->nb_symb; i < vlc->nb_nodes; i++) {
    child0 = vlc_get_child0 (vlc, i);
    child1 = vlc_get_child1 (vlc, i);

    if (child0 == DUMM_NODE || child1 == DUMM_NODE || probas[child0] == 0
	|| probas[child1] == 0 || probas[i] == 0)
      H[i] = 0;
    else {
      p0 = probas[child0] / probas[i];
      p1 = probas[child1] / probas[i];
      H[i] = -p0 * log2 (p0) - p1 * log2 (p1);
    }
  }

  vec_delete (probas);
  return H;
}


/*------------------------------------------------------------------*/
/* Return the respective variance of the nodes                      */

/* process the contribution of the leaves of the node node to the expectation of 
   the mean error corresponding to the reconstruction value rec_value )         */
static double vlc_node_variance_tmp (const vlc_t * vlc,
				     int node,
				     double rec_value,
				     vec pdf, vec symb, double *p_sum_pdf)
{
  double r = 0;
  int  child0, child1;

  if (node < vlc->nb_symb) {
    *p_sum_pdf += pdf[node];
    return pdf[node] * (symb[node] - rec_value) * (symb[node] - rec_value);
  }

  /* At this point, we know that we do not consider a leaf */
  child0 = vlc_get_child0 (vlc, node);
  child1 = vlc_get_child1 (vlc, node);

  if (child0 != DUMM_NODE)
    r += vlc_node_variance_tmp (vlc, child0, rec_value, pdf, symb, p_sum_pdf);

  if (child1 != DUMM_NODE)
    r += vlc_node_variance_tmp (vlc, child1, rec_value, pdf, symb, p_sum_pdf);

  return r;
}


double vlc_node_variance (const vlc_t * vlc, int node, vec pdf, vec symb)
{
  double sum_pdf = 0;
  double v =
    vlc_node_variance_tmp (vlc, node, symb[node], pdf, symb, &sum_pdf);

  if (sum_pdf > 0)
    return v / sum_pdf;
  else
    return 0;
}


vec vlc_nodes_variance (const vlc_t * vlc, vec pdf, vec symb)
{
  vec  E = vlc_nodes_expectation (vlc, pdf, symb);
  int  i;
  vec  variances = vec_new_zeros (vlc->nb_nodes);
  for (i = vlc->nb_symb; i < vlc->nb_nodes; i++)
    variances[i] = vlc_node_variance (vlc, i, pdf, E);

  vec_delete (E);
  return variances;
}


/*------------------------------------------------------------------*/
/* Expectation of the decrease of MSE corresponding to internal nodes */
vec vlc_nodes_delta_energy (const vlc_t * vlc, vec pdf, vec symb)
{
  int  i, child0, child1;
  double v0, v1, p0, p1;
  vec  P = vlc_nodes_pdf (vlc, pdf);
  vec  V = vlc_nodes_variance (vlc, pdf, symb);
  vec  D = vec_new_zeros (vlc->nb_nodes);

  for (i = vlc->nb_symb; i < vlc->nb_nodes; i++) {
    child0 = vlc_get_child0 (vlc, i);
    child1 = vlc_get_child1 (vlc, i);

    if (child0 == DUMM_NODE) {
      v0 = 0;
      p0 = 0;
    }
    else {
      v0 = V[child0];
      p0 = P[child0];
    }

    if (child1 == DUMM_NODE) {
      v1 = 0;
      p1 = 0;
    }
    else {
      v1 = V[child1];
      p1 = P[child1];
    }
    D[i] = V[i] - (p0 * v0 + p1 * v1) / P[i];
  }

  vec_delete (P);
  vec_delete (V);
  return D;
}


/*------------------------------------------------------------------*/
void vlc_print (const vlc_t * vlc)
{
  int  i, c;
  assert (vlc);

  printf ("{");
  for (i = 0; i < vlc->nb_symb; i++) {
    int  code;

    if (i > 0)
      printf (" ");

    code = vlc_get_cwd (vlc, i);

    if (vlc->cwd_length[i] == 0)	/* Should not occur every day... */
      printf ("/");
    else
      for (c = 0; c < vlc->cwd_length[i]; c++) {
	char bit = (((code % 2) > 0) ? 1 : 0);
	code >>= 1;
	printf ("%d", bit);
      }
  }
  printf ("}");
}


/*------------------------------------------------------------------*/
vlc_t *vlc_read (const char *svlc)
{
  int  node, symb, bit, child_node;
  char *s_start = strpbrk (svlc, "{");
  char *s_end = strpbrk (svlc, "}");
  char *s;
  vlc_t *vlc = (vlc_t *) malloc (sizeof (vlc_t));

  assert (vlc);
  if (s_start == NULL || s_end == NULL || s_start >= s_end)
    it_error ("String " "%s" " does not represent a valid vlc\n", s_start);

  s = strpbrk (s_start, "01");
  if (s == NULL)
    it_error ("String " "%s" " does not represent a valid vlc\n", s_start);

  /* First count the number of symbols */
  vlc->nb_symb = 0;

  while (s < s_end && s != NULL) {
    s = strpbrk (s, " ,}");
    s = strpbrk (s, "01");
    vlc->nb_symb++;
  }

  /* By default, only the root node is allocated */
  vlc->nb_nodes = vlc->nb_symb + 1;
  vlc->node_root = vlc->nb_symb;

  vlc->cwd_length = ivec_new_alloc (vlc->nb_symb + 1, vlc->nb_symb + 1);
  vlc->cwd = ivec_new_alloc (vlc->nb_symb + 1, vlc->nb_symb + 1);
  vlc->map[0] = ivec_new_alloc (vlc->nb_symb + 1, vlc->nb_symb + 1);
  vlc->map[1] = ivec_new_alloc (vlc->nb_symb + 1, vlc->nb_symb + 1);

  vlc->cwd_length[vlc->node_root] = 0;
  vlc->cwd[vlc->node_root] = 0;
  vlc->map[0][vlc->node_root] = DUMM_NODE;
  vlc->map[1][vlc->node_root] = DUMM_NODE;

  symb = 0;
  s = strpbrk (s_start, "01");
  while (s < s_end && s != NULL) {

    node = vlc->node_root;

    while (*s == '0' || *s == '1') {
      bit = (int) (*s - '0');

      /* if node does not exists => allocated a new node */
      if (vlc->map[bit][node] == DUMM_NODE) {
	if (s[1] == '0' || s[1] == '1') {	/* Internal node is added */
	  child_node = ivec_length (vlc->cwd);
	  ivec_push (vlc->cwd,
		     vlc->cwd[node] + (bit << vlc->cwd_length[node]));
	  ivec_push (vlc->cwd_length, vlc->cwd_length[node] + 1);
	  ivec_push (vlc->map[0], DUMM_NODE);
	  ivec_push (vlc->map[1], DUMM_NODE);
	  vlc->map[bit][node] = child_node;
	  vlc->nb_nodes++;
	}
	else {			/* Symbol */
	  child_node = symb++;
	  vlc->cwd[child_node] = vlc->cwd[node]
	    + (bit << vlc->cwd_length[node]);
	  vlc->cwd_length[child_node] = vlc->cwd_length[node] + 1;
	  vlc->map[bit][node] = child_node;
	}
      }

      node = vlc->map[bit][node];	/* Go tho child node */
      s++;
    }
    s = strpbrk (s, "01");
  }

  vlc_affect_codewords (vlc);

  return vlc;
}


/*------------------------------------------------------------------*/
ivec vlc_energy_order (const vlc_t * vlc, vec pdf, vec symb)
{
  ivec NO = ivec_new_alloc (0, vlc->nb_nodes - vlc->nb_symb);	/* Node order */
  ivec NI = ivec_new_alloc (0, vlc->nb_nodes - vlc->nb_symb);	/* Temporary */
  vec  RDNI, HNI;

  vec  D = vlc_nodes_delta_energy (vlc, pdf, symb);
  vec  H = vlc_nodes_entropy (vlc, pdf);

  int  n, node, child0, child1;	/* The node with the highest energy and its child */

  ivec_push (NI, vlc->node_root);

  while (ivec_length (NI) > 0) {
    RDNI = vec_index_by (D, NI);
    HNI = vec_index_by (H, NI);
    vec_div (RDNI, HNI);

    n = vec_max_index (RDNI);
    node = NI[n];

    child0 = vlc_get_child0 (vlc, node);
    child1 = vlc_get_child1 (vlc, node);

    ivec_push (NO, node);
    ivec_del (NI, n);

    if (child0 >= vlc->nb_symb)
      ivec_push (NI, child0);

    if (child1 >= vlc->nb_symb)
      ivec_push (NI, child1);

    vec_delete (RDNI);
    vec_delete (HNI);
  }

  ivec_delete (NI);
  vec_delete (D);
  vec_delete (H);

  return NO;
}


/*------------------------------------------------------------------*/
void vlc_print_all (const vlc_t * vlc, vec pdf, vec symb)
{
  int  i, c;
  vec  P = vlc_nodes_pdf (vlc, pdf);
  vec  E = vlc_nodes_expectation (vlc, pdf, symb);
  vec  V = vlc_nodes_variance (vlc, pdf, symb);
  vec  H = vlc_nodes_entropy (vlc, pdf);
  vec  D = vlc_nodes_delta_energy (vlc, pdf, symb);

  for (i = 0; i < vlc->nb_nodes; i++) {
    int  code;
    if (vlc_is_leaf (vlc, i))
      printf ("Symbol %d\t: %.3f\t%.7f \t\t\t\t\t\t", i, E[i], P[i]);
    else
      printf ("Node %d\t: %.3f\t%.7f\t%.5f\t%.3f\t%.4f\t%.3f\t",
	      i, E[i], P[i], H[i], V[i], D[i], D[i] / H[i]);

    code = vlc_get_cwd (vlc, i);

    if (vlc->cwd_length[i] == 0)
      printf ("/");
    else
      for (c = 0; c < vlc->cwd_length[i]; c++) {
	char bit = (((code % 2) > 0) ? 1 : 0);
	code >>= 1;
	printf ("%d", bit);
      }
    printf ("\n");
  }

  vec_delete (E);
  vec_delete (V);
  vec_delete (P);
  vec_delete (H);
  vec_delete (D);
}


/*------------------------------------------------------------------------------*/
int vlc_nb_bits_required (const vlc_t * vlc, ivec S)
{
  idx_t i;
  int  n = 0;
  for (i = 0; i < ivec_length (S); i++)
    n += vlc_get_cwd_length (vlc, S[i]);
  return n;
}


/*------------------------------------------------------------------------------*/
bvec vlc_encode_concat (const vlc_t * vlc, ivec S)
{
  idx_t t;			/* symbol clock */
  int  c, cwd;			/* codeword clock, codeword  */
  bvec E = bvec_new_alloc (0, ivec_length (S) * vlc_maxh (vlc));

  for (t = 0; t < ivec_length (S); t++) {
    cwd = vlc_get_cwd (vlc, S[t]);

    for (c = 0; c < vlc_get_cwd_length (vlc, S[t]); c++) {
      bvec_push (E, (char) (cwd & 1));
      cwd >>= 1;
    }
  }
  return E;
}


/*------------------------------------------------------------------------------*/
ivec vlc_decode_concat (const vlc_t * vlc, bvec E)
{
  idx_t b;			/* bit clock */
  ivec D = ivec_new_alloc (0, bvec_length (E) / vlc_maxh (vlc));
  int  n = vlc->node_root;	/* node */

  for (b = 0; b < bvec_length (E); b++) {
    n = vlc->map[E[b]][n];

    if (vlc_is_leaf (vlc, n)) {	/* The symbol is decoded */
      ivec_push (D, n);
      n = vlc->node_root;
    }
  }
  return D;
}


ivec vlc_decode_concat_N (const vlc_t * vlc, bvec E, idx_t N)
{
  idx_t b, t = 0;		/* bit and symbol clocks */
  ivec D = ivec_new (N);
  int  n = vlc->node_root;	/* node */
  ivec_set (D, vlc->node_root);

  for (b = 0; b < bvec_length (E) && t < N; b++) {
    n = vlc->map[E[b]][n];

    if (vlc_is_leaf (vlc, n)) {	/* The symbol is decoded */
      D[t++] = n;
      n = vlc->node_root;
    }
  }
  return D;
}

