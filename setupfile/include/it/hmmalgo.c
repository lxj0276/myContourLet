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
  Hidden Markov Model : estimation algorithms for Markov chains
  Copyright (C) 2005 Herve Jegou, Vivien Chappelier
*/

#include <it/hmmalgo.h>

/*---------------------------------------------------------------------------*/
/* This algorithm finds the a posteriori marginal distribution for each symbol
   of a sequence given the observation 'x' on each branch for each instant,
   the state model 'next_state' and transition probabilities 'next_state_pt'.
   The initial and ending state probabilities are given by alpha_0
   and beta_end respectively.
*/
mat bcjr (mat x, imat next_state, mat next_state_pt, vec alpha_0,
	  vec beta_end)
{
  int  m, mm, br;
  /* Number of states in the automaton */
  int  nb_states = imat_width (next_state);
  /* Number of branches from a node */
  int  nb_branch = imat_height (next_state);
  int  K = mat_width (x), t;
  double gamma;

  /* BCJR variables start at index 1 */
  mat  alpha = mat_new_zeros (K + 1, nb_states);
  mat  beta = mat_new_zeros (K + 1, nb_states);
  mat  pbranch = mat_new_zeros (nb_branch, K);
  vec  pbranch_sum;

  /* BCJR initialization */
  vec_copy (alpha[0], alpha_0);
  vec_copy (beta[K], beta_end);

  /* Forward pass */
  for (t = 1; t <= K; t++) {

    for (br = 0; br < nb_branch; br++) {
      for (mm = 0; mm < nb_states; mm++) {

	/* Check if the transition exists in the automaton */
	m = next_state[br][mm];
	if (m == DUMM_NODE)
	  continue;

	gamma = x[br][t - 1] * next_state_pt[br][mm];	/* Symbol clock differs by 1 */
	alpha[t][m] += alpha[t - 1][mm] * gamma;
      }
    }
    vec_normalize (alpha[t], 1);
  }

  /* Backward pass */
  for (t = K - 1; t > 0; t--) {
    for (br = 0; br < nb_branch; br++)
      for (m = 0; m < nb_states; m++) {
	mm = next_state[br][m];
	if (mm == DUMM_NODE)
	  continue;

	/* Here, this variable is used for symbol clock t+1 */
	gamma = x[br][t] * next_state_pt[br][m];
	beta[t][m] += beta[t + 1][mm] * gamma;
      }
    vec_normalize (beta[t], 1);
  }

  /* Process the metrics and project them on branchs */
  for (t = 1; t <= K; t++) {
    for (br = 0; br < nb_branch; br++)
      for (mm = 0; mm < nb_states; mm++) {
	m = next_state[br][mm];
	if (m == DUMM_NODE)
	  continue;

	gamma = x[br][t - 1] * next_state_pt[br][mm];
	pbranch[br][t - 1] += alpha[t - 1][mm] * beta[t][m] * gamma;
      }
  }

  pbranch_sum = mat_cols_sum (pbranch);

  for (t = 0; t < K; t++)
    mat_col_div_by (pbranch, t, pbranch_sum[t]);

  mat_delete (alpha);
  mat_delete (beta);
  vec_delete (pbranch_sum);

  return pbranch;
}


/*---------------------------------------------------------------------------*/
/* Maximum A Posteriori Sequence Decoding                                    */
/* This algorithm find the most probable sequence given the observation 'x'
   on each branch for each instant, the state model 'next_state' and
   transition probabilities 'next_state_pt'. The initial and ending state
   probabilities are given by alpha_0 and beta_end respectively.
*/
ivec viterbi (mat x, imat next_state, mat next_state_pt, vec alpha_0,
	      vec beta_end)
{
  int  m, mm, br;
  /* Number of states in the automaton */
  int  nb_states = imat_width (next_state);
  /* Number of branches from a node */
  int  nb_branch = imat_height (next_state);
  int  K = mat_width (x), t, cur_state;
  double rt, logalphacur, logp;

  /* BCJR variables start at index 1 */
  vec  logalpha_curr = vec_new (nb_states);	/* current value of logalpha    */
  vec  logalpha_next = vec_new (nb_states);	/* next value of logalpha       */
  vec  tmp;			/* use double-buffering for logalpha    */
  double gamma;			/* gamma needs only be computed locally */
  imat prev_branch = imat_new_set (DUMM_NODE, K + 1, nb_states);
  imat prev_state = imat_new_set (DUMM_NODE, K + 1, nb_states);
  ivec seq = ivec_new (K);

  /* Viterbi initialization */
  vec_copy (logalpha_curr, alpha_0);
  vec_apply_function (logalpha_curr, IT_FUNCTION (log), NULL);

  /* Forward pass */
  for (t = 1; t <= K; t++) {

    vec_set (logalpha_next, -HUGE_VAL);

    /* Compute the value gamma. Assume the loop on m' (represented by mm) */
    for (mm = 0; mm < nb_states; mm++)
      for (br = 0; br < nb_branch; br++) {

	/* Check if the transition exists in the automaton */
	m = next_state[br][mm];
	if (m == DUMM_NODE)
	  continue;

	rt = x[br][t - 1];	/* Symbol clock differs by 1 */

	gamma = rt * next_state_pt[br][mm];

	/* Check if the value of this path is better than
	   all the ones previously tested. */
	logalphacur = logalpha_curr[mm] + log (gamma);

	if (logalphacur > logalpha_next[m]) {
	  logalpha_next[m] = logalphacur;
	  prev_branch[t][m] = br;
	  prev_state[t][m] = mm;
	}
      }
    /* swap the buffers */
    tmp = logalpha_curr;
    logalpha_curr = logalpha_next;
    logalpha_next = tmp;
  }

  /* Find the ending state.
     For this purpose, the values of logalpha are modified
     with the a priori probability of the finishing states (beta_end).
     Process the corresponding probabilities of the paths              */
  logp = -HUGE_VAL;
  cur_state = -1;

  for (m = 0; m < nb_states; m++)
    if (logalpha_curr[m] + log (beta_end[m]) > logp) {
      logp = logalpha_curr[m] + log (beta_end[m]);
      cur_state = m;
    }

  /* Backward pass for reconstructing the sequence of branches */
  for (t = K - 1; t >= 0; t--) {
    seq[t] = prev_branch[t + 1][cur_state];
    cur_state = prev_state[t + 1][cur_state];
  }

  vec_delete (logalpha_curr);
  vec_delete (logalpha_next);
  imat_delete (prev_state);
  imat_delete (prev_branch);

  return seq;
}


/*----------------------------------------------------------------*/
/* same as viterbi() except all probabilities are expressed in log-domain */
/* This algorithm is the one that is usually refered to as
   "The Viterbi Algorithm" which finds the path of minimal cummulated metric.*/
ivec logviterbi (mat logx, imat next_state, mat next_state_logpt,
		 vec logalpha_0, vec logbeta_end)
{
  int  m, mm, br;
  /* Number of states in the automaton */
  int  nb_states = imat_width (next_state);
  /* Number of branches from a node */
  int  nb_branch = imat_height (next_state);
  int  K = mat_width (logx), t, cur_state;
  double rt, logalphacur, logp;

  /* BCJR variables start at index 1 */
  vec  logalpha_curr = vec_new (nb_states);	/* current value of logalpha    */
  vec  logalpha_next = vec_new (nb_states);	/* next value of logalpha       */
  vec  tmp;			/* use double-buffering for logalpha    */
  double loggamma;		/* gamma needs only be computed locally */
  imat prev_branch = imat_new_set (DUMM_NODE, K + 1, nb_states);
  imat prev_state = imat_new_set (DUMM_NODE, K + 1, nb_states);
  ivec seq = ivec_new (K);

  /* Viterbi initialization */
  vec_copy (logalpha_curr, logalpha_0);

  /* Forward pass */
  for (t = 1; t <= K; t++) {

    vec_set (logalpha_next, -HUGE_VAL);

    /* Compute the value gamma. Assume the loop on m' (represented by mm) */
    for (mm = 0; mm < nb_states; mm++)
      for (br = 0; br < nb_branch; br++) {

	/* Check if the transition exists in the automaton */
	m = next_state[br][mm];
	if (m == DUMM_NODE)
	  continue;

	rt = logx[br][t - 1];	/* Symbol clock differs by 1 */

	loggamma = rt + next_state_logpt[br][mm];

	/* Check if the value of this path is better than
	   all the ones previously tested. */
	logalphacur = logalpha_curr[mm] + loggamma;

	if (logalphacur > logalpha_next[m]) {
	  logalpha_next[m] = logalphacur;
	  prev_branch[t][m] = br;
	  prev_state[t][m] = mm;
	}
      }
    /* swap the buffers */
    tmp = logalpha_curr;
    logalpha_curr = logalpha_next;
    logalpha_next = tmp;
  }

  /* Find the ending state.
     For this purpose, the values of logalpha are modified
     with the a priori probability of the finishing states (beta_end).
     Process the corresponding probabilities of the paths              */
  logp = -HUGE_VAL;
  cur_state = -1;

  for (m = 0; m < nb_states; m++)
    if (logalpha_curr[m] + logbeta_end[m] > logp) {
      logp = logalpha_curr[m] + logbeta_end[m];
      cur_state = m;
    }

  /* Backward pass for reconstructing the sequence of branches */
  for (t = K - 1; t >= 0; t--) {
    seq[t] = prev_branch[t + 1][cur_state];
    cur_state = prev_state[t + 1][cur_state];
  }

  vec_delete (logalpha_curr);
  vec_delete (logalpha_next);
  imat_delete (prev_state);
  imat_delete (prev_branch);

  return seq;
}


/*-----------------------------------------------------------------------------*/
/* Maximum A Posteriori Sequence Decoding with side information at the decoder */
/* This algorithm is the same as viterbi, but exploit some side information 
   represented by the a priori probabilities sideinfo of each state 
   for all positions described by the vector of index sideinfo_pos. 
   The positions of the side information that are given
   in sideinfo_pos are given between 1 and the length of the sequence.  
*/
ivec viterbi_side (mat x, imat next_state, mat next_state_pt,
		   vec alpha_0, ivec sideinfo_pos, mat sideinfo)
{
  int  m, mm, br;
  /* Number of states in the automaton */
  int  nb_states = imat_width (next_state);
  /* Number of branches from a node */
  int  nb_branch = imat_height (next_state);
  int  K = mat_width (x), t, cur_state;
  int  sideinfo_idx = 0;	/* The index position in the vector of side information */
  double rt, logalphacur, logp;

  /* variables start at index 1 */
  vec  logalpha_curr = vec_new (nb_states);	/* current value of logalpha    */
  vec  logalpha_next = vec_new (nb_states);	/* next value of logalpha       */
  vec  tmp;			/* use double-buffering for logalpha    */
  double gamma;			/* gamma needs only be computed locally */
  imat prev_branch = imat_new_set (DUMM_NODE, K + 1, nb_states);
  imat prev_state = imat_new_set (DUMM_NODE, K + 1, nb_states);
  ivec seq = ivec_new (K);

  /* Viterbi initialization */
  vec_copy (logalpha_curr, alpha_0);
  vec_apply_function (logalpha_curr, IT_FUNCTION (log), NULL);

  /* Forward pass */
  for (t = 1; t <= K; t++) {
    vec_set (logalpha_next, -HUGE_VAL);

    /* Compute the value gamma. Assume the loop on m' (represented by mm) */
    for (mm = 0; mm < nb_states; mm++)
      for (br = 0; br < nb_branch; br++) {

	/* Check if the transition exists in the automaton */
	m = next_state[br][mm];
	if (m == DUMM_NODE)
	  continue;

	rt = x[br][t - 1];	/* Symbol clock differs by 1 */

	gamma = rt * next_state_pt[br][mm];

	/* Check if the value of this path is better than
	   all the ones previously tested. */
	logalphacur = logalpha_curr[mm] + log (gamma);

	if (logalphacur > logalpha_next[m]) {
	  logalpha_next[m] = logalphacur;
	  prev_branch[t][m] = br;
	  prev_state[t][m] = mm;
	}
      }

    /* Check if there is an additional side information to exploit */
    if (sideinfo_idx < ivec_length (sideinfo_pos))
      while (sideinfo_pos[sideinfo_idx] <= t) {
	for (mm = 0; mm < nb_states; mm++)
	  logalpha_next[mm] += log (sideinfo[sideinfo_idx][mm]);

	sideinfo_idx++;
	if (sideinfo_idx >= ivec_length (sideinfo_pos))
	  break;
      }

    /* swap the buffers */
    tmp = logalpha_curr;
    logalpha_curr = logalpha_next;
    logalpha_next = tmp;
  }

  /* Find the ending state.
     For this purpose, the values of logalpha are modified
     with the a priori probability of the finishing states (beta_end).
     Process the corresponding probabilities of the paths              */
  logp = -HUGE_VAL;
  cur_state = -1;

  /* Select the best terminating path */
  for (m = 0; m < nb_states; m++)
    if (logalpha_curr[m] > logp) {
      logp = logalpha_curr[m];
      cur_state = m;
    }

  /* Backward pass for reconstructing the sequence of branches */
  for (t = K - 1; t >= 0; t--) {
    seq[t] = prev_branch[t + 1][cur_state];
    cur_state = prev_state[t + 1][cur_state];
  }

  vec_delete (logalpha_curr);
  vec_delete (logalpha_next);
  imat_delete (prev_state);
  imat_delete (prev_branch);

  return seq;
}


/*---------------------------------------------------------------------------*/
/* BCJR algorithm with side information                                      */
mat bcjr_side (mat x, imat next_state, mat next_state_pt,
	       vec alpha_0, ivec sideinfo_pos, mat sideinfo)
{
  int  m, mm, br;
  /* Number of states in the automaton */
  int  nb_states = imat_width (next_state);
  /* Number of branches from a node */
  int  nb_branch = imat_height (next_state);
  int  K = mat_width (x), t;
  int  sideinfo_idx;
  double gamma;

  /* BCJR variables start at index 1 */
  mat  alpha = mat_new_zeros (K + 1, nb_states);
  mat  beta = mat_new_ones (K + 1, nb_states);
  mat  pbranch = mat_new_zeros (nb_branch, K);
  vec  pbranch_sum;

  /* Ending condition for beta */
  vec_set (beta[K], 1. / nb_states);

  /* BCJR initialization */
  vec_copy (alpha[0], alpha_0);

  /* Forward pass */
  sideinfo_idx = 0;
  for (t = 1; t <= K; t++) {
    for (br = 0; br < nb_branch; br++) {
      for (mm = 0; mm < nb_states; mm++) {
	/* Check if the transition exists in the automaton */
	m = next_state[br][mm];
	if (m == DUMM_NODE)
	  continue;

	gamma = x[br][t - 1] * next_state_pt[br][mm];	/* Symbol clock differs by 1 */
	alpha[t][m] += alpha[t - 1][mm] * gamma;
      }
    }

    /* Check if there is an additional side information to exploit */
    if (sideinfo_idx < ivec_length (sideinfo_pos))
      while (sideinfo_pos[sideinfo_idx] == t) {
	for (mm = 0; mm < nb_states; mm++)
	  alpha[t][mm] *= sideinfo[sideinfo_idx][mm];

	sideinfo_idx++;
	if (sideinfo_idx >= ivec_length (sideinfo_pos))
	  break;
      }
    vec_normalize (alpha[t], 1);
  }

  /* Backward pass */
  for (sideinfo_idx = ivec_length (sideinfo_pos) - 1;
       sideinfo_idx >= 0; sideinfo_idx--) {
    if (sideinfo_pos[sideinfo_idx] != K)
      break;

    for (mm = 0; mm < nb_states; mm++)
      beta[K][mm] *= sideinfo[sideinfo_idx][mm];
  }

  for (t = K - 1; t > 0; t--) {

    for (br = 0; br < nb_branch; br++)
      for (m = 0; m < nb_states; m++) {
	mm = next_state[br][m];
	if (mm == DUMM_NODE)
	  continue;

	/* Here, this variable is used for symbol clock t+1 */
	gamma = x[br][t] * next_state_pt[br][m];
	beta[t][m] += beta[t + 1][mm] * gamma;
      }

    /* Check if there is an additional side information to exploit */
    if (sideinfo_idx >= 0)
      while (sideinfo_pos[sideinfo_idx] == t) {
	for (mm = 0; mm < nb_states; mm++)
	  beta[t][mm] *= sideinfo[sideinfo_idx][mm];

	sideinfo_idx--;
	if (sideinfo_idx < 0)
	  break;
      }

    vec_normalize (beta[t], 1);
  }

  /* Process the metrics and project them on branchs */
  for (t = 1; t <= K; t++) {
    for (br = 0; br < nb_branch; br++)
      for (mm = 0; mm < nb_states; mm++) {
	m = next_state[br][mm];
	if (m == DUMM_NODE)
	  continue;

	gamma = x[br][t - 1] * next_state_pt[br][mm];
	pbranch[br][t - 1] += alpha[t - 1][mm] * beta[t][m] * gamma;
      }
  }

  pbranch_sum = mat_cols_sum (pbranch);

  for (t = 0; t < K; t++)
    mat_col_div_by (pbranch, t, pbranch_sum[t]);

  mat_delete (alpha);
  mat_delete (beta);
  vec_delete (pbranch_sum);

  return pbranch;
}
