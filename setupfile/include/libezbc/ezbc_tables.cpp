/* libezbc - EZBC subband coding/decoding
      Copyright (C) 2004 Vivien Chappelier

      This library is inspired from code originaly developped by Yongjun Wu Jan in MC-EZBC. This rewritten version may still contain code similar to his original implementation.

      Reference paper: "Embedded Image Coding Using Zeroblocks Of Subband/Wavelet Coefficient And Context Modeling", Shih-Ta Hsiang and John W. Woods, Data Compression Conference (DCC '01), March 27 - 29, 2001, Snowbird, Utah 
*/

#include <string.h>
#include "ezbc_types.h"
#include "ezbc_tables.h"

ezbc_byte_t LIP_zc_cxts_main[ZC_MASK+1];
ezbc_byte_t node_zc_cxts_main[ZC_MASK+1];
ezbc_byte_t LSP_zc_cxts_main[ZC_MASK+1];
ezbc_byte_t LSP2_zc_cxts_main[ZC_MASK+1];
ezbc_byte_t sc_cxts_main[SC_MASK+1];

ezbc_byte_t jsig_00_zc_cxts_main[ZC_MASK+1];
ezbc_byte_t jsig_01_zc_cxts_main[ZC_MASK+1];
ezbc_byte_t jsig_10_zc_cxts_main[ZC_MASK+1];
ezbc_byte_t jsig_11_zc_cxts_main[ZC_MASK+1];
ezbc_byte_t jsig0_00_zc_cxts_main[ZC_MASK+1];
ezbc_byte_t jsig0_01_zc_cxts_main[ZC_MASK+1];
ezbc_byte_t jsig0_10_zc_cxts_main[ZC_MASK+1];
ezbc_byte_t jsig0_11_zc_cxts_main[ZC_MASK+1];

int LIP_zc_cxts_main_count = 4;
int node_zc_cxts_main_count = 10;
int LSP_zc_cxts_main_count = 5;
int sc_cxts_main_count = 9;
int LSP2_zc_cxts_main_count = 3;
int jsig_00_zc_cxts_main_count = 10;
int jsig_01_zc_cxts_main_count = 10;
int jsig_10_zc_cxts_main_count = 30;
int jsig_11_zc_cxts_main_count = 11;
int jsig0_00_zc_cxts_main_count = 10;
int jsig0_01_zc_cxts_main_count = 10;
int jsig0_10_zc_cxts_main_count = 30;
int jsig0_11_zc_cxts_main_count = 11;


ezbc_byte_t LIP_zc_cxts_diag[ZC_MASK+1];
ezbc_byte_t node_zc_cxts_diag[ZC_MASK+1];
ezbc_byte_t LSP_zc_cxts_diag[ZC_MASK+1];
ezbc_byte_t LSP2_zc_cxts_diag[ZC_MASK+1];
ezbc_byte_t sc_cxts_diag[SC_MASK+1];

ezbc_byte_t jsig_00_zc_cxts_diag[ZC_MASK+1];
ezbc_byte_t jsig_01_zc_cxts_diag[ZC_MASK+1];
ezbc_byte_t jsig_10_zc_cxts_diag[ZC_MASK+1];
ezbc_byte_t jsig_11_zc_cxts_diag[ZC_MASK+1];
ezbc_byte_t jsig0_00_zc_cxts_diag[ZC_MASK+1];
ezbc_byte_t jsig0_01_zc_cxts_diag[ZC_MASK+1];
ezbc_byte_t jsig0_10_zc_cxts_diag[ZC_MASK+1];
ezbc_byte_t jsig0_11_zc_cxts_diag[ZC_MASK+1];

int LIP_zc_cxts_diag_count = 3;
int node_zc_cxts_diag_count = 9;
int LSP_zc_cxts_diag_count = 4;
int sc_cxts_diag_count = 9;
int LSP2_zc_cxts_diag_count = 3;
int jsig_00_zc_cxts_diag_count = 10;
int jsig_01_zc_cxts_diag_count = 10;
int jsig_10_zc_cxts_diag_count = 30;
int jsig_11_zc_cxts_diag_count = 11;
int jsig0_00_zc_cxts_diag_count = 10;
int jsig0_01_zc_cxts_diag_count = 10;
int jsig0_10_zc_cxts_diag_count = 30;
int jsig0_11_zc_cxts_diag_count = 11;

ezbc_byte_t LIP_zc_cxts_hv[ZC_MASK+1];
ezbc_byte_t node_zc_cxts_hv[ZC_MASK+1];
ezbc_byte_t LSP_zc_cxts_hv[ZC_MASK+1];
ezbc_byte_t LSP2_zc_cxts_hv[ZC_MASK+1];
ezbc_byte_t sc_cxts_hv[SC_MASK+1];

ezbc_byte_t jsig_00_zc_cxts_hv[ZC_MASK+1];
ezbc_byte_t jsig_01_zc_cxts_hv[ZC_MASK+1];
ezbc_byte_t jsig_10_zc_cxts_hv[ZC_MASK+1];
ezbc_byte_t jsig_11_zc_cxts_hv[ZC_MASK+1];
ezbc_byte_t jsig0_00_zc_cxts_hv[ZC_MASK+1];
ezbc_byte_t jsig0_01_zc_cxts_hv[ZC_MASK+1];
ezbc_byte_t jsig0_10_zc_cxts_hv[ZC_MASK+1];
ezbc_byte_t jsig0_11_zc_cxts_hv[ZC_MASK+1];

int LIP_zc_cxts_hv_count = 3;
int node_zc_cxts_hv_count = 9;
int LSP_zc_cxts_hv_count = 4;
int sc_cxts_hv_count = 9;
int LSP2_zc_cxts_hv_count = 3;
int jsig_00_zc_cxts_hv_count = 10;
int jsig_01_zc_cxts_hv_count = 10;
int jsig_10_zc_cxts_hv_count = 30;
int jsig_11_zc_cxts_hv_count = 11;
int jsig0_00_zc_cxts_hv_count = 10;
int jsig0_01_zc_cxts_hv_count = 10;
int jsig0_10_zc_cxts_hv_count = 30;
int jsig0_11_zc_cxts_hv_count = 11;

ezbc_byte_t LIP_zc_cxts_other[ZC_MASK+1];
ezbc_byte_t node_zc_cxts_other[ZC_MASK+1];
ezbc_byte_t LSP_zc_cxts_other[ZC_MASK+1];
ezbc_byte_t LSP2_zc_cxts_other[ZC_MASK+1];
ezbc_byte_t sc_cxts_other[SC_MASK+1];

ezbc_byte_t jsig_00_zc_cxts_other[ZC_MASK+1];
ezbc_byte_t jsig_01_zc_cxts_other[ZC_MASK+1];
ezbc_byte_t jsig_10_zc_cxts_other[ZC_MASK+1];
ezbc_byte_t jsig_11_zc_cxts_other[ZC_MASK+1];
ezbc_byte_t jsig0_00_zc_cxts_other[ZC_MASK+1];
ezbc_byte_t jsig0_01_zc_cxts_other[ZC_MASK+1];
ezbc_byte_t jsig0_10_zc_cxts_other[ZC_MASK+1];
ezbc_byte_t jsig0_11_zc_cxts_other[ZC_MASK+1];

int LIP_zc_cxts_other_count = 3;
int node_zc_cxts_other_count = 9;
int LSP_zc_cxts_other_count = 4;
int sc_cxts_other_count = 9;
int LSP2_zc_cxts_other_count = 3;
int jsig_00_zc_cxts_other_count = 10;
int jsig_01_zc_cxts_other_count = 10;
int jsig_10_zc_cxts_other_count = 30;
int jsig_11_zc_cxts_other_count = 11;
int jsig0_00_zc_cxts_other_count = 10;
int jsig0_01_zc_cxts_other_count = 10;
int jsig0_10_zc_cxts_other_count = 30;
int jsig0_11_zc_cxts_other_count = 11;

ezbc_byte_t zc_cxts_none[ZC_MASK+1];
ezbc_byte_t sc_cxts_none[SC_MASK+1];

void ezbc_attach_cxts_main(ezbc_context_tables_t *tables)
{
  tables->LIP_zc_cxts = LIP_zc_cxts_main;
  tables->node_zc_cxts = node_zc_cxts_main;
  tables->jsig0_zc_cxts[0] = jsig0_00_zc_cxts_main;
  tables->jsig0_zc_cxts[1] = jsig0_01_zc_cxts_main;
  tables->jsig0_zc_cxts[2] = jsig0_10_zc_cxts_main;
  tables->jsig0_zc_cxts[3] = jsig0_11_zc_cxts_main;
  tables->jsig_zc_cxts[0] = jsig_00_zc_cxts_main;
  tables->jsig_zc_cxts[1] = jsig_01_zc_cxts_main;
  tables->jsig_zc_cxts[2] = jsig_10_zc_cxts_main;
  tables->jsig_zc_cxts[3] = jsig_11_zc_cxts_main;
  //  tables->LSP_zc_cxts = LSP_zc_cxts_main;
  tables->LSP_zc_cxts = LSP2_zc_cxts_main;
  tables->sc_cxts = sc_cxts_main;

  tables->LIP_zc_cxts_count = LIP_zc_cxts_main_count;
  tables->node_zc_cxts_count = node_zc_cxts_main_count;
  tables->jsig0_zc_cxts_count[0] = jsig0_00_zc_cxts_main_count;
  tables->jsig0_zc_cxts_count[1] = jsig0_01_zc_cxts_main_count;
  tables->jsig0_zc_cxts_count[2] = jsig0_10_zc_cxts_main_count;
  tables->jsig0_zc_cxts_count[3] = jsig0_11_zc_cxts_main_count;
  tables->jsig_zc_cxts_count[0] = jsig_00_zc_cxts_main_count;
  tables->jsig_zc_cxts_count[1] = jsig_01_zc_cxts_main_count;
  tables->jsig_zc_cxts_count[2] = jsig_10_zc_cxts_main_count;
  tables->jsig_zc_cxts_count[3] = jsig_11_zc_cxts_main_count;
  //  tables->LSP_zc_cxts_count = LSP_zc_cxts_main_count;
  tables->LSP_zc_cxts_count = LSP2_zc_cxts_main_count;
  tables->sc_cxts_count = sc_cxts_main_count;
}

void ezbc_attach_cxts_diag(ezbc_context_tables_t *tables)
{
  tables->LIP_zc_cxts = LIP_zc_cxts_diag;
  tables->node_zc_cxts = node_zc_cxts_diag;
  tables->jsig0_zc_cxts[0] = jsig0_00_zc_cxts_diag;
  tables->jsig0_zc_cxts[1] = jsig0_01_zc_cxts_diag;
  tables->jsig0_zc_cxts[2] = jsig0_10_zc_cxts_diag;
  tables->jsig0_zc_cxts[3] = jsig0_11_zc_cxts_diag;
  tables->jsig_zc_cxts[0] = jsig_00_zc_cxts_diag;
  tables->jsig_zc_cxts[1] = jsig_01_zc_cxts_diag;
  tables->jsig_zc_cxts[2] = jsig_10_zc_cxts_diag;
  tables->jsig_zc_cxts[3] = jsig_11_zc_cxts_diag;
  //  tables->LSP_zc_cxts = LSP_zc_cxts_diag;
  tables->LSP_zc_cxts = LSP2_zc_cxts_diag;
  tables->sc_cxts = sc_cxts_diag;

  tables->LIP_zc_cxts_count = LIP_zc_cxts_diag_count;
  tables->node_zc_cxts_count = node_zc_cxts_diag_count;
  tables->jsig0_zc_cxts_count[0] = jsig0_00_zc_cxts_diag_count;
  tables->jsig0_zc_cxts_count[1] = jsig0_01_zc_cxts_diag_count;
  tables->jsig0_zc_cxts_count[2] = jsig0_10_zc_cxts_diag_count;
  tables->jsig0_zc_cxts_count[3] = jsig0_11_zc_cxts_diag_count;
  tables->jsig_zc_cxts_count[0] = jsig_00_zc_cxts_diag_count;
  tables->jsig_zc_cxts_count[1] = jsig_01_zc_cxts_diag_count;
  tables->jsig_zc_cxts_count[2] = jsig_10_zc_cxts_diag_count;
  tables->jsig_zc_cxts_count[3] = jsig_11_zc_cxts_diag_count;
  //  tables->LSP_zc_cxts_count = LSP_zc_cxts_diag_count;
  tables->LSP_zc_cxts_count = LSP2_zc_cxts_diag_count;
  tables->sc_cxts_count = sc_cxts_diag_count;
}

void ezbc_attach_cxts_hv(ezbc_context_tables_t *tables)
{
  tables->LIP_zc_cxts = LIP_zc_cxts_hv;
  tables->node_zc_cxts = node_zc_cxts_hv;
  tables->jsig0_zc_cxts[0] = jsig0_00_zc_cxts_hv;
  tables->jsig0_zc_cxts[1] = jsig0_01_zc_cxts_hv;
  tables->jsig0_zc_cxts[2] = jsig0_10_zc_cxts_hv;
  tables->jsig0_zc_cxts[3] = jsig0_11_zc_cxts_hv;
  tables->jsig_zc_cxts[0] = jsig_00_zc_cxts_hv;
  tables->jsig_zc_cxts[1] = jsig_01_zc_cxts_hv;
  tables->jsig_zc_cxts[2] = jsig_10_zc_cxts_hv;
  tables->jsig_zc_cxts[3] = jsig_11_zc_cxts_hv;
  //  tables->LSP_zc_cxts = LSP_zc_cxts_hv;
  tables->LSP_zc_cxts = LSP2_zc_cxts_hv;
  tables->sc_cxts = sc_cxts_hv;

  tables->LIP_zc_cxts_count = LIP_zc_cxts_hv_count;
  tables->node_zc_cxts_count = node_zc_cxts_hv_count;
  tables->jsig0_zc_cxts_count[0] = jsig0_00_zc_cxts_hv_count;
  tables->jsig0_zc_cxts_count[1] = jsig0_01_zc_cxts_hv_count;
  tables->jsig0_zc_cxts_count[2] = jsig0_10_zc_cxts_hv_count;
  tables->jsig0_zc_cxts_count[3] = jsig0_11_zc_cxts_hv_count;
  tables->jsig_zc_cxts_count[0] = jsig_00_zc_cxts_hv_count;
  tables->jsig_zc_cxts_count[1] = jsig_01_zc_cxts_hv_count;
  tables->jsig_zc_cxts_count[2] = jsig_10_zc_cxts_hv_count;
  tables->jsig_zc_cxts_count[3] = jsig_11_zc_cxts_hv_count;
  //  tables->LSP_zc_cxts_count = LSP_zc_cxts_hv_count;
  tables->LSP_zc_cxts_count = LSP2_zc_cxts_hv_count;
  tables->sc_cxts_count = sc_cxts_hv_count;
}

void ezbc_attach_cxts_other(ezbc_context_tables_t *tables)
{
  tables->LIP_zc_cxts = LIP_zc_cxts_other;
  tables->node_zc_cxts = node_zc_cxts_other;
  tables->jsig0_zc_cxts[0] = jsig0_00_zc_cxts_other;
  tables->jsig0_zc_cxts[1] = jsig0_01_zc_cxts_other;
  tables->jsig0_zc_cxts[2] = jsig0_10_zc_cxts_other;
  tables->jsig0_zc_cxts[3] = jsig0_11_zc_cxts_other;
  tables->jsig_zc_cxts[0] = jsig_00_zc_cxts_other;
  tables->jsig_zc_cxts[1] = jsig_01_zc_cxts_other;
  tables->jsig_zc_cxts[2] = jsig_10_zc_cxts_other;
  tables->jsig_zc_cxts[3] = jsig_11_zc_cxts_other;
  //  tables->LSP_zc_cxts = LSP_zc_cxts_other;
  tables->LSP_zc_cxts = LSP2_zc_cxts_other;
  tables->sc_cxts = sc_cxts_other;

  tables->LIP_zc_cxts_count = LIP_zc_cxts_other_count;
  tables->node_zc_cxts_count = node_zc_cxts_other_count;
  tables->jsig0_zc_cxts_count[0] = jsig0_00_zc_cxts_other_count;
  tables->jsig0_zc_cxts_count[1] = jsig0_01_zc_cxts_other_count;
  tables->jsig0_zc_cxts_count[2] = jsig0_10_zc_cxts_other_count;
  tables->jsig0_zc_cxts_count[3] = jsig0_11_zc_cxts_other_count;
  tables->jsig_zc_cxts_count[0] = jsig_00_zc_cxts_other_count;
  tables->jsig_zc_cxts_count[1] = jsig_01_zc_cxts_other_count;
  tables->jsig_zc_cxts_count[2] = jsig_10_zc_cxts_other_count;
  tables->jsig_zc_cxts_count[3] = jsig_11_zc_cxts_other_count;
  //  tables->LSP_zc_cxts_count = LSP_zc_cxts_other_count;
  tables->LSP_zc_cxts_count = LSP2_zc_cxts_other_count;
  tables->sc_cxts_count = sc_cxts_other_count;
}

void ezbc_attach_cxts_none(ezbc_context_tables_t *tables)
{
  memset(zc_cxts_none, 0, (ZC_MASK+1) * sizeof(ezbc_byte_t));
  memset(sc_cxts_none, 0, (SC_MASK+1) * sizeof(ezbc_byte_t));

  tables->LIP_zc_cxts = zc_cxts_none;
  tables->node_zc_cxts = zc_cxts_none;
  tables->jsig0_zc_cxts[0] = zc_cxts_none;
  tables->jsig0_zc_cxts[1] = zc_cxts_none;
  tables->jsig0_zc_cxts[2] = zc_cxts_none;
  tables->jsig0_zc_cxts[3] = zc_cxts_none;
  tables->jsig_zc_cxts[0] = zc_cxts_none;
  tables->jsig_zc_cxts[1] = zc_cxts_none;
  tables->jsig_zc_cxts[2] = zc_cxts_none;
  tables->jsig_zc_cxts[3] = zc_cxts_none;
  //  tables->LSP_zc_cxts = LSP_zc_cxts_none;
  tables->LSP_zc_cxts = zc_cxts_none;
  tables->sc_cxts = sc_cxts_none;

  tables->LIP_zc_cxts_count = 1;
  tables->node_zc_cxts_count = 1;
  tables->jsig0_zc_cxts_count[0] = 1;
  tables->jsig0_zc_cxts_count[1] = 1;
  tables->jsig0_zc_cxts_count[2] = 1;
  tables->jsig0_zc_cxts_count[3] = 1;
  tables->jsig_zc_cxts_count[0] = 1;
  tables->jsig_zc_cxts_count[1] = 1;
  tables->jsig_zc_cxts_count[2] = 1;
  tables->jsig_zc_cxts_count[3] = 1;
  //  tables->LSP_zc_cxts_count = LSP_zc_cxts_other_count;
  tables->LSP_zc_cxts_count = 1;
  tables->sc_cxts_count = 1;
}

static void initialize_LIP_zc_cxts(void)
{
  int index;
  int h, v, d;
  ezbc_byte_t main_context, diag_context, hv_context, other_context;

  for(index = 0; index <= ZC_MASK; index++) {
      /* First, form the context map for the horizontal and vertical bands.
         These both use the same context map, because the horizontally
         high-pass band is physically transposed before encoding. */

      /* d v d */
      /* h   h */
      /* d v d */

      h = ((index >> CL_POS) & 1) + ((index >> CR_POS) & 1);
      v = ((index >> TC_POS) & 1) + ((index >> BC_POS) & 1);
      d = ((index >> TL_POS) & 1) + ((index >> TR_POS) & 1) +
          ((index >> BL_POS) & 1) + ((index >> BR_POS) & 1);

      if ((h == 0) && (v == 0) && (d == 1)) {
        main_context = 0; diag_context = 0;
      }	else if ((h == 0) && (v == 0)) {
        main_context = 1; diag_context = 1;
      } else if ((h == 0) && (v == 1))  {
        main_context = 1; diag_context = 1;
      } else if ((h == 1) && (v == 0)) {
        main_context = 2; diag_context = 1;
      } else if(((h == 0) && (v == 2)) || ((h == 2) && (v == 0))){
        main_context = 3; diag_context = 2;
      } else if((h+v) > 2) {
        main_context = 3; diag_context = 2;
      } else { 
        main_context = 3; diag_context = 2;
      }

      LIP_zc_cxts_main[index] = main_context;
      LIP_zc_cxts_diag[index] = diag_context;

      /* hv context */
      if (d == 0 && h+v == 1) {
        hv_context = 0;
      }	else if (d == 0) {
        hv_context = 1;
      } else if (d == 1) {
        hv_context = 1;
      } else { 
        hv_context = 2;
      }
      LIP_zc_cxts_hv[index] = hv_context;

      /* other context */
      other_context = 0;//diag_context;
      LIP_zc_cxts_other[index] = other_context;
  }
}

static void initialize_node_zc_cxts(void)
{
  int index;
  int h, v, dp, dm, d, p;
  ezbc_byte_t main_context, diag_context, hv_context, other_context;

  for(index = 0; index <= ZC_MASK; index++) {

    p = ((index>>PA_POS)&1);
    h = ((index>>CL_POS)&1) + ((index>>CR_POS)&1);
    v = ((index>>TC_POS)&1) + ((index>>BC_POS)&1);
    dm = ((index>>TL_POS)&1) + ((index>>BR_POS)&1);
    dp = ((index>>TR_POS)&1) + ((index>>BL_POS)&1);
    d = dp + dm;

    if(p) { /* parent exists */
      if(h == 0 && v == 0 && d == 1) {
        main_context = 0; diag_context = 0;
      } else if(h == 0 && v == 0) {
        main_context = 1; diag_context = 1;
      } else if(h == 0 && v == 1) {
        main_context = 1; diag_context = (d < 3) ? 1 : 2;
      } else if(h == 1 && v == 0){
        main_context = 2; diag_context = (d < 3) ? 1 : 2;
      } else if((h == 0 && v == 2) || (h == 2 && v == 0)) {
        main_context = 3; diag_context = 2;
      } else if((h + v) > 2) {
        main_context = 4; diag_context = 3;
      } else {
        main_context = 3; diag_context = (d < 2) ? 1 : 2;
      }
      
      /* hv context */
      if(d == 0 && h+v == 1) {
        hv_context = 0;
      } else if(d == 0) {
        hv_context = 1;
      } else if(d == 1) {
        hv_context = (h+v < 3) ? 1 : 2;
      } else if(dm == 1 && dp == 1) {
        hv_context = (h+v < 2) ? 1 : 2;
      } else if(d == 2) {
        hv_context = 2;
      } else {
        hv_context = 3;
      }
    } else {
      if(h == 0 && v == 0) {
        main_context = 5; diag_context = 4;
      } else if(h + v == 1) {
        main_context = (d < 3) ? 5 : 6; diag_context = (d < 3) ? 4 : 5;
      } else if((h == 0 && v == 2) ||(v == 2 && h == 0)) {
        main_context = 6; diag_context = 5;
      } else if(h + v > 2) {
        main_context = 7; diag_context = 6;
      } else {
        main_context = 6; diag_context = 5;
      }

      /* hv context */
      if(d == 0) {
        hv_context = 4;
      } else if(d == 1) {
        hv_context = (h+v < 3) ? 4 : 5;
      } else if(d == 2) {
        hv_context = 5;
      } else {
        hv_context = 6;
      }
    }

    /* other context */
    other_context = p?1:0;
    
    node_zc_cxts_main[index] = main_context;
    node_zc_cxts_diag[index] = diag_context;
    node_zc_cxts_hv[index] = hv_context;
    node_zc_cxts_other[index] = other_context;
  }
}

void ezbc_initialize_sc_cxts(void)
{
  int index;
  int vpos, vneg, hpos, hneg, h, v;
  ezbc_byte_t main_context, main_predict;
  ezbc_byte_t diag_context, diag_predict;
  ezbc_byte_t hv_context, hv_predict;
  ezbc_byte_t other_context, other_predict;
  int nwpos, nwneg, nepos, neneg, nw, ne;

  for(index = 0; index <= SC_MASK; index++) {

    vpos = (index >> V_PVE_BIT_POS) & 1;
    vneg = (index >> V_NVE_BIT_POS) & 1;
    hpos = (index >> H_PVE_BIT_POS) & 1;
    hneg = (index >> H_NVE_BIT_POS) & 1;

    h = hpos - hneg;
    v = vpos - vneg;

    if(v > 0) {
      diag_predict = main_predict = SIGN_PRED_NEGATIVE;
    } else {
      diag_predict = main_predict = SIGN_PRED_POSITIVE;
      if(v < 0) {
	h = -h;
	v = -v;
      }
    }

    if(v == 0) {
      if(h < 0) {
	main_predict =  SIGN_PRED_NEGATIVE;
	main_context = 1;
	diag_context = 1;
      } else if(h > 0) {
        diag_predict = SIGN_PRED_NEGATIVE;
        main_context = 1;
	diag_context = 1;
      } else {
        main_context = 0;
	diag_context = 0;
      }
    } else {
      main_context = 3 + h;
      diag_context = 3 + h;
    }

    if(main_context == 0) { /* v == 0 && h == 0 */
      nwpos = (index >> NW_PVE_BIT_POS) & 1;
      nwneg = (index >> NW_NVE_BIT_POS) & 1;
      nepos = (index >> NE_PVE_BIT_POS) & 1;
      neneg = (index >> NE_NVE_BIT_POS) & 1;

      nw = nwpos - nwneg;
      ne = nepos - neneg;

      if(nw < 0){
        nw = -nw;
	ne = -ne;
        diag_predict = main_predict = SIGN_PRED_NEGATIVE;
      }

      if(nw == 0) {
        if(ne < 0) {
          diag_predict = SIGN_PRED_NEGATIVE;
          main_context = 5;
	  diag_context = 5;
        } else if(ne > 0) {
          main_predict = SIGN_PRED_NEGATIVE;
          main_context = 5;
	  diag_context = 5;
        } else {
          main_context = 0;
	  diag_context = 0;
        }
      } else {
        main_context = 5;
	diag_context = 5;
      }
    }

    sc_cxts_main[index] = main_context | main_predict;
    sc_cxts_diag[index] = diag_context | diag_predict;
  
    /* hv context */
    nwpos = (index >> NW_PVE_BIT_POS) & 1;
    nwneg = (index >> NW_NVE_BIT_POS) & 1;
    nepos = (index >> NE_PVE_BIT_POS) & 1;
    neneg = (index >> NE_NVE_BIT_POS) & 1;
    
    nw = nwpos - nwneg;
    ne = nepos - neneg;

    if(ne > 0) {
      hv_predict = SIGN_PRED_NEGATIVE;
    } else {
      hv_predict = SIGN_PRED_POSITIVE;
      if(ne < 0) {
	nw = -nw;
	ne = -ne;
      }
    }

    if(ne == 0) {
      if(nw < 0) {
	hv_context = 1;
      } else if(nw > 0) {
        hv_predict = SIGN_PRED_NEGATIVE;
	hv_context = 1;
      } else {
	hv_context = 0;
      }
    } else {
      hv_context = 3 + nw;
    }

    if(ne == 0 && nw == 0) {
      vpos = (index >> V_PVE_BIT_POS) & 1;
      vneg = (index >> V_NVE_BIT_POS) & 1;
      hpos = (index >> H_PVE_BIT_POS) & 1;
      hneg = (index >> H_NVE_BIT_POS) & 1;

      h = hpos - hneg;
      v = vpos - vneg;

      if(h < 0) {
        h = -h;
	v = -v;
        hv_predict = SIGN_PRED_NEGATIVE;
      }

      if(h == 0) {
        if(v < 0) {
          hv_predict = SIGN_PRED_NEGATIVE;
	  hv_context = 5;
        } else if(v > 0) {
	  hv_context = 5;
        } else {
	  hv_context = 0;
        }
      } else {
	hv_context = 5;
      }
    }

    sc_cxts_hv[index] = hv_context | hv_predict;

    /* other context */
    other_context = 0;
    other_predict = SIGN_PRED_POSITIVE;
    sc_cxts_other[index] = other_context | other_predict;
  }
}

static void initialize_LSP_zc_cxts(void)
{
  int index;
  int h, v, d, dm, dp;
  ezbc_byte_t main_context, diag_context, hv_context, other_context;

  for(index = 0; index <= ZC_MASK; index++) {
      /* First, form the context map for the horizontal and vertical bands.
         These both use the same context map, because the horizontally
         high-pass band is physically transposed before encoding. */

      h = ((index>>CL_POS)&1) + ((index>>CR_POS)&1);
      v = ((index>>TC_POS)&1) + ((index>>BC_POS)&1);
      dm = ((index>>TL_POS)&1) + ((index>>BR_POS)&1);
      dp = ((index>>TR_POS)&1) + ((index>>BL_POS)&1);

      if((h == 0) && (v == 0)) {
        main_context = 0; diag_context = 0;
      } else if((h == 0)) {
        main_context = 1; diag_context = 1;
      } else if((v == 0)) {
        main_context = 2; diag_context = 1;
      } else if((h == 2) && (v == 2)) {
        main_context = 4; diag_context = 3;
      } else {
        main_context = 3; diag_context = 2;
      }

      /* hv context */
      if(d == 0) {
        hv_context = 0;
      } else if(dm == 0 || dp == 0) {
        hv_context = 1;
      } else if(dm == 2 && dp == 2) {
        hv_context = 3;
      } else {
        hv_context = 2;
      }

      LSP_zc_cxts_hv[index] = hv_context;


      /* other context */
      other_context = 0;//diag_context;
      LSP_zc_cxts_other[index] = other_context;
  }
}

static void initialize_LSP2_zc_cxts(void)
{
  int index;
  ezbc_byte_t main_context, diag_context, hv_context, other_context;
  int h, v, dm, dp, hf, vf;

  for(index = 0; index <= ZC_MASK; index++) {
      /* First, form the context map for the horizontal and vertical bands.
         These both use the same context map, because the horizontally
         high-pass band is physically transposed before encoding. */

      hf = ((index>>CL2_POS)&1) + ((index>>CR2_POS)&1);
      vf = ((index>>TC2_POS)&1) + ((index>>BC2_POS)&1);

      h = ((index>>CL_POS)&1) + ((index>>CR_POS)&1);
      v = ((index>>TC_POS)&1) + ((index>>BC_POS)&1);
      dm = ((index>>TL_POS)&1) + ((index>>BR_POS)&1);
      dp = ((index>>TR_POS)&1) + ((index>>BL_POS)&1);

      main_context = 0; diag_context = 0;

      if((h == 0) && (v == 0)) {
        main_context = 0; diag_context = 0;
      } else if((hf == 0) && (vf == 0) && (dm+dp < 2)) {
        main_context = 0; diag_context = 0;
      } else if((hf == 0) && (vf == 0)) {
        main_context = 0; diag_context = 0;
      } else if((hf == 0)) {
        main_context = 2; diag_context = 1;
      } else if((vf == 0)) {
        main_context = 2; diag_context = 1;
      } else if ((hf == 2) && (vf == 2)) {
        main_context = 1; diag_context = 1;
      } else if ((hf == 2)) {
        main_context = 1; diag_context = 1;
      } else if ((vf == 2)) {
        main_context = 1; diag_context = 1;
      } else {
        main_context = 2; diag_context = 1;
      }

      LSP2_zc_cxts_main[index] = main_context;
      LSP2_zc_cxts_diag[index] = diag_context;

      /* hv context: TEMP // BORROWED FROM LSP */
      if(dm+dp == 0) {
        hv_context = 0;
      } else if(dm == 0 || dp == 0) {
        hv_context = 1;
      } else if(dm == 2 && dp == 2) {
        hv_context = 3;
      } else {
        hv_context = 2;
      }

      LSP_zc_cxts_hv[index] = hv_context;

      /* other context */
      other_context = 0;
      LSP_zc_cxts_other[index] = other_context;
  }
}

static void initialize_jsig_zc_cxts(void)
{
  int index;
  ezbc_byte_t main_context, diag_context;
  int h, v, d, p, l, r, t, bc, tr;

  for(index = 0; index <= ZC_MASK; index++) {

    p = ((index>>PA_POS)&1);
    l = ((index>>CL_POS)&1);
    t = ((index>>TC_POS)&1);
    d = ((index>>TL_POS)&1) + ((index>>TR_POS)&1) + ((index>>BL_POS)&1);

    if(p == 0) {
      if(l == 1 && t == 1) {
        main_context = 3; diag_context = 3;
      } else if(l == 1 || t == 1) {
        main_context = 1; diag_context = 1;
      } else {
        main_context = 0; diag_context = 0;
      }
    } else {
      if(l == 1 && t == 1) {
        main_context = 3; diag_context = 3;
      } else if(l == 1 || t == 1) {
        main_context = 2; diag_context = 2;
      } else if(d > 1) {
        main_context = 1; diag_context = 2;
      } else {
        main_context = 1; diag_context = 1;
      }
    }

    jsig_00_zc_cxts_main[index] = main_context;
    jsig_00_zc_cxts_diag[index] = diag_context;
    jsig_00_zc_cxts_hv[index] = main_context;
    jsig_00_zc_cxts_other[index] = p?1:0;

    r = ((index>>CR_POS)&1);
    d = ((index>>TL_POS)&1) + ((index>>TR_POS)&1) + ((index>>BR_POS)&1);

    if(l == 0) {//sig
      if(p == 0) {
	if(r == 1 && t == 1) {
	  main_context = 2; diag_context = 2;
	} else if(r == 1 || t == 1){
	  main_context = 1; diag_context = 1;
	} else if(d > 1){
	  main_context = 0; diag_context = 1; 
	} else {
	  main_context = 0; diag_context = 0;
	}
      } else {
	if(r || t) {
	  main_context = 2; diag_context = 2;
	} else if(d > 1) {
	  main_context = 1; diag_context = 2;
	} else {
	  main_context = 1; diag_context = 1;
	}
      }
    } else {
      if(p == 0) {
	if(r == 1 && t == 1) {
	  main_context = 7; diag_context = 7;
	} else if(r || t) {
	  main_context = 4; diag_context = 4;
	} else if(d > 1) {
	  main_context = 3; diag_context = 4; 
	} else {
	  main_context = 3; diag_context = 3;
	}
      }else{
	if(r && t) {
	  main_context = 7; diag_context = 7;
	} else if(r || t) {
	  main_context = 6; diag_context = 6;
	} else if(d > 1) {
	  main_context = 5; diag_context = 6;
	} else {
	  main_context = 5; diag_context = 5;
	}
      }
    }

    jsig_01_zc_cxts_main[index] = main_context;
    jsig_01_zc_cxts_diag[index] = diag_context;
    jsig_01_zc_cxts_hv[index] = main_context;
    jsig_01_zc_cxts_other[index] = p?1:0;

    bc = ((index>>BC_POS)&1);
    tr = ((index>>TR_POS)&1);
    v = ((index>>TC_POS)&1) + ((index>>BC_POS)&1);
    d = ((index>>TL_POS)&1) + ((index>>TR_POS)&1) +
        ((index>>BL_POS)&1) + ((index>>BR_POS)&1);

    if((t == 0) && (tr == 0)) {
      if(p == 0) {
	if(bc && l) {
	  main_context = 3; diag_context = 3;
	} else if(bc || l) {
	  main_context = 1; diag_context = 2;
	} else if(d) {
	  main_context = 0; diag_context = (d > 1) ? 1 : 0;
	} else {
	  main_context = 0; diag_context = 0;
	}
      } else {
	if(l) {
	  main_context = 3; diag_context = 3;
	} else if(bc) {
	  main_context = 2; diag_context = 3;
	} else if(d > 1) {
	  main_context = 1; diag_context = 3;
	} else {
	  main_context = 1; diag_context = 2;
	}
      }
    } else {
      if(p == 0){
	if(t) {
	  if(l && bc) {
	    main_context = 8; diag_context = 9;
	  } else if(l || bc){
	    main_context = 7; diag_context = 8;
	  } else if(d) {
	    main_context = 5; diag_context = (d > 1)? 8 : 7;
	  } else{
	    main_context = 5; diag_context = 5;
	  }
	}else{
	  if(l || bc) {
	    main_context = 6; diag_context = 6;
	  } else if(d > 1) {
	    main_context = 5; diag_context = 6;
	  } else{
	    main_context = 4; diag_context = 4;
	  }
	}
      } else {
	if(t) {
	  if(l && bc) {
	    main_context = 12; diag_context = 13;
	  } else if(l || bc) {
	    main_context = 11; diag_context = 12;
	  } else if(d) {
	    main_context = 9; diag_context = (d > 1)? 12 : 10;
	  } else {
	    main_context = 9; diag_context = 10;
	  }
	}else{
	  if(l && bc) {
	    main_context = 12; diag_context = 13;
	  } else if(l || bc) {
	    main_context = 10; diag_context = 11;
	  } else if(d > 1) {
	    main_context = 9; diag_context = 11;
	  } else {
	    main_context = 9; diag_context = 10;
	  }
	}
      }
    }

    jsig_10_zc_cxts_main[index] = main_context;
    jsig_10_zc_cxts_diag[index] = diag_context;
    jsig_10_zc_cxts_hv[index] = main_context;
    jsig_10_zc_cxts_other[index] = p?1:0;

    h = ((index>>CL_POS)&1) + ((index>>CR_POS)&1);

    if(p) {
      if((h == 0) && (v == 0) && (d == 1)) {
	main_context = 0; diag_context = 0;
      } else if((h == 0) && (v == 0)) {
	main_context = 1; diag_context = 1;
      } else if((h == 0) && (v == 1)) {
	main_context = 1; diag_context = 1;
      } else if((h == 1) && (v == 0)) {
	main_context = 2; diag_context = 1;
      } else if((h + v) > 2) {
	main_context = 4; diag_context = 3;
      } else {//(v1 == 1) && (v2 == 1)
	main_context = 3; diag_context = 2;
      }
    } else {
      if((h == 0) && (v == 0) && (d == 1)) {
	main_context = 5; diag_context = 4;
      } else if((h == 0) && (v == 0)) {
	main_context = 6; diag_context = 5;
      } else if((h == 0) && (v == 1)) {
        main_context = 6; diag_context = 5;
      } else if((h == 1) && (v == 0)) {
        main_context = 6; diag_context = 5;
      } else if((h + v) > 2) {
        main_context = 8; diag_context = 7;
      } else {//(h == 1) && (v == 1)
        main_context = 7; diag_context = 6;
      }
    }

    jsig_11_zc_cxts_main[index] = main_context;
    jsig_11_zc_cxts_diag[index] = diag_context;
    jsig_11_zc_cxts_hv[index] = main_context;
    jsig_11_zc_cxts_other[index] = p?1:0;
  }
}

static void initialize_jsig0_00_zc_cxts(void)
{
  int index;
  ezbc_byte_t main_context, diag_context;
  int p, h, v, d;

  for(index = 0; index <= ZC_MASK; index++) {

    p = ((index>>PA_POS)&1);
    h = ((index>>CL_POS)&1);
    v = ((index>>TC_POS)&1);
    d = ((index>>TL_POS)&1) + ((index>>TR_POS)&1) + ((index>>BL_POS)&1);

    if(h && v){
      main_context = 5; diag_context = 3;
    }else if(h){
      main_context = (p) ? 4 : 2; diag_context = 2;
    }else if(v){
      main_context = (p) ? 3 : 2; diag_context = 2;
    }else if(d > 1){
      main_context = 1; diag_context = 2;
    }else if(d){
      main_context = 1; diag_context = 1;
    }else{
      main_context = 0; diag_context = 0;
    }

    jsig0_00_zc_cxts_main[index] = main_context;
    jsig0_00_zc_cxts_diag[index] = diag_context;
    jsig0_00_zc_cxts_hv[index] = main_context;
    jsig0_00_zc_cxts_other[index] = p?1:0;
  }
}

static void initialize_jsig0_01_zc_cxts(void)
{
  int index;
  ezbc_byte_t main_context, diag_context;
  int p, l, r, t, d;

  for(index = 0; index <= ZC_MASK; index++) {
    p = ((index>>PA_POS)&1);
    l = ((index>>CL_POS)&1);
    r = ((index>>CR_POS)&1);
    t = ((index>>TC_POS)&1);
    d = ((index>>TL_POS)&1) + ((index>>TR_POS)&1) + ((index>>BR_POS)&1);

    if(l == 0) {
      if(r && t) {
        main_context = 4; diag_context = 3;
      } else if(r) {
        main_context = (p) ? 3 : 2; diag_context = 2;
      } else if(t) {
        main_context = (p) ? 3 : 2; diag_context = 2;
      } else if(d > 1) {
        main_context = 1; diag_context = 2;
      } else if(d) {
        main_context = 1; diag_context = 1;
      } else {
        main_context = 0; diag_context = 0;
      }
    } else {
      if(r && t){
        main_context = 7; diag_context = 6;
      } else if(r || t) {
        main_context = (p) ? 7 : 6;
        diag_context = (p) ? 6 : 5;
      } else if(d > 1) {
        main_context = 5; diag_context = (p) ? 6 : 5;
      } else {
        main_context = 5; diag_context = 4;
      }
    }

    jsig0_01_zc_cxts_main[index] = main_context;
    jsig0_01_zc_cxts_diag[index] = diag_context;
    jsig0_01_zc_cxts_hv[index] = main_context;
    jsig0_01_zc_cxts_other[index] = p?1:0;
  }
}

static void initialize_jsig0_10_zc_cxts(void)
{
  int index;
  ezbc_byte_t main_context, diag_context;
  int p, l, t, tr, bc, v, d;

  for(index = 0; index <= ZC_MASK; index++) {
 
    p = ((index>>PA_POS)&1);
    l = ((index>>CL_POS)&1);
    t = ((index>>TC_POS)&1);
    bc = ((index>>BC_POS)&1);
    tr = ((index>>TR_POS)&1);

    v = ((index>>TC_POS)&1) + ((index>>BC_POS)&1);
    d = ((index>>TL_POS)&1) + ((index>>TR_POS)&1) +
        ((index>>BL_POS)&1) + ((index>>BR_POS)&1);
 
    if((t == 0) && (tr == 0)){
      if(p == 0){
	if(bc && l) {
	  main_context = 5; diag_context = 4;
	} else if(bc || l){
	  main_context = 1; diag_context = 1;
	} else if(d) {
	  main_context = 2; diag_context = (d > 1)? 1 : 2;
	} else {
	  main_context = 0; diag_context = 0;
	}
      }else{
	if(bc && l){
	  main_context = 5; diag_context = 4;
	} else if(l){
	  main_context = 4; diag_context = 3;
	} else if(bc){
	  main_context = 3; diag_context = 3;
	} else if(d > 1){
	  main_context = 2; diag_context = 3;
	} else if(d){
	  main_context = 2; diag_context = 2;
	} else {
	  main_context = 0; diag_context = 0;
	}
      }
    } else {
      if(p == 0) {
	if(t) {
	  if(l && bc) {
	    main_context = 19; diag_context = 18;
	  } else if(l) {
	    main_context = 17; diag_context = 17;	    	    
	  } else if(bc) {
	    main_context = 17; diag_context = 17;
	  } else if(d) {
	    main_context = 15; diag_context = (d > 1) ? 17 : 15;
	  } else{
	    main_context = 15; diag_context = 14;	    
	  }
	} else {
	  if(l && bc) {
	    main_context = 12; diag_context = 12;
	  } else if(l) {
	    main_context = 12; diag_context = 12;
	  } else if(bc) {
	    main_context = 12; diag_context = 12;
	  } else if(d > 1) {
	    main_context = 15; diag_context = 12;
	  } else {
	    main_context = 10; diag_context = 10;
	  }
	}
      } else {
	if(t) {
	  if(l && bc) {
	    main_context = 24; diag_context = 23;
	  } else if(l) {
	    main_context = 27; diag_context = 27;
	  } else if(bc){
	    main_context = 27; diag_context = 27;
	  } else if(d){
	    main_context = 25; diag_context = (d > 1)? 27 : 24;
	  } else{
	    main_context = 25; diag_context = 24;
          }
        }else{
          if(l && bc) {
            main_context = 24; diag_context = 23;
          } else if(l) {
            main_context = 22; diag_context = 22;
          } else if(bc) {
            main_context = 22; diag_context = 22;
          } else if(d > 1) {
            main_context = 25; diag_context = (d == 2)? 22 : 22;
          } else {
            main_context = 25; diag_context = 24;
          }
        }
      }
    }

    jsig0_10_zc_cxts_main[index] = main_context;
    jsig0_10_zc_cxts_diag[index] = diag_context;
    jsig0_10_zc_cxts_hv[index] = main_context;
    jsig0_10_zc_cxts_other[index] = p?1:0;
  }
}

static void initialize_jsig0_11_zc_cxts(void)
{
  int index;
  ezbc_byte_t main_context, diag_context;
  int p, h, v, d;

  for(index = 0; index <= ZC_MASK; index++) {

    p = ((index>>PA_POS)&1);
    h = ((index>>CL_POS)&1) + ((index>>CR_POS)&1);
    v = ((index>>TC_POS)&1) + ((index>>BC_POS)&1);
    d = ((index>>TL_POS)&1) + ((index>>TR_POS)&1) +
        ((index>>BL_POS)&1) + ((index>>BR_POS)&1);

    if((h == 0) && (v == 0) && (d == 1)) {
      main_context = 0; diag_context = 0;
    } else if((h == 0) && (v == 0)) {
      main_context = 1; diag_context = 1;
    } else if((h == 0) && (v == 1)) {
      main_context = 1; diag_context = 1;
    } else if((h == 1) && (v == 0)) {
      main_context = 2; diag_context = 1;
    } else if((h == 0) && (v == 2)) {
      main_context = 3; diag_context = 2;
    } else if((h == 2) && (v == 0)) {
      main_context = 3; diag_context = 2;
    } else if((h+v) > 2) {
      main_context = 4; diag_context = 3;
    } else {//(h == 1) && (v == 1)
      main_context = 3; diag_context = 2;
    }

    jsig0_11_zc_cxts_main[index] = main_context;
    jsig0_11_zc_cxts_diag[index] = diag_context;
    jsig0_11_zc_cxts_hv[index] = main_context;
    jsig0_11_zc_cxts_other[index] = p?1:0;
  }
}

/* initialize zero coding lookup tables */
void ezbc_initialize_zc_cxts(void)
{
  initialize_LIP_zc_cxts();

  initialize_node_zc_cxts();

  initialize_LSP_zc_cxts();

  initialize_LSP2_zc_cxts();

  initialize_jsig_zc_cxts();

  initialize_jsig0_00_zc_cxts();
  initialize_jsig0_01_zc_cxts();
  initialize_jsig0_10_zc_cxts();
  initialize_jsig0_11_zc_cxts();
}
