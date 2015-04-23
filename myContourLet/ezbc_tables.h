/* libezbc - EZBC subband coding/decoding
      Copyright (C) 2004 Vivien Chappelier

      This library is inspired from code originaly developped by Yongjun Wu Jan in MC-EZBC. This rewritten version may still contain code similar to his original implementation.

      Reference paper: "Embedded Image Coding Using Zeroblocks Of Subband/Wavelet Coefficient And Context Modeling", Shih-Ta Hsiang and John W. Woods, Data Compression Conference (DCC '01), March 27 - 29, 2001, Snowbird, Utah 
*/

#ifndef __EZBC_TABLES
#define __EZBC_TABLES

#define SIGN_PRED_NEGATIVE   (1 << 7)
#define SIGN_PRED_POSITIVE   0
#define EZBC_SIGN_PREDICTOR(x) ((x) >> 7)
#define EZBC_SIGN_CONTEXT(x) ((x) & 0x7f)

typedef struct {
  int LIP_zc_cxts_count;         /* contexts for List of Independent Pixels */
  ezbc_byte_t *LIP_zc_cxts;      
  int node_zc_cxts_count;        /* contexts for quad-tree nodes */
  ezbc_byte_t *node_zc_cxts;     
  int jsig_zc_cxts_count[4];     /* contexts for just significant nodes */
  ezbc_byte_t *jsig_zc_cxts[4];  
  int jsig0_zc_cxts_count[4];    /* contexts for just significant nodes ?? */
  ezbc_byte_t *jsig0_zc_cxts[4];
  int LSP_zc_cxts_count;         /* contexts for List of Significant Pixels */
  ezbc_byte_t *LSP_zc_cxts;      
  int sc_cxts_count;             /* contexts for sign coding */
  ezbc_byte_t *sc_cxts;          
} ezbc_context_tables_t;

void ezbc_initialize_zc_cxts(void);
void ezbc_initialize_sc_cxts(void);

void ezbc_attach_cxts_main(ezbc_context_tables_t *tables);
void ezbc_attach_cxts_diag(ezbc_context_tables_t *tables);
void ezbc_attach_cxts_hv(ezbc_context_tables_t *tables);
void ezbc_attach_cxts_other(ezbc_context_tables_t *tables);
void ezbc_attach_cxts_none(ezbc_context_tables_t *tables);

#endif
