/* libezbc - EZBC subband coding/decoding
      Copyright (C) 2004 Vivien Chappelier

      This library is inspired from code originaly developped by Yongjun Wu Jan in MC-EZBC. This rewritten version may still contain code similar to his original implementation.

      Reference paper: "Embedded Image Coding Using Zeroblocks Of Subband/Wavelet Coefficient And Context Modeling", Shih-Ta Hsiang and John W. Woods, Data Compression Conference (DCC '01), March 27 - 29, 2001, Snowbird, Utah 
*/

#ifndef __EZBC_ENCODER_H
#define __EZBC_ENCODER_H

#include "ezbc_codec.h"

/* encode the subband using the EZBC algorithm */
int ezbc_subband_encode(ezbc_subband_codec_t *codec, int bit, int *truncation_list);

#endif
