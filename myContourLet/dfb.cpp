/*
   contourlet - Implementation of the contourlet transform for image coding
   Copyright (C) 2005 Vivien Chappelier - IRISA/University of Rennes 1

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public
   License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU General Public
   License along with this program; if not, write to the Free
   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
#include "stdafx.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "it/mat.h"
#include "it/io.h"
#include "dfb.h"

#define OLD_COORD

#define max(x,y) (((x)>(y))?(x):(y))

#define Q(y,x) Q[((y)<<1)+x]
#define R(y,x) R[((y)<<1)+x]
#define B(y,x) B[((y)<<1)+x]
//fq add begin*************************************************************
/* set some elements of the matrix to value val*/
/*#define mat_set_between( m, r1, c1, r2, c2, val ) do {\ 
    int ce = c2;\
    int re = r2;\
    int x, y;   \
    assert( m );  \
    mat_end_col_param( m, ce );     \
    mat_end_row_param( m, re );     \
    for( y = r1 ; y <= re ; y++ )   \
for( x = c1 ; x <= ce ; x++ ) \
  m[ y ][ x ] = val;    \
} while( 0 );*/
/* to use the entity end in a mat function, this macro is called */
#define mat_set_between( m, r1, c1, r2, c2, val ) do {\
	int ce = c2;\
	int re = r2;\
    int x, y;   \
    assert( m );  \
	MAT_END_COL_PARAM( m, ce );     \
    MAT_END_ROW_PARAM( m, re );     \
	for( y = r1 ; y <= re ; y++ )   \
		for( x = c1 ; x <= ce ; x++ ) \
			m[ y ][ x ] = val;    \
	}while(0);

static inline bmat bmat_set_between (bmat m, int r1, int c1, int r2,
 int c2, byte val) {
    mat_set_between (m, r1, c1, r2, c2, val);
    return m;
}
//fq add end*************************************************************/
static int Q0[2][2] = {
  {  1, -1 },
  {  1,  1 },
};
static int Q1[2][2] = {
  {  1,  1 },
  { -1,  1 },
};

static int R0[2][2] = {
  {  1,  1 },
  {  0,  1 },
};
static int R1[2][2] = {
  {  1,  -1 },
  {  0,  1 },
};
static int R2[2][2] = {
  {  1,  0 },
  {  1,  1 },
};
static int R3[2][2] = {
  {  1,  0 },
  { -1,  1 },
};

static int I[2][2] = {
  {  1,  0 },
  {  0,  1 },
};



#if 1

#define N 6
static double const B[N] =
  { 0.6290, 0.1914, 0.0958, 0.0520, 0.0300, 0.0114 };

#else

#define N 16

static double const B[N] = {
  0.6337676,
  0.2037873,
  0.1137364,
  0.0728123,
  0.0488452,
  0.0331130,
  0.0222485,
  0.0146264,
  0.0093114,
  0.0056841,
  0.0032908,
  0.0017822,
  0.0008854,
  0.0003910,
  0.0001448,
  0.0000397
};
#endif

#define symetrize(n, rl, rh) table_symetrize[N+rl-1][N+rh-1][N+n-1]
static int table_symetrize[2*N][2*N][2*N];

#define ODD_SYMETRY
#ifdef ODD_SYMETRY

/*
[ . . . a b c . . ]
[ B C B a b c B A ]
*/

void generate_symetrize()
{
  int n, rl, rh;

  for(rh = -N+1; rh <= N; rh++)
    for(rl = -N+1; rl <= N; rl++)
      for(n = -N+1; n <= N; n++) {
	if(rh == rl)
	  symetrize(n, rl, rh) = rh;
	else
	  symetrize(n, rl, rh) = rl+abs(abs(((n) + rh - 2*rl - 1) % (2*rh - 2*rl)) - (rh - rl - 1));
      }
}
#else

/*
[ . . . a b c . . ]
[ C B A a b c C B ]
*/

void generate_symetrize()
{
  int n, rl, rh;

  for(rh = -N+1; rh <= N; rh++)
    for(rl = -N+1; rl <= N; rl++)
      for(n = -N+1; n <= N; n++) {
	if(rh == rl)
	  symetrize(n, rl, rh) = rh;
	else {
	  int n_, rh_, rl_;

	  rl_ = rl + N - 1;
	  rh_ = rh + N - 1;
	  if(rh_ < rl_) { n_ = rl_; rl_ = rh_; rh_ = n_; }
	  n_ = n + N - 1;

	  while(n_ < rl_ || n_ > rh_) {
	    n_ -= rl_;
	    if(n_ < 0) n_ = ~n_;
	    n_  = -n_;
	    n_  += rh_ - rl_;
	    if(n_ < 0) n_ = ~n_;
	    n_  -= rh_ - rl_;
	    n_  = -n_;
	    n_  += rl_;
	  }
	  n_ -= N - 1;
	  symetrize(n, rl, rh) = n_;
	}
      }
}
#endif

#define FILTER(sys, src, X, Y)						\
do {									\
  o = 0; v = 1;                                                         \
  do {                                                                  \
    /* find valid range */						\
    rl = rh = 2*N;							\
    for(n = N-1; n > -N; n--) {						\
      n_ = n+1;								\
      k = X;								\
      l = Y;								\
      if(is_inside(sys, k, l)) { rh = n+1; break; }	                \
    }									\
    for(n = N-1; n > -N; n--) {						\
      n_ = -n;								\
      k = X;								\
      l = Y;								\
      if(is_inside(sys, k, l)) { rl = -n; break; }	                \
    }									\
    if(rl == 2*N && rh == 2*N) switch(o) {                              \
	  case 0: o = +2; break;                                        \
	  case 2: o = -2; break;                                        \
	  default: v = 0; fprintf(stderr, "Fatal: could not find an interior point for filtering at (%d,%d)\n", X, Y);  /*(*(char*)0=0); */break;\
    }                                                                   \
  } while(rl == 2*N && rh == 2*N && v);                                      \
  if(!v) break; \
  if(rl == 2*N) break; \
  if(rh == 2*N) break; \
   v = 0;								\
  for(n = 0; n < N; n++) {						\
    u = 0;								\
    n_ = symetrize(n+1, rl, rh);					\
    k = X;								\
    l = Y;								\
    u += (src)[l][k];							\
    									\
    n_ = symetrize(-n, rl, rh);						\
    k = X;								\
    l = Y;								\
    u -= (src)[l][k];							\
    									\
    v += u*B[n];							\
  }									\
									\
} while(0)

/*----------------------------- geometry ------------------------------------*/

/* tell if the point x,y is inside the mask */
static inline int is_inside(bmat mask, int x, int y)
{
  int w, h;

  if(x < 0) return(0);
  if(y < 0) return(0);
  w = bmat_width(mask);
  if(x >= w) return(0);
  h = bmat_height(mask);
  if(y >= h) return(0);
  return(mask[y][x]);
}

/* find the top right corner of a rectangular mask */
void find_origin(bmat mask, int *ox, int *oy)
{
  int i;
  int w, h;
  int l;
  
  w = bmat_width(mask);
  h = bmat_height(mask);
  l = (max(w, h) + 1) & (~1);

  if(!is_inside(mask, l/2, l/2))
    fprintf(stderr, "Fatal: central point not in mask\n");

  /* assumes rectangular region */
  for(i = l/2; i >= 0 && is_inside(mask, i, l/2); i--);
  *ox = i + 1;
  for(i = l/2; i >= 0 && is_inside(mask, l/2, i); i--);
  *oy = i + 1;
}

/* find the bottom left corner of a rectangular mask */
void find_bound(bmat mask, int *rx, int *ry)
{
  int i;
  int w, h;
  int l;
  
  w = bmat_width(mask);
  h = bmat_height(mask);
  l = (max(w, h) + 1) & (~1);

  if(!is_inside(mask, l/2, l/2))
    fprintf(stderr, "Fatal: central point not in mask\n");

  /* assumes rectangular region */
  for(i = l/2; i < l && is_inside(mask, i, l/2); i++);
  *rx = i;
  for(i = l/2; i < l && is_inside(mask, l/2, i); i++);
  *ry = i;
}

/* resample the subband to a rectangular support region */
/* input: subband to resample          */
/* tmp: temporary buffer               */
/* B: resampling matrix                */
/* w: width                            */
/* h: height                           */
void resample_mask(bmat input, bmat output, int *B, int w, int h)
{
  int x, y, X, Y;
  int l = (max(w,h) + 1) & ~1;

  for(y = 0; y < h; y ++) {
    for(x = 0; x < w; x ++) {
      X = B(0,0)*(x - l/2) + B(0,1)*(y - l/2) + l/2;
      Y = B(1,0)*(x - l/2) + B(1,1)*(y - l/2) + l/2;
      if(X >= 0 && X < w && Y >= 0 && Y < h)
	output[y][x] = input[Y][X];
      else
	output[y][x] = 0;
    }
  }
}

void forward_mask(bmat hmask,
		  bmat vmask,
		  bmat mask,
		  int *Q, int *R,
		  int w, int h)
{
  bmat tmp;
  int x, y;
  int X, Y;
  int l;

  w = bmat_width(mask);
  h = bmat_height(mask);
  l = (max(w, h) + 1) & ~1;
  tmp = bmat_new_zeros(l, l);

  /* rotate input */
  for(y = 0; y < h; y++) {
    for(x = 0; x < w; x++) {
#ifdef OLD_COORD
      X =  R(0,0)*x + R(0,1)*y - R(0,1)*l/2;
      Y =  R(1,0)*x + R(1,1)*y - R(1,0)*l/2;
#else
      X =  R(0,0)*(x-l/2) + R(0,1)*(y-l/2) + l/2;
      Y =  R(1,0)*(x-l/2) + R(1,1)*(y-l/2) + l/2;
#endif
      if(X >= 0 && X < w &&
	 Y >= 0 && Y < h)
	tmp[y][x] = mask[Y][X];
    }
  }

  /* split/pack hband/vband */
  bmat_zeros(vmask);
  bmat_zeros(hmask);

  for(y = 0; y < h; y ++) {
    for(x = 0; x < w; x ++) {
#ifdef OLD_COORD
      X = (Q(0,0)*x + Q(1,0)*y + (Q(0,0)+Q(0,1))*l/2)/2;
      Y = (Q(0,1)*x + Q(1,1)*y + (Q(1,0)+Q(1,1))*l/2)/2;
#else
      X = ( Q(1,1)*(x-l/2) - Q(0,1)*(y-l/2)) / 2 + l/2;
      Y = (-Q(1,0)*(x-l/2) + Q(0,0)*(y-l/2)) / 2 + l/2;
#endif
      // TEMP
      //      fprintf(stderr, "x,y=%d,%d X,Y=%d,%d\n", x, y, X, Y);
      if(!((x + y) & 1))
        hmask[Y][X] = tmp[y][x];
      else
        vmask[Y][X] = tmp[y][x];
    }
  }

  bmat_delete(tmp);
}

/*-------------------------- sampling ---------------------------------------*/

void forward_update(mat hband,
		    mat vband,
		    mat input,
		    int *Q, int *R,
		    int w, int h)
{
  mat tmp = hband;
  int x, y;
  int X, Y;
  int l = (max(w,h) + 1) & ~1;

  /* rotate input */
  mat_zeros(tmp);
  for(y = 0; y < h; y++) {
    for(x = 0; x < w; x++) {
#ifdef OLD_COORD
      X =  R(0,0)*x + R(0,1)*y - R(0,1)*l/2;
      Y =  R(1,0)*x + R(1,1)*y - R(1,0)*l/2;
#else
      X =  R(0,0)*(x-l/2) + R(0,1)*(y-l/2) + l/2;
      Y =  R(1,0)*(x-l/2) + R(1,1)*(y-l/2) + l/2;
#endif
      if(X >= 0 && X < w &&
	 Y >= 0 && Y < h)
      tmp[y][x] = input[Y][X];
    }
  }

  for(y = 0; y < h; y ++) {
    for(x = 0; x < w; x ++) {
      input[y][x] = tmp[y][x];
    }
  }

  /* split/pack hband/vband */
  for(y = 0; y < h; y ++) {
    for(x = 0; x < w; x ++) {
#ifdef OLD_COORD
      X = (Q(0,0)*x + Q(1,0)*y + (Q(0,0)+Q(0,1))*l/2)/2;
      Y = (Q(0,1)*x + Q(1,1)*y + (Q(1,0)+Q(1,1))*l/2)/2;
#else
      X = ( Q(1,1)*(x-l/2) - Q(0,1)*(y-l/2)) / 2 + l/2;
      Y = (-Q(1,0)*(x-l/2) + Q(0,0)*(y-l/2)) / 2 + l/2;
#endif

      if(!((x + y) & 1))
        hband[Y][X] = input[y][x];
      else
        vband[Y][X] = input[y][x];
    }
  }
}

void backward_update(mat hband,
		     mat vband,
		     mat input,
		     int *Q, int *R,
		     int w, int h)
{
  mat tmp = hband;
  int x, y;
  int X, Y;
  int l = (max(w,h) + 1) & ~1;

  /* unpack hband/vband */
  for(y = 0; y < h; y ++) {
    for(x = 0; x < w; x ++) {

#ifdef OLD_COORD
      X = (Q(0,0)*x + Q(1,0)*y + (Q(0,0)+Q(0,1))*l/2)/2;
      Y = (Q(0,1)*x + Q(1,1)*y + (Q(1,0)+Q(1,1))*l/2)/2;
#else
      X = ( Q(1,1)*(x-l/2) - Q(0,1)*(y-l/2)) / 2 + l/2;
      Y = (-Q(1,0)*(x-l/2) + Q(0,0)*(y-l/2)) / 2 + l/2;
#endif
      if(!((x + y) & 1))
        input[y][x] = hband[Y][X];
      else {
	//	fprintf(stderr, "%d,%d %d,%d\n", x, y, X, Y);
        input[y][x] = vband[Y][X];
      }
    }
  }

  for(y = 0; y < h; y ++) {
    for(x = 0; x < w; x ++) {
      tmp[y][x] = input[y][x];
    }
  }

  /* rotate output */
  mat_zeros(input);
  for(y = 0; y < h; y ++) {
    for(x = 0; x < w; x ++) {
#ifdef OLD_COORD
      X = R(0,0)*x - R(0,1)*y + R(0,1)*l/2;
      Y =-R(1,0)*x + R(1,1)*y + R(1,0)*l/2;
#else
      X = R(1,1)*(x-l/2) - R(0,1)*(y-l/2) + l/2;
      Y =-R(1,0)*(x-l/2) + R(0,0)*(y-l/2) + l/2;
#endif
      if(X >= 0 && X < w &&
	 Y >= 0 && Y < h)
	input[y][x] = tmp[Y][X];
    }
  }

}

/*------------------------- filtering ---------------------------------------*/

void direct_poly_filter_downsample_pass_A(mat tmp,
					  dfb_subband_t *hband,
					  dfb_subband_t *vband,
					  int w, int h)
{
  double v, u;
  int x, y, k, l, n;
  int rl, rh, o, n_;

  /* separable filtering for channel A */
  for(y = 0; y < h; y ++) {
    for(x = 0; x < w; x ++) {
      if(is_inside(hband->mask, x, y)) {
	FILTER(hband->mask, hband->image, x + n_, y + o);
	tmp[y][x] = v;
      }
    }
  }

  for(y = 0; y < h; y ++) {
    for(x = 0; x < w; x ++) {
      if(is_inside(vband->mask, x, y)) {
	FILTER(hband->mask, tmp, x + o, y + n_);
	vband->image[y][x] += v;
      }
    }
  }
}

void direct_poly_filter_downsample_pass_B(mat tmp,
					  dfb_subband_t *hband,
					  dfb_subband_t *vband,
					  int w, int h)
{
  double v, u;
  int x, y, k, l, n;
  int rl, rh, n_, o;

  /* separable filtering for channel B */
  for(y = 0; y < h; y ++) {
    for(x = 0; x < w; x ++) {
      if(is_inside(vband->mask, x, y)) {
	FILTER(vband->mask, vband->image, x + n_, y + o);
	tmp[y][x] = v*0.5;
      }
    }
  }

  for(y = 0; y < h; y ++) {
    for(x = 0; x < w; x ++) {
      if(is_inside(hband->mask, x, y)) {
	FILTER(vband->mask, tmp, x-1 + o, y-1 + n_);
	hband->image[y][x] -= v;
      }
    }
  }
}


/* polyphase directional filter bank analysis */
void direct_poly_filter_downsample(dfb_subband_t *hband,
				   dfb_subband_t *vband,
				   dfb_subband_t *input,
				   int *Q, int *R,
				   int w, int h)
{
  mat tmp;

  /* x(n0,n1)  -> !M -> (z0z1)^-N -+--------------------> A */
  /* |                             |                |       */
  /* |                         B(z0)B(z1)       -B(z0)B(z1) */
  /* |                             |                |       */
  /* |-> z0^-1 -> !M ---------------->(z0z1)^-2N+1 -+---> B */

  forward_update(hband->image, vband->image, input->image, Q, R, w, h);

  tmp = input->image;
#ifndef BYPASS
  /* separable filtering for channel A */
  direct_poly_filter_downsample_pass_A(tmp, hband, vband, w, h);

  /* separable filtering for channel B */
  direct_poly_filter_downsample_pass_B(tmp, hband, vband, w, h);
#endif
}

void direct_poly_interpolate_upsample_pass_A(mat tmp,
					     dfb_subband_t *hband,
					     dfb_subband_t *vband,
					     int w, int h)
{
  double v, u;
  int x, y, k, l, n;
  int rl, rh, o, n_;

  /* separable filtering for channel A */
  for(y = 0; y < h; y ++) {
    for(x = 0; x < w; x ++) {
      if(is_inside(vband->mask, x, y)) {
	FILTER(vband->mask, vband->image, x + n_, y + o);
	tmp[y][x] = v*0.5;
      }
    }
  }

  for(y = 0; y < h; y ++) {
    for(x = 0; x < w; x ++) {
      if(is_inside(hband->mask, x, y)) {
	FILTER(vband->mask, tmp, x-1 + o, y-1 + n_);
	hband->image[y][x] += v;
      }
    }
  }
}

void direct_poly_interpolate_upsample_pass_B(mat tmp,
					     dfb_subband_t *hband,
					     dfb_subband_t *vband,
					     int w, int h)
{
  double v, u;
  int x, y, k, l, p, n;
  int rl, rh, n_, o;

  p = w;

  /* separable filtering for channel A */
  for(y = 0; y < h; y ++) {
    for(x = 0; x < w; x ++) {
      if(is_inside(hband->mask, x, y)) {
	FILTER(hband->mask, hband->image, x + n_, y + o);
	tmp[y][x] = v;
      }
    }
  }

  for(y = 0; y < h; y ++) {
    for(x = 0; x < w; x ++) {
      if(is_inside(vband->mask, x, y)) {
	FILTER(hband->mask, tmp, x + o, y + n_);
	vband->image[y][x] -= v;
      }
    }
  }
}

/* polyphase directional filter bank synthesis */
void direct_poly_interpolate_upsample(dfb_subband_t *hband,
				      dfb_subband_t *vband,
				      dfb_subband_t *input,
				      int *Q, int *R,
				      int w, int h)
{
  mat tmp;

  /* separable filtering for channel A */
#ifndef BYPASS
  tmp = input->image;
  direct_poly_interpolate_upsample_pass_A(tmp, hband, vband, w, h);

  /* separable filtering for channel B */
  tmp = input->image;
  direct_poly_interpolate_upsample_pass_B(tmp, hband, vband, w, h);
#endif
  backward_update(hband->image, vband->image, input->image, Q, R, w, h);
}

/* resample the subband to a rectangular support region */
/* input: subband to resample          */
/* tmp: temporary buffer               */
/* B: resampling matrix                */
/* w: width                            */
/* h: height                           */
void forward_resample(dfb_subband_t *input, mat tmp, int *B, int w, int h)
{
  int x, y, X, Y;
  mat band = input->image;
  int l = (max(w,h) + 1) & ~1;

  /* resample */
  mat_copy(tmp, band);

  for(y = 0; y < h; y ++) {
    for(x = 0; x < w; x ++) {
      X = B(0,0)*(x - l/2) + B(0,1)*(y - l/2) + l/2;
      Y = B(1,0)*(x - l/2) + B(1,1)*(y - l/2) + l/2;
      if(X >= 0 && X < l && Y >= 0 && Y < l)
	band[y][x] = tmp[Y][X];
    }
  }
}

void backward_resample(dfb_subband_t *input, mat tmp, int *B, int w, int h)
{
  int x, y, X, Y;
  mat band = input->image;
  int l = (max(w,h) + 1) & ~1;

  /* resample */
  mat_copy(tmp, band);

  for(y = 0; y < h; y ++) {
    for(x = 0; x < w; x ++) {
      X = B(0,0)*(x - l/2) + B(0,1)*(y - l/2) + l/2;
      Y = B(1,0)*(x - l/2) + B(1,1)*(y - l/2) + l/2;
      if(X >= 0 && X < l && Y >= 0 && Y < l)
	band[Y][X] = tmp[y][x];
    }
  }
}

void dfb_init()
{
  generate_symetrize();
}

dfb_subband_t *dfb_subband_new()
{
  dfb_subband_t *subband;

  subband = (dfb_subband_t *) malloc(sizeof(dfb_subband_t));
  subband->image = NULL;
  subband->mask = NULL;
  memset(subband->Q, 0, sizeof(subband->Q));
  memset(subband->R, 0, sizeof(subband->R));
  return(subband);
}

void dfb_subband_delete(dfb_subband_t *subband)
{
  bmat_delete(subband->mask);
  free(subband);
}

/* create the structure for a new directional filterbank */
/* width: width of the image(s) to transform             */
/* height: height of the image(s) to transform           */
/* levels: (maximum) number of decomposition levels      */
dfb_t *dfb_new(int width, int height, int levels)
{
  dfb_t *dfb;
  int i;
  int w, h, l;
  int B[2][2];
  int dfb_levels;
  int dl;
  dfb_subband_t *hband, *vband, *pband;
  bmat mask;

  dfb_init();
  
  dfb = (dfb_t *) malloc(sizeof(dfb_t));

  dfb->width = w = width;
  dfb->height = h = height;
  dfb->levels = dfb_levels = levels;
  l = (max(w, h) + 1) & ~1; /* round up to nearest multiple of two */

  /* allocate a temporary mask */
  mask = bmat_new(l, l);

  /* allocate dfb binary subband tree */
  dfb->bands = (dfb_subband_t **)
    calloc(2 << dfb_levels, sizeof(dfb_subband_t *));

  /* allocate subband origins */
  dfb->ox = ivec_new(1 << dfb_levels);
  dfb->oy = ivec_new(1 << dfb_levels);
  dfb->rx = ivec_new(1 << dfb_levels);
  dfb->ry = ivec_new(1 << dfb_levels);

  /* start with a full mask */
  dfb_band(dfb, 0, 0) = dfb_subband_new();
  dfb_band(dfb, 0, 0)->image = mat_new(l, l);
  dfb_band(dfb, 0, 0)->mask = bmat_new_zeros(l, l);
  bmat_set_between(dfb_band(dfb, 0, 0)->mask,
		   (l-h) / 2, (l-w) / 2, (l+h) / 2 - 1, (l+w) / 2 - 1, 1);

  /* directional decomposition */
  for(dl = 0; dl < dfb_levels; dl++) {
    for(i = 0; i  < (1 << dl); i++) {

      pband = dfb_band(dfb, dl, i);

#define assign_band_matrixes(b, q, r) \
      do { memcpy(b->R, r, sizeof(r)); memcpy(b->Q, q, sizeof(q)); } while(0)

      if(dl > 1) {
	if(i >= (1 << (dl-1))) {
	  if(i & 1) assign_band_matrixes(pband, Q1, R2);
	  else      assign_band_matrixes(pband, Q0, R3);
	} else {
	  if(i & 1) assign_band_matrixes(pband, Q0, R0);
	  else      assign_band_matrixes(pband, Q1, R1);
	}
      } else {
	if(dl & 1) assign_band_matrixes(pband, Q1, I);
	else       assign_band_matrixes(pband, Q0, I);
      }

      /* split tree */
      hband = dfb_band(dfb, dl+1, 2*i) = dfb_subband_new();
      vband = dfb_band(dfb, dl+1, 2*i+1) = dfb_subband_new();

      /* allocate the subbands */
      hband->image = mat_new(l, l);
      vband->image = mat_new(l, l);
      hband->mask = bmat_new(l, l);
      vband->mask = bmat_new(l, l);

      /* resampling */
      if(dl == dfb_levels - 1) {
	int c = i & ((1 << (dfb_levels - 2)) - 1);

	B[0][0] = B[1][1] = 1;
	B[0][1] = B[1][0] = 0;

	if(dfb_levels > 1) {
	  if(i >= (1 << (dfb_levels - 2)))
	    B[0][1] = -(c << 1) + (1 << (dfb_levels - 2)) - 1;
	  else
	    B[1][0] = -(c << 1) + (1 << (dfb_levels - 2)) - 1;
	}

	assign_band_matrixes(hband, I, B);
	assign_band_matrixes(vband, I, B);	 
      }

      /* compute mask */
      forward_mask(hband->mask, vband->mask, pband->mask,
		   (int *)pband->Q, (int *)pband->R, w, h);

      if(dl == levels - 1) {
	/* resample the subband mask to a rectangular region */
	resample_mask(hband->mask, mask, (int *) hband->R, w, h);

	/* find the origin */
	find_origin(mask, &dfb->ox[2*i+0], &dfb->oy[2*i+0]);
	find_bound(mask, &dfb->rx[2*i+0], &dfb->ry[2*i+0]);
	// TEMP
	//	fprintf(stderr, "rect %d,%d - %d,%d\n",
	//		dfb->ox[2*i+0], dfb->oy[2*i+0],
	//		dfb->rx[2*i+0], dfb->ry[2*i+0]);

	/* resample the subband mask to a rectangular region */
	resample_mask(vband->mask, mask, (int *) vband->R, w, h);

	/* find the origin */
	find_origin(mask, &dfb->ox[2*i+1], &dfb->oy[2*i+1]);
	find_bound(mask, &dfb->rx[2*i+1], &dfb->ry[2*i+1]);
	// TEMP
	//	fprintf(stderr, "rect %d,%d - %d,%d\n",
	//		dfb->ox[2*i+1], dfb->oy[2*i+1],
	//		dfb->rx[2*i+1], dfb->ry[2*i+1]);
      }
    }
  }

  /* assign filter norms */
  dfb_levels = dfb->levels;
  for(dl = 0; dl < dfb_levels; dl++) {
    for(i = 0; i  < (1 << dl); i++) {
      hband = dfb_band(dfb, dl+1, 2*i);
      vband = dfb_band(dfb, dl+1, 2*i+1);
    }
  }

  return(dfb);
}

void dfb_delete(dfb_t *dfb)
{
  int i;
  int dl;
  dfb_subband_t *hband, *vband;

  /* free band images */
  for(dl = dfb->levels-1; dl >= 0; dl--) {
    for(i = 0; i  < (1 << dl); i++) {
      hband = dfb_band(dfb, dl+1, 2*i);
      vband = dfb_band(dfb, dl+1, 2*i+1);
      mat_delete(hband->image);
      mat_delete(vband->image);
    }
  }

  ivec_delete(dfb->ox);
  ivec_delete(dfb->oy);

  for(i = 1; i  < (2 << dfb->levels); i++)
    dfb_subband_delete(dfb->bands[i]);

  free(dfb->bands);
}


void dfb_transform(dfb_t *dfb, mat *bands, mat image, int levels)
{
  int w = mat_width(image);
  int h = mat_height(image);
  int l = (max(w, h) + 1) & (~1);
  dfb_subband_t *hband, *vband, *pband;
  int dl, i;

  mat_set_submatrix(dfb_band(dfb, 0, 0)->image, image, (l-h) / 2, (l-w) / 2);

  for(dl = 0; dl < levels; dl++) {
    for(i = 0; i  < (1 << dl); i++) {
      /* split tree */
      pband = dfb_band(dfb, dl, i);
      hband = dfb_band(dfb, dl+1, 2*i);
      vband = dfb_band(dfb, dl+1, 2*i+1);

      /* compute subbands */
      direct_poly_filter_downsample(hband, vband, pband,
				    (int *) pband->Q, (int *) pband->R, l, l);

      if(dl == levels - 1) {
	/* resample and normalize the subband */
	forward_resample(hband, pband->image, (int *) hband->R, l, l);
	forward_resample(vband, pband->image, (int *) vband->R, l, l);
      }
    }
  }

  for(i = 0; i < (1 << levels); i++) {
    int x, y, ox, oy;
    int ws, hs;
    mat src, dst;

    ox = dfb->ox[i];
    oy = dfb->oy[i];
    ws = dfb->rx[i] - ox;
    hs = dfb->ry[i] - oy;
    //TEMP
    //    fprintf(stderr, "extracting %d,%d - %d,%d\n", ox, oy, ox+ws, oy+hs);

    bands[i] = mat_new(hs, ws);
    src = dfb_band(dfb, levels, i)->image;
    dst = bands[i];
    for(y = 0; y < hs; y++) {
      for(x = 0; x < ws; x++) {
	dst[y][x] = src[oy + y][ox + x];// * normalization[levels-1][i];
      }
    }
  }
}

void dfb_itransform(dfb_t *dfb, mat *bands, mat image, int levels)
{
  int w = mat_width(image);
  int h = mat_height(image);
  int l = (max(w, h) + 1) & (~1);
  dfb_subband_t *hband, *vband, *pband;
  int dl, i;
  mat tmp;

  for(i = 0; i < (1 << levels); i++) {
    int x, y, ox, oy;
    int ws, hs;
    mat src, dst;

    ox = dfb->ox[i];
    oy = dfb->oy[i];
    ws = dfb->rx[i] - ox;
    hs = dfb->ry[i] - oy;
    dst = dfb_band(dfb, levels, i)->image;
    src = bands[i];
    for(y = 0; y < hs; y++) {
      for(x = 0; x < ws; x++) {
	dst[oy + y][ox + x] = src[y][x];// / normalization[levels-1][i];
      }
    }
    mat_delete(bands[i]);
    bands[i] = NULL;
  }

  for(dl = dfb->levels-1; dl >= 0; dl--) {
    for(i = 0; i  < (1 << dl); i++) {

      pband = dfb_band(dfb, dl, i);
      hband = dfb_band(dfb, dl+1, 2*i);
      vband = dfb_band(dfb, dl+1, 2*i+1);

      if(dl == dfb->levels - 1) {
	backward_resample(hband, pband->image, (int *) hband->R, l, l);
	backward_resample(vband, pband->image, (int *) vband->R, l, l);
      }

      direct_poly_interpolate_upsample(hband, vband, pband, (int *)pband->Q, (int *)pband->R, l, l);
    }
  }

  tmp = mat_get_submatrix(dfb_band(dfb, 0, 0)->image,
			  (l-h) / 2, (l-w) / 2,
			  (l+h) / 2 - 1, (l+w) / 2 - 1);
  mat_copy(image, tmp);
  mat_delete(tmp);
}

void dfb_flatten(mat *high, mat image, int levels)
{
  int h = mat_height(image);
  int w = mat_width(image);
  int ws, hs;
  int X, Y, x, y;
  int i, j, n;

  /* flatten high bands */
  ws = w;
  hs = h;
  n = (1 << levels);

  for(i = 0; i < n; i++) {
    ws = mat_width(high[i]);
    hs = mat_height(high[i]);
    // TEMP
    //    fprintf(stderr, "ws = %d hs = %d\n", ws, hs);

    if(i < n / 2) {
      X = 0;
      for(j = 0; j < i; j++)
	X += mat_width(high[j]);
      Y = 0;
    } else {
      Y = mat_height(high[0]);
      if(i < 3 * n / 4) {
	X = 0;
	for(j = n / 2; j < i; j++)
	  Y += mat_height(high[j]);
      } else {
	for(j = 3 * n / 4; j < i; j++)
	  Y += mat_height(high[j]);
	X = mat_width(high[n / 2]);
      }
    }

    for(y = 0; y < hs; y++)
      for(x = 0; x < ws; x++)
	image[Y+y][X+x] = high[i][y][x];
  }
}
