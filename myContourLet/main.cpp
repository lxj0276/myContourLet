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
#include "cv.h"
#include "highgui.h"
#include <time.h>
#include "it/io.h"
#include "it/distance.h"
#include <stdio.h>
#include <math.h>
#include "it/wavelet2D.h"
#include "it/mat.h"
#include "contourlet.h"
#include "dfb.h"
#include "ezbc.h"


#define BUFFER_SIZE (1*1024*1024)

// 9/7 contourlet low subband norm [level]
double norm_low[6] = {
  1.000000,
  0.982948,
  1.030575,
  1.051979,
  1.058014,
  1.058312
};

// 9/7 contourlet high subbands norms [level][dfb_levels][subband]
double norm_high[6][5][16] = {

  // DFB
  {
    {1.000000},
    {1.338955, 0.768281},
    {1.788734, 1.031742, 1.031699, 0.588007},
    {2.350204, 1.388625, 1.473061, 0.755227, 1.521047, 0.718018, 0.760509, 0.466449},
    {2.990107, 1.859578, 2.009466, 0.993439, 2.153701, 1.040220, 1.071638, 0.565028,
     2.310007, 1.015735, 1.043946, 0.511974, 1.108749, 0.539535, 0.580279, 0.383226}
  },

  // Highest frequencies
  {
    {0.759782},
    {1.068118, 0.710115},
    {1.557636, 0.922336, 0.885625, 0.513870},
    {2.066849, 1.199964, 1.312638, 0.679034, 1.314328, 0.611044, 0.667065, 0.406727},
    {2.591734, 1.650462, 1.726335, 0.866216, 1.933047, 0.919827, 0.946988, 0.519211,
     2.004114, 0.871491, 0.880225, 0.441114, 0.979973, 0.469240, 0.495562, 0.338999}
  },

  {
    {0.709848},
    {1.006673, 0.691288},
    {1.505208, 0.880912, 0.857108, 0.490243},
    {2.004624, 1.154857, 1.268061, 0.637940, 1.281535, 0.585212, 0.641248, 0.383733},
    {2.461666, 1.619596, 1.693704, 0.813626, 1.870790, 0.884116, 0.917510, 0.470455,
     1.949415, 0.851772, 0.859531, 0.412020, 0.943970, 0.448550, 0.476647, 0.313942}
  },

  {
    {0.753806},
    {1.067996, 0.730151},
    {1.591337, 0.929780, 0.908668, 0.518957},
    {2.142259, 1.210500, 1.324341, 0.681005, 1.369928, 0.613891, 0.670235, 0.410106},
    {2.646395, 1.724389, 1.779619, 0.850579, 1.943830, 0.929273, 0.975566, 0.504108,
     2.089051, 0.907208, 0.903755, 0.431151, 0.981607, 0.471800, 0.507208, 0.336407}
  },

  {
    {0.775910},
    {1.098225, 0.747825},
    {1.631322, 0.953360, 0.932293, 0.532706},
    {2.204255, 1.236997, 1.351154, 0.702051, 1.409039, 0.627751, 0.684257, 0.422790},
    {2.703611, 1.755877, 1.805906, 0.865382, 1.966184, 0.947038, 0.995714, 0.521175,
     2.134306, 0.924937, 0.917630, 0.439011, 0.993215, 0.481139, 0.517305, 0.348460}
  },

  // Lowest frequencies
  {
    {0.782607},
    {1.107343, 0.753079},
    {1.643244, 0.959912, 0.939323, 0.536721},
    {2.203433, 1.234178, 1.349443, 0.704742, 1.411630, 0.628675, 0.683316, 0.427423},
    {2.846164, 2.127241, 1.898790, 1.298897, 2.080579, 1.299276, 1.085920, 0.808484,
     2.180500, 1.169607, 0.982991, 0.594266, 1.034878, 0.563380, 0.533165, 0.393168}
  }
};

mat load_mat(IplImage* img)

{

       int i,j;

       int width=img->width;//矩阵列数为图像宽

       int height=img->height;  //矩阵行数为图像高

       mat m = mat_new (height, width);

       for(i=0;i<height;i++)

       {

              for(j=0;j<width;j++)

              {

                     m[i][j]=cvGetReal2D(img, i, j);

              }

       }

       return m;

}

IplImage* save_mat(mat m)
{ 
	int i,j;
	int width=mat_height(m); //矩阵列数为图像宽
	int height=mat_width(m); //矩阵行数为图像高

	IplImage* img = cvCreateImage(cvSize(height,width),8,1);

	for(i=0;i<height;i++)
	{
		for(j=0;j<width;j++)
			{
				cvSetReal2D(img, i,j,m[i][j] );
			}
	}
		return img;
}

     



int main()
{
  int i, j;
  int w, h;
  int levels, ct_levels, wt_levels,level_init;
  double rate,rate_init;
  ivec dfb_levels;
  mat source, dest;
  contourlet_t *contourlet;
  mat wavelet;
  int length;
  unsigned char *buffer;

  //初始化参数
  int argc=6;
      rate_init=2;
	  level_init=5;


#define LEVELS 5
#define IMPULSE 100.
	  //读入图片
	 IplImage *img=cvLoadImage("E:\\C练习\\Contourlet_use\\Contourlet_demo_src\\1.pgm");
	  
	 //灰度化
	 IplImage* src = cvCreateImage(cvGetSize(img),img->depth,1);
	  cvCvtColor(img,src,CV_BGR2GRAY);	 
	  //归一化
	  cvNormalize(src,src,1,0,CV_MINMAX);
	 
	  //改变尺寸
	    CvSize dst_size;
		dst_size.height = 128;
		dst_size.width = 128;
		IplImage* dst= cvCreateImage( dst_size,8, 1 );
		cvResize(src,dst,CV_INTER_LINEAR );
  
//IplImage转矩阵mat
		 source=load_mat(dst);
 
// source = mat_pgm_read("1.pgm");//程序本来的读入图片方式

//IplImage* high_img=save_mat(source);//矩阵转图片
//cvSaveImage("outP.jpg",high_img); //保存

  h = mat_height(source);
  w = mat_width(source);
  dest = mat_new(w, h);
  rate = rate_init * w * h;
  levels = level_init;//5
  ct_levels = argc - 5;              /* contourlet levels 1    argc=6 */
  wt_levels = 3;    /* wavelet levels 3 */
  dfb_levels = ivec_new(ct_levels);
  for(i = 0; i < ct_levels; i++)
    dfb_levels[i] = 2;//每层的方向数


  buffer = bvec_new_zeros(BUFFER_SIZE);

  contourlet = contourlet_new(ct_levels, dfb_levels);
  contourlet->wt_levels = wt_levels;
//Contourlet分解
  contourlet_transform(contourlet, source);
  //小波分解
  wavelet = it_dwt2D(contourlet->low, it_wavelet_lifting_97, wt_levels);
  contourlet->dwt = it_wavelet2D_split(wavelet, wt_levels);

  /* normalize the subbands */
  for(i = 0; i < ct_levels; i++)
    for(j = 0; j < (1 << dfb_levels[i]); j++)
      mat_mul_by(contourlet->high[i][j], norm_high[1+i][dfb_levels[i]][j]);
  mat_mul_by(contourlet->low, norm_low[ct_levels]);

  /* make flat images */
  mat_pgm_write("dwt.pgm", wavelet);
  for(i = 0; i < ct_levels; i++) {
    char filename[256];

    mat dfb_rec = mat_new((h >> i) + 1, (w >> i) + 1);
    if(dfb_levels[i])
      dfb_flatten(contourlet->high[i], dfb_rec, dfb_levels[i]);
    else
      mat_set_submatrix(dfb_rec, contourlet->high[i][0], 0, 0);
    mat_incr(dfb_rec, 128);
    sprintf(filename, "dfb%d.pgm", i);
    mat_pgm_write(filename, dfb_rec);
    mat_decr(dfb_rec, 128);
    mat_delete(dfb_rec);
  }

  /* EZBC encoding */
  length = ezbc_encode(contourlet, buffer, BUFFER_SIZE, rate);

  /* EZBC decoding */
  ezbc_decode(contourlet, buffer, BUFFER_SIZE, rate);

  mat_pgm_write("rec_low.pgm", contourlet->dwt[0]);

  /* make flat images */
  for(i = 0; i < ct_levels; i++) {
    char filename[256];

    mat dfb_rec = mat_new((h >> i) + 1, (w >> i) + 1);
    if(dfb_levels[i])
      dfb_flatten(contourlet->high[i], dfb_rec, dfb_levels[i]);
    else
      mat_set_submatrix(dfb_rec, contourlet->high[i][0], 0, 0);
    mat_incr(dfb_rec, 128);
    sprintf(filename, "rec_dfb%d.pgm", i);
    mat_pgm_write(filename, dfb_rec);
    mat_decr(dfb_rec, 128);
    mat_delete(dfb_rec);
  }

  /* normalize the subbands */
  for(i = 0; i < ct_levels; i++)
    for(j = 0; j < (1 << dfb_levels[i]); j++)
      mat_div_by(contourlet->high[i][j], norm_high[1+i][dfb_levels[i]][j]);
  mat_div_by(contourlet->low, norm_low[ct_levels]);


  //  mat_pgm_write("rec_low.pgm", contourlet->dwt[0]);

  /* TODO: fix this in libit */
  if(wt_levels)
    wavelet = it_wavelet2D_merge(contourlet->dwt, wt_levels);
  else
    mat_copy(wavelet, contourlet->dwt[0]);

  mat_pgm_write("rec_dwt.pgm", wavelet);

  contourlet->low = it_idwt2D(wavelet, it_wavelet_lifting_97, wt_levels);

  //Contourlet重构
  contourlet_itransform(contourlet, dest);
   IplImage* result_img=save_mat(dest);
 
  //释放变量 
   contourlet_delete(contourlet);

  mat_pgm_write("rec.pgm", dest);
  IplImage* rec=save_mat(dest);
cvSaveImage("rec0.jpg",rec); 
  printf("rate = %f PSNR = %f\n", length * 8. / (w*h), 10*log10(255*255/mat_distance_mse(source, dest, 0)));

  ivec_delete(dfb_levels);
  mat_delete(dest);
  mat_delete(source);
  bvec_delete(buffer);

  return(0);
}
