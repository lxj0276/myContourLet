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
#include "load_mat.h"
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