#include"stdafx.h"
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
#include "save_mat.h"
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