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

mat load_mat(IplImage* img);
