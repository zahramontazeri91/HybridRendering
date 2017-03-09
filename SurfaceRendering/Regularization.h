#ifndef REGULARIZATION
#define REGULARIZATION
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <iostream>
#include <vector>

using namespace cv;
using namespace std;

void regularization(vector<Mat>& warpsReg, vector<Mat>& weftsReg);
void regularization(vector<Mat> morphed_patches, vector<Mat>& warpsReg, vector<Mat>& weftsReg);
Mat warpRegression(Mat im, int startCol, int endCol);
Mat weftRegression(Mat im, int startCol, int endCol);

#endif