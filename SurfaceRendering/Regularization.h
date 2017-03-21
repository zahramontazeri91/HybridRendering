#ifndef REGULARIZATION
#define REGULARIZATION
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <iostream>
#include <vector>

using namespace cv;
using namespace std;

vector<Mat> regularization(vector<Mat> morphed_patches, int padding );
vector<Mat> regularization(int padding);

#endif