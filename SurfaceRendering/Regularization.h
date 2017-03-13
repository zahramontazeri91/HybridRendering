#ifndef REGULARIZATION
#define REGULARIZATION
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <iostream>
#include <vector>

using namespace cv;
using namespace std;

vector<Mat> regularization(vector<Mat> morphed_patches, vector < vector<Point> > fixedPoints, int padding );
vector<Mat> regularization(vector < vector<Point> > fixedPoints, int padding);

#endif