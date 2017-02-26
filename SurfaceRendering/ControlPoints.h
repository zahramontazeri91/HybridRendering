#ifndef _CONTROL_POINTS
#define _CONTROL_POINTS

#include "opencv2/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>
#include <vector>

using namespace cv;
using namespace std;

void ControlPoints(Mat movingImage, Mat fixedImage, vector<Point>& movingPoints, vector<Point>& fixedPoints);

#endif 