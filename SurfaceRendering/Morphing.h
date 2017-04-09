#ifndef MORPHING
#define MORPHING
#include "opencv2/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace cv;

void slowLaplace(MatrixXd Dx, MatrixXd Dy, MatrixXd temp, MatrixXd& Lx, MatrixXd& Ly);
void morphing(Mat& img_new, Mat img, vector<Point> movingPoints, vector<Point> fixedPoints);
void fastLaplace(MatrixXd Dx, MatrixXd Dy, MatrixXd id);

#endif

