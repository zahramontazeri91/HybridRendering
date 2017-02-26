#ifndef MORPHING
#define MORPHING
#include "opencv2/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace cv;

void fastLaplace(MatrixXd Dx, MatrixXd Dy, MatrixXd temp, MatrixXd& Lx, MatrixXd& Ly);
#endif

