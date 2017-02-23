#include "ImageProcessing.h"
#include "opencv2/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using namespace cv;
using namespace std;

void ControlPoints(Mat movingImage, Mat fixedImage, vector<Point>& movingPoints, vector<Point>& fixedPoints)
{
	movingImage = ZeroPadding(movingImage);
	fixedImage = ZeroPadding(fixedImage);

	Mat edge_moving = EdgeDetector(movingImage);
	Mat edge_fixed = EdgeDetector(fixedImage);
	vector<Point> alignCorner;

	//namedWindow("image", CV_WINDOW_AUTOSIZE);
	//imshow("image", edge_moving);

	/// capturing the 4 cornels of aligned mask :
	CornerDetector(fixedImage, alignCorner);
	fixedPoints = alignCorner;
	CornerDetector(movingImage, movingPoints);

	return;
}
