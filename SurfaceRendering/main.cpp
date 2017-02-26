#include "opencv2/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
//#include <opencv2/core/eigen.hpp>
#include <iostream>
#include "ControlPoints.h"
#include "Morphing.h"
#include <iostream>
#include <Eigen/Dense>
#include <vector>

using Eigen::MatrixXd;
using namespace cv;
using namespace std;

int main(int argc, char** argv)
{
	Mat grayImage;
	vector<Mat> alignedMasks;
	vector<Mat> masks;

	grayImage = imread("height.exr", IMREAD_GRAYSCALE); // Read the height map
	if (!grayImage.data) // Check for invalid input
	{
		cout << "Could not open or find the image" << std::endl;
		return -1;
	}

	///***************************segmentation section
	//segmentation();
	for (int i = 0; i < 30; i++) {
		alignedMasks.push_back( imread("masks/aligned/patch_"+ to_string(i+1)+ ".png", CV_LOAD_IMAGE_GRAYSCALE) );
	}

	for (int i = 0; i < 30; i++) {
		masks.push_back(imread("masks/manually/patch_" + to_string(i + 1) + ".png", CV_LOAD_IMAGE_GRAYSCALE));
	}

	//namedWindow("image", CV_WINDOW_AUTOSIZE);
	//imshow("image", masks[1]);

	///*************************Control Point section
	for (int i = 0; i < 30; i++) {
		vector<Point> movingPoints;
		vector<Point> fixedPoints;
		ControlPoints(masks[i], alignedMasks[i], movingPoints, fixedPoints);
		///************************do morphing section using the points
	}

	//MatrixXd Dx = MatrixXd::Zero(2,2);
	//MatrixXd Dy = MatrixXd::Zero(2,2);
	//MatrixXd temp = MatrixXd::Zero(2,2);
	//MatrixXd Lx(2,2), Ly(2,2);
	//fastLaplace(Dx, Dy, temp, Lx, Ly);

	waitKey(0); // Wait for a keystroke in the window
	return 0;
}