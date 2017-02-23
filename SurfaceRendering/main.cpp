#include "opencv2/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>
#include "ControlPoints.h"
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
	vector<vector<Point>> moving_CPs;
	vector<vector<Point>> fixed_CPs;
	vector<Point> movingPoints;
	vector<Point> fixedPoints;

	grayImage = imread("height.exr", IMREAD_GRAYSCALE); // Read the height map
	if (!grayImage.data) // Check for invalid input
	{
		cout << "Could not open or find the image" << std::endl;
		return -1;
	}

	///***************************segmentation section
	//segmentation();
	for (int i = 0; i < 3; i++) {
		alignedMasks.push_back( imread("masks/align_warp_"+ to_string(i+1)+ ".png", CV_LOAD_IMAGE_GRAYSCALE) );
	}
	//namedWindow("image", CV_WINDOW_AUTOSIZE);
	//imshow("image", listOfMasks[1]);

	for (int i = 0; i < 3; i++) {
		masks.push_back(imread("masks/warp_" + to_string(i + 1) + ".png", CV_LOAD_IMAGE_GRAYSCALE));
	}

	///*************************Control Point section
	for (int i = 0; i < 3; i++) {
		ControlPoints(masks[i], alignedMasks[i], movingPoints, fixedPoints);
		moving_CPs.push_back(movingPoints);
		fixed_CPs.push_back(fixedPoints);
	}
	//cout << moving_CPs[0];

	///************************morphing section

	waitKey(0); // Wait for a keystroke in the window
	return 0;
}