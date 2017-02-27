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
#include "ImageProcessing.h"

using Eigen::MatrixXd;
using namespace cv;
using namespace std;

int main(int argc, char** argv)
{
	Mat grayImage;
	vector<Mat> alignedMasks;
	vector<Mat> masks;
	vector<Mat> morphed_patches;
	vector<Mat> patches;
	int padding = 30;

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
	///mask the input to get all the patches 
	int n = grayImage.rows + 2*padding;
	int m = grayImage.cols + 2*padding;
	for (int i = 0; i < 30; i++) {
		Mat temp;
		grayImage.copyTo(temp, masks[i]);
		patches.push_back(temp);
		morphed_patches.push_back(Mat::zeros(n, m, CV_32F));
	}
	//imshow("image", patches[12]);

	///*************************Control Point section
	Mat grayImage_pad = ZeroPadding(grayImage, padding);
	for (int i = 0; i < 30; i++) {
		vector<Point> movingPoints;
		vector<Point> fixedPoints;
		masks[i] = ZeroPadding(masks[i], padding);
		alignedMasks[i] = ZeroPadding(alignedMasks[i], padding);
		patches[i] = ZeroPadding(patches[i], padding);

		ControlPoints(masks[i], alignedMasks[i], movingPoints, fixedPoints);
		///************************do morphing section using the control points for each patch
		morphing(morphed_patches[i], grayImage_pad, movingPoints, fixedPoints);
		Mat temp;
		morphed_patches[i].copyTo(temp, alignedMasks[i]);
		imwrite("Morphed Patches/morphed_patch_" + std::to_string(i) + ".png", temp);
	}

	/// To check FastLaplace() is working
	//MatrixXd Dx(4,5);
	//Dx << 8.9, 8.9, 8.9, 8.9, 8.9,
	//	8.4, 0, 0, 0, 9.2,
	//	7.2, 0, 10.2, 0, 9.4,
	//	6.1, 6.8, 7.7, 8.7, 6.1;
	//MatrixXd Dy(4,5);
	//Dy << 8.9, 8.9, 8.9, 8.9, 8.9,
	//	8.4, 0, 0, 0, 9.2,
	//	7.2, 0, 0, 0, 9.4,
	//	6.1, 6.8, 7.7, 8.7, 6.1;
	//MatrixXd temp(4,5);
	//temp << 8.9, 8.9, 8.9, 8.9, 8.9,
	//	8.4, 0, 0, 0, 9.2,
	//	7.2, 0, 1, 0, 9.4,
	//	6.1, 6.8, 7.7, 8.7, 6.1;
	//MatrixXd Lx(4,5), Ly(4,5);

	waitKey(0); // Wait for a keystroke in the window
	return 0;
}