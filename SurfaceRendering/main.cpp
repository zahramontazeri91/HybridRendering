#include "opencv2/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
//#include <opencv2/core/eigen.hpp>
#include <iostream>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include "ImageProcessing.h"
#include "Regularization.h"
#include "ControlPoints.h"
#include "Morphing.h"
#include <time.h>

using Eigen::MatrixXd;
using namespace cv;
using namespace std;

int main(int argc, char** argv)
{
	clock_t init, final;
	init = clock();

	Mat grayImage;
	vector<Mat> alignedMasks;
	vector<Mat> alignedMasks_pad;
	vector<Mat> masks;
	vector<Mat> masks_pad;
	vector<Mat> patches;
	vector<Mat> patches_pad;
	vector<Mat> morphed_patches;
	vector<Mat> reg_patches;
	int padding = 30;
	vector < vector<Point> > movingPoints;
	vector < vector<Point> > fixedPoints;

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
	for (int i = 0; i < 30; i++) {
		Mat temp;
		grayImage.copyTo(temp, masks[i]);
		patches.push_back(temp);
		morphed_patches.push_back(Mat::zeros(grayImage.rows, grayImage.cols, CV_32F));
		//imwrite("Patches/patch_" + std::to_string(i) + ".png", temp);
	}
	//imshow("image", patches[12]);

	///*************************Control Point section
	for (int i = 0; i < 30; i++) {
		masks_pad.push_back( ZeroPadding(masks[i], padding));
		alignedMasks_pad.push_back(ZeroPadding(alignedMasks[i], padding));
		patches_pad.push_back(ZeroPadding(patches[i], padding));

		vector<Point> movingPoints_tmp;
		vector<Point> fixedPoints_tmp;
		ControlPoints(masks_pad[i], alignedMasks_pad[i], movingPoints_tmp, fixedPoints_tmp);
		movingPoints.push_back(movingPoints_tmp);
		fixedPoints.push_back(fixedPoints_tmp);
	}

	///************************ morphing section using the control points for each patch
	//for (int i = 0; i < 30; i++) {
	//	Mat temp;
	//	Mat grayImage_pad = ZeroPadding(grayImage, padding);
	//	morphing(temp, grayImage_pad, movingPoints[i], fixedPoints[i]);

	//	///change the size of the results back, before padding
	//	// Setup a rectangle to define the region of interest with width and height
	//	cv::Rect myROI(padding, padding, grayImage.cols, grayImage.rows);
	//	//// Crop the full image to that image contained by the rectangle myROI
	//	cv::Mat croppedTemp = temp(myROI);
	//	croppedTemp.copyTo(morphed_patches[i], alignedMasks[i]);
	//	imwrite("Morphed Patches/morphed_patch_" + std::to_string(i) + ".png", morphed_patches[i]);
	//}

	///************************Regression section for each patch
	vector<Mat> warpsReg;
	vector<Mat> weftsReg;
	regularization(warpsReg, weftsReg);
	/// now let's separate out the patches from regularized yarns
	Mat temp;
	weftsReg[0].copyTo(temp, alignedMasks[0]);
	reg_patches.push_back(temp);
	warpsReg[0].copyTo(temp, alignedMasks[1]);
	reg_patches.push_back(temp);
	weftsReg[3].copyTo(temp, alignedMasks[2]);
	reg_patches.push_back(temp);
	warpsReg[0].copyTo(temp, alignedMasks[3]);
	reg_patches.push_back(temp);
	warpsReg[1].copyTo(temp, alignedMasks[4]);
	reg_patches.push_back(temp);
	weftsReg[2].copyTo(temp, alignedMasks[5]);
	reg_patches.push_back(temp);
	warpsReg[1].copyTo(temp, alignedMasks[6]);
	reg_patches.push_back(temp);
	weftsReg[5].copyTo(temp, alignedMasks[7]);
	reg_patches.push_back(temp);
	warpsReg[2].copyTo(temp, alignedMasks[8]);
	reg_patches.push_back(temp);
	weftsReg[1].copyTo(temp, alignedMasks[9]);
	reg_patches.push_back(temp);
	warpsReg[2].copyTo(temp, alignedMasks[10]);
	reg_patches.push_back(temp);
	weftsReg[4].copyTo(temp, alignedMasks[11]);
	reg_patches.push_back(temp);
	warpsReg[2].copyTo(temp, alignedMasks[12]);
	reg_patches.push_back(temp);
	weftsReg[0].copyTo(temp, alignedMasks[13]);
	reg_patches.push_back(temp);
	warpsReg[3].copyTo(temp, alignedMasks[14]);
	reg_patches.push_back(temp);
	weftsReg[4].copyTo(temp, alignedMasks[15]);
	reg_patches.push_back(temp);
	warpsReg[3].copyTo(temp, alignedMasks[16]);
	reg_patches.push_back(temp);
	warpsReg[4].copyTo(temp, alignedMasks[17]);
	reg_patches.push_back(temp);
	weftsReg[2].copyTo(temp, alignedMasks[18]);
	reg_patches.push_back(temp);
	warpsReg[4].copyTo(temp, alignedMasks[19]);
	reg_patches.push_back(temp);
	weftsReg[5].copyTo(temp, alignedMasks[20]);
	reg_patches.push_back(temp);
	warpsReg[5].copyTo(temp, alignedMasks[21]);
	reg_patches.push_back(temp);
	weftsReg[1].copyTo(temp, alignedMasks[22]);
	reg_patches.push_back(temp);
	warpsReg[5].copyTo(temp, alignedMasks[23]);
	reg_patches.push_back(temp);
	weftsReg[4].copyTo(temp, alignedMasks[24]);
	reg_patches.push_back(temp);
	warpsReg[5].copyTo(temp, alignedMasks[25]);
	reg_patches.push_back(temp);
	weftsReg[0].copyTo(temp, alignedMasks[26]);
	reg_patches.push_back(temp);
	warpsReg[6].copyTo(temp, alignedMasks[27]);
	reg_patches.push_back(temp);
	weftsReg[3].copyTo(temp, alignedMasks[28]);
	reg_patches.push_back(temp);
	warpsReg[6].copyTo(temp, alignedMasks[29]);
	reg_patches.push_back(temp);

	////******************************* Undo_Morphing section
	//we can first do the unmorphing and then masking

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
	final = clock() - init;
	cout << "the program is finished in "<< (double)final / ((double)CLOCKS_PER_SEC);

	waitKey(0); // Wait for a keystroke in the window
	return 0;
}

